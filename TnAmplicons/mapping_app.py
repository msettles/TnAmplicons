
import os
import sys
import traceback
import time
import signal
import re
from collections import Counter
from itertools import groupby

from subprocess import Popen
from subprocess import PIPE

from TnAmplicons import misc


def sp_bowtie2_index(ref, overwrite=False):
    if os.path.isfile(ref):
        if os.path.isfile(ref + '.rev.2.bt2') and not overwrite:
            sys.stderr.write('Found existing bowtie2 index for %s\n' % ref)
            return 0
        else:
            FNULL = open(os.devnull, 'w')
            call = 'bowtie2-build'
            call = call + ' ' + ref + ' ' + ref
            p = Popen(['bowtie2-build', ref, ref],
                      stdout=FNULL,
                      stderr=FNULL,
                      bufsize=-1,
                      preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
            p.communicate()
            if p.returncode:
                sys.stderr.write('Something in bowtie2-build went wrong\n')
                raise
            # system call, check for return
            sys.stderr.write('Successfully indexed %s\n' % ref)
            return 0
    else:
        sys.stderr.write("%s Reference file not found\n" % ref)
        return 1
    sys.stderr.write('Something in bowtie2-build went wrong\n')
    raise


# Template
# bowtie2 -x caplanStuff -U <(zcat ../../CaplanShit/00-RawData/Sample_GCCAAT/GCCAAT_R1.fastq.gz| sed 's, ,_,g')
def sp_bowtie2_screen(pe1, pe2, se, ref, overwrite=False, sensitivity=0, procs=1, minins=0, maxins=1000):
    # build the call,
    # each file must first go through awk to replace spaces with a parsable character
    if sp_bowtie2_index(ref, overwrite) != 0:
        sys.exit(1)

    sensitivity_switch = ['--very-fast', '--fast', '--sensitive', '--very-sensitive']
    call = 'bowtie2 -I ' + str(minins) + ' -X ' + str(maxins) + ' ' + sensitivity_switch[sensitivity] + ' -p ' + str(procs) + ' -x ' + ref
    if ((pe1 is not None) and (pe2 is not None) and (len(pe1) == len(pe2))):
        pe1_gz = "gunzip -c"
        pe2_gz = "gunzip -c"
        pe1_gz_true = False
        pe2_gz_true = False
        pe1_ngz = "cat"
        pe2_ngz = "cat"
        pe1_ngz_true = False
        pe2_ngz_true = False
        for pe_read in pe1:
            if pe_read.split(".")[-1] == "gz":
                pe1_gz = pe1_gz + " " + pe_read
                pe1_gz_true = True
            else:
                pe1_ngz = pe1_ngz + " " + pe_read
                pe1_ngz_true = True
        for pe_read in pe2:
            if pe_read.split(".")[-1] == "gz":
                pe2_gz = pe2_gz + " " + pe_read
                pe2_gz_true = True
            else:
                pe2_ngz = pe2_ngz + " " + pe_read
                pe2_ngz_true = True
        if pe1_gz_true is True:
            call = call + " -1 <(" + pe1_gz + ")"
        if pe2_gz_true is True:
            call = call + " -2 <(" + pe2_gz + ")"
        if pe1_ngz_true is True:
            call = call + " -1 <(" + pe1_ngz + ")"
        if pe2_ngz_true is True:
            call = call + " -2 <(" + pe2_ngz + ")"
    if (se is not None):
        se_gz = "gunzip -c"
        se_gz_true = False
        se_ngz = "cat"
        se_ngz_true = False
        for se_read in se:
            if se_read.split(".")[-1] == "gz":
                se_gz = se_gz + " " + se_read
                se_gz_true = True
            else:
                se_ngz = se_ngz + " " + se_read
                se_ngz_true = True
        if se_gz_true is True:
            call = call + " -U <(" + se_gz + ")"
        if se_ngz_true is True:
            call = call + " -U <(" + se_gz + ")"
    sys.stdout.write(call + '\n')
    p = Popen(call,
              stdout=PIPE,
              stderr=None,
              bufsize=-1,
              shell=True,
              executable='/bin/bash',
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'a': 't', 'A': 'T', 'c': 'g', 'C': 'G', 'g': 'c', 'G': 'C', 't': 'a', 'T': 'A', 'n': 'n', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])


def reverse(s):
    """
    given a sequence return the reverse
    """
    letters = list(s)
    return ''.join(letters[::-1])


class cigarString:
    """ Class to parse and handle sam formatted cigar strings """
    pattern = re.compile('([MIDNSHPX=])')

    def __init__(self, cigar):
        values = self.pattern.split(cigar)[:-1]
        self.paired = (values[n:n+2] for n in xrange(0, len(values), 2))  # pair values by twos

    def getAlignmentLength(self):
        g = 0
        for pair in self.paired:
            l = int(pair[0])
            t = pair[1]
            if t == 'M':
                g += l
            elif t == 'I':
                pass
            elif t == 'D':
                g += l
            elif t == 'N':
                pass
            elif t == 'S':
                pass
            elif t == 'H':
                pass
            elif t == 'P':
                pass
            elif t == 'X':
                pass
            elif t == '=':
                pass
            else:
                sys.stderr.write("encountered unhandled CIGAR character %s\n" % t)
                pass
        return g


def getUniqueKey(r1_line, r2_line=None):

    if r2_line is not None:
        if r1_line[0] != r2_line[0]:
            sys.stderr.write("Something went wrong building unique key Read 1 and Read 2 names do not match\n")
            sys.exit(1)
        if r1_line[2] != r2_line[2]:
            sys.stderr.write("Something went wrong building unique key Read 1 and Read 2 chrom do not match\n")
            sys.exit(1)
        cigar2 = cigarString(r2_line[5])

    # J00113:93:H5GNNBBXX:6:1101:28605:1261*NCGGAACA|Linker*Sol_AP1-1*1|CAAGTAAAA|Primer*Sol_Mar_4b*1*Adapter_2.2_Bar_B*0
    id1 = r1_line[0].split("|")  # 0) originial name + bc 1) Linker (3) 2) Molecular index (1) 3) Interior primers (5)
    flag1 = int(r1_line[1])
    cigar1 = cigarString(r1_line[5])
    # barcode, Molecular Index, P5 Primer, P7 Primer
    # uniquekey = ','.join([id1[0].split('*')[1], id1[2], id1[3].split('*')[1], id1[1].split('*')[1]])
    uniquekey = ','.join([id1[0].split('*')[1], id1[2]])  # does not include primers

    if (flag1 & 0x10):  # if RevComp
        if r2_line is not None:
            uniquekey = '_*_'.join([uniquekey, 'R', r1_line[2], str(int(r1_line[3])+cigar1.getAlignmentLength()-2), r2_line[3]])
        else:
            uniquekey = '_*_'.join([uniquekey, 'R', r1_line[2], str(int(r1_line[3])+cigar1.getAlignmentLength()-2), r1_line[3]])
    else:
        if r2_line is not None:
            uniquekey = '_*_'.join([uniquekey, 'F', r1_line[2], r1_line[3], str(int(r2_line[3])+cigar2.getAlignmentLength()-1)])
        else:
            uniquekey = '_*_'.join([uniquekey, 'F', r1_line[2], r1_line[3], str(int(r1_line[3])+cigar1.getAlignmentLength()-1)])
    return uniquekey


def getISites(reference, insertionSite='TA'):
    # given a reference file and insertion site, identify all possible insertion site locations on both strands
    sites = {}
    ref = open(reference)
    faiter = (list(x[1]) for x in groupby(ref, lambda line: line[0] == ">"))
    for g in faiter:
        # drop the ">"
        header = g[0][1:].strip().split(' ')[0]
        # join all sequence lines to one.
        seq = str("".join(s.strip() for s in faiter.next())).upper()
        sites[header] = [str(m.start() + 1) for m in re.finditer(insertionSite, seq)]
    # sites keyed by 'contig' ID then vector of sites (RC sites are +1bp away)
    return sites


class mappingApp:

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, reference, insertionSite, overwrite, sensitivity, output_prefix, minins, maxins, procs, mapq, dedup_reads=True, uncompressed=False, verbose=True, debug=False):
        """
            screen reads against a reference fasta file
        """
        try:
            mapped_pairs_count = 0
            mapped_pairs_lowqual = 0
            unmapped_pairs_count = 0
            mapped_singles_count = 0
            mapped_singles_lowqual = 0
            unmapped_singles_count = 0
            secondary_alignment = 0
            missing_IS = 0  # missing the insersion site
            duplicate_count = 0  # duplcicates (OF mapped!)
            count = Counter()

            # Set up output
            misc.make_sure_path_exists(os.path.dirname(output_prefix))
            run_out = {}
            run_out["mappedsam"] = open(output_prefix + '.sam', 'w')

            read_profile = open(output_prefix + '.TnRead_profile.txt', 'w')

            # 0x1 template having multiple segments in sequencing
            # 0x2 each segment properly aligned according to the aligner
            # 0x4 segment unmapped
            # 0x8 next segment in the template unmapped
            # 0x10 SEQ being reverse complemented
            # 0x20 SEQ of the next segment in the template being reversed
            # 0x40 the first segment in the template
            # 0x80 the last segment in the template
            # 0x100 secondary alignment
            # 0x200 not passing quality controls
            # 0x400 PCR or optical duplicate
            PE1 = {}
            PE2 = {}

            # Read in reference and identify all potential insertion sites
            site_counts_F = {}
            site_counts_R = {}
            sites = getISites(reference, insertionSite)
            for contig in sites.keys():
                sys.stdout.write("Found %s insertion sites for contig %s\n" % (len(sites[contig]), contig))
                site_counts_F[contig] = Counter()
                site_counts_R[contig] = Counter()

            i = 0
            lasttime = time.time()

            for line in sp_bowtie2_screen(fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, procs, minins, maxins):
                if i % 10000 == 0 and i > 0:
                    read_profile.write("%s\t%s\t%s\t%s\t%s\n" % (i, (mapped_pairs_count+mapped_singles_count), (mapped_pairs_count+mapped_singles_count+duplicate_count), sum(len(site_counts_F[c]) for c in site_counts_F.keys()), sum(len(site_counts_R[c]) for c in site_counts_R.keys())))
                if i % 100000 == 0 and i > 0 and verbose:
                        sys.stderr.write("Processed: %s, PE in ref: %s, SE in ref: %s in %s minutes\n" % (i, mapped_pairs_count, mapped_singles_count, round((time.time()-lasttime)/(60), 2)))
                if line[0] == "@":  # header line
                    # write out to sam
                    run_out["mappedsam"].write(line)
                else:
                    i += 1

                    line2 = line.strip().split()
                    flag = int(line2[1])
                    # Secondary alignment
                    if (flag & 0x100):
                        secondary_alignment += 1
                        continue

                    # check for insersion site before any further processing
                    if (not (flag & 0x1)) or (flag & 0x40):
                        if (flag & 0x10):
                            IS = reverseComplement(line2[9])[0:2] == insertionSite
                        else:
                            IS = line2[9][0:2] == insertionSite
                        if not IS:
                            missing_IS += 1
                            continue

                    # Handle SE:
                    # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1 if mapped
                    if not (flag & 0x1):  # SE READ
                        if not (flag & 0x4):  # MAPPED
                            if int(line2[4]) >= mapq:  # check mapq
                                key = getUniqueKey(line2)
                                if dedup_reads:
                                    count[key] += 1
                                    if count[key] == 1:
                                        if key.split('_*_')[1] == 'F':
                                            site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                        else:
                                            site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                        mapped_singles_count += 1
                                        run_out["mappedsam"].write('\t'.join(line2) + '\n')
                                    else:
                                        duplicate_count += 1
                                else:
                                    if key.split('_*_')[1] == 'F':
                                        site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                    else:
                                        site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                    mapped_singles_count += 1
                                    run_out["mappedsam"].write('\t'.join(line2) + '\n')
                            else:  # MAPPED BUT LOW QUAL
                                mapped_singles_lowqual += 1
                        else:  # UNMAPPED
                            unmapped_singles_count += 1
                        continue
                    # Handle PE:
                    # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
                    if (flag & 0x1):  # PE READ
                        if ((not (flag & 0x4) and not (flag & 0x8)) and (flag & 0x2)):  # both pair mapped and concordant
                            if int(line2[4]) >= mapq:  # check mapq
                                if (flag & 0x40):  # is this PE1 (first segment in template)
                                    # PE1 read, check that PE2 is in dict and write out
                                    ID = line2[0]
                                    if ID in PE2:
                                        if (flag & 0x10):  # reverse complement
                                            line2[1] = str(flag - 0x1 - 0x2 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                        else:  # forward complements
                                            line2[1] = str(flag - 0x1 - 0x2 - 0x20 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                        key = getUniqueKey(line2, PE2[ID])
                                        if dedup_reads:
                                            count[key] += 1
                                            if count[key] == 1:
                                                if key.split('_*_')[1] == 'F':
                                                    site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                                else:
                                                    site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                                mapped_pairs_count += 1
                                                run_out["mappedsam"].write('\t'.join(line2) + '\n')
                                            else:
                                                duplicate_count += 1
                                        else:
                                            if key.split('_*_')[1] == 'F':
                                                site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                            else:
                                                site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                            mapped_pairs_count += 1
                                            run_out["mappedsam"].write('\t'.join(line2) + '\n')

                                        del PE2[ID]
                                    else:
                                        PE1[ID] = line2
                                elif (flag & 0x80):  # is this PE2 (last segment in template)
                                    # PE2 read, check that PE1 is in dict and write out
                                    ID = line2[0]
                                    if ID in PE1:
                                        if (int(PE1[ID][1]) & 0x10):  # reverse complement
                                            PE1[ID][1] = str(int(PE1[ID][1]) - 0x1 - 0x2 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                        else:  # forward complements
                                            PE1[ID][1] = str(int(PE1[ID][1]) - 0x1 - 0x2 - 0x20 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                        key = getUniqueKey(PE1[ID], line2)
                                        if dedup_reads:
                                            count[key] += 1
                                            if count[key] == 1:
                                                if key.split('_*_')[1] == 'F':
                                                    site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                                else:
                                                    site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                                mapped_pairs_count += 1
                                                run_out["mappedsam"].write('\t'.join(PE1[ID]) + '\n')
                                            else:
                                                duplicate_count += 1
                                        else:
                                            if key.split('_*_')[1] == 'F':
                                                site_counts_F[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                            else:
                                                site_counts_R[key.split('_*_')[2]][key.split('_*_')[3]] += 1
                                            mapped_pairs_count += 1
                                            run_out["mappedsam"].write('\t'.join(PE1[ID]) + '\n')

                                        del PE1[ID]
                                    else:
                                        PE2[ID] = line2
                            else:
                                mapped_pairs_lowqual += 1
                        else:  # an 'unmapped' pair (both pairs unmapped, one of pair unmapped, or both mapped but discordant)
                            unmapped_pairs_count += 1

            sys.stderr.write("Processed: %s, PE in ref: %s, SE in ref: %s in %s minutes\n" % (i, mapped_pairs_count, mapped_singles_count, round((time.time()-lasttime)/(60), 2)))

            read_profile.write("%s\t%s\t%s\t%s\t%s\n" % (i, (mapped_pairs_count+mapped_singles_count), (mapped_pairs_count+mapped_singles_count+duplicate_count), sum(len(site_counts_F[c]) for c in site_counts_F.keys()), sum(len(site_counts_R[c]) for c in site_counts_R.keys())))
            read_profile.close()

            site_output = open(output_prefix + '.sites', 'w')
            wig_output = open(output_prefix + '.wig', 'w')
            wig_output.write("# output:%s fastq1:%s fastq2:%s fastq3:%s\n" % (output_prefix, fastq_file1, fastq_file2, fastq_file3))
            for contig in sites.keys():
                wig_output.write("variableStep chrom=%s\n" % contig)
                for site in sites[contig]:
                    site_output.write("%s\t%s\t%s\t%s\n" % (contig, site, site_counts_F[contig][site], site_counts_R[contig][site]))
                    wig_output.write("%s %s\n" % (site, str(int(site_counts_F[contig][site]) + int(site_counts_R[contig][site]))))
            site_output.close()
            wig_output.close()

            if debug:
                data = open("reverse.txt", 'w')
                for contig in sites.keys():
                    countF = 0
                    sitesF = 0
                    countR = 0
                    sitesR =0
                    for site in site_counts_F[contig].keys():
                        countF += site_counts_F[contig][site]
                        sitesF += 1
                    for site in site_counts_R[contig].keys():
                        countR += site_counts_R[contig][site]
                        data.write("%s\t%s\n" % (site, site_counts_R[contig][site]))
                        sitesR += 1
                    print countF
                    print sitesF
                    print countR
                    print sitesR
                data.close()

            sys.stdout.write("total records: %s\n" % i)
            sys.stdout.write("secondary alignments: %s\n" % secondary_alignment)
            sys.stdout.write("pairs: %s\n" % (mapped_pairs_count + mapped_pairs_lowqual + unmapped_pairs_count))
            sys.stdout.write("\tmapped pairs: %s\n" % mapped_pairs_count)
            sys.stdout.write("\tlowqual pairs: %s\n" % mapped_pairs_lowqual)
            sys.stdout.write("\tunmapped pairs: %s\n" % unmapped_pairs_count)
            sys.stdout.write("singles: %s\n" % (mapped_singles_count + mapped_singles_lowqual + unmapped_singles_count))
            sys.stdout.write("\tmapped singles: %s\n" % mapped_singles_count)
            sys.stdout.write("\tlowqual singles: %s\n" % mapped_singles_lowqual)
            sys.stdout.write("\tunmapped singles: %s\n" % unmapped_singles_count)
            sys.stdout.write("missing insertion site: %s\n" % missing_IS)
            sys.stdout.write("PCR duplicate count: %s\n" % duplicate_count)

            self.clean()
            return 0

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1
        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")
        try:
            for key in self.run_out:
                self.run_out[key].close()
            self.read_profile.close()
            self.site_output.close()
            pass
        except:
            pass

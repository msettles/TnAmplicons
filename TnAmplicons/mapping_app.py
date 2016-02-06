
import os
import sys
import traceback
import time
import signal

from subprocess import Popen
from subprocess import PIPE

from dbcAmplicons import misc


def sp_bowtie2_index(ref, overwrite=False):
    if os.path.isfile(ref):
        if os.path.isfile(ref + '.rev.2.bt2') and not overwrite:
            print 'Found existing bowtie2 index for %s' % ref
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
                print 'Something in bowtie2-build went wrong'
                raise
            # system call, check for return
            print 'Successfully indexed %s' % ref
            return 0
    else:
        print "%s Reference file not found" % ref
        return 1
    print 'Something in bowtie2-build went wrong'
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
    print call
    FNULL = open(os.devnull, 'w')
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
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
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


class mappingApp:

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, output_prefix, minins, maxins, procs, mapq, uncompressed=False, verbose=True, debug=False):
        """
            screen reads against a reference fasta file
        """
        self.verbose = verbose
        try:
            mapped_pairs_count = 0
            mapped_pairs_lowqual = 0
            unmapped_pairs_count = 0
            mapped_singles_count = 0
            mapped_singles_lowqual = 0
            unmapped_singles_count = 0
            secondary_alignment = 0

            # Set up output
            misc.make_sure_path_exists(os.path.dirname(output_prefix))
            self.run_out = {}
            self.run_out["mappedsam"] = open(output_prefix + '.sam', 'w')

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

            i = 0
            for line in sp_bowtie2_screen(fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, procs, minins, maxins):
                if i % 100000 == 0 and i > 0:
                    if self.verbose:
                        print "Processed: %s, PE in ref: %s, SE in ref: %s" % (i, mapped_pairs_count, mapped_singles_count)
                if line[0] == "@":  # header line
                    # write out to sam
                    self.run_out["mappedsam"].write(line)
                else:
                    i += 1

                    line2 = line.strip().split()
                    flag = int(line2[1])
                    # Secondary alignment
                    if (flag & 0x100):
                        secondary_alignment += 1
                        continue

                    # Handle SE:
                    # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1 if mapped
                    if not (flag & 0x1):  # SE READ
                        if not (flag & 0x4):  # MAPPED
                            if int(line2[4] >= mapq):  # check mapq
                                if (flag & 0x10):  # reverse complement
                                    mapped_singles_count += 1
                                    self.run_out["mappedsam"].write('\t'.join(line2) + '\n')
                                else:  # forward complement
                                    mapped_singles_count += 1
                                    self.run_out["mappedsam"].write('\t'.join(line2) + '\n')
                            else:  # MAPPED BUT LOW QUAL
                                mapped_singles_lowqual += 1
                        else:  # UNMAPPED
                            unmapped_singles_count += 1
                        continue
                    # Handle PE:
                    # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
                    if (flag & 0x1):  # PE READ
                        if ((not (flag & 0x4) and not (flag & 0x8)) and (flag & 0x2)):  # both pair mapped and concordant
                            if int(line2[4] >= mapq):  # check mapq
                                if (flag & 0x40):  # is this PE1 (first segment in template)
                                    # PE1 read, check that PE2 is in dict and write out
                                    ID = line2[0]
                                    if (flag & 0x10):  # reverse complement
                                        line2[1] = str(flag - 0x1 - 0x2 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                    else:  # forward complements
                                        line2[1] = str(flag - 0x1 - 0x2 - 0x20 - 0x40)  # modify read1 flag (remove read2 assoc flags)
                                    if ID in PE2:
                                        mapped_pairs_count += 1
                                        self.run_out["mappedsam"].write('\t'.join(line2) + '\n')
                                        del PE2[ID]
                                    else:
                                        PE1[ID] = line2
                                elif (flag & 0x80):  # is this PE2 (last segment in template)
                                    # PE2 read, check that PE1 is in dict and write out
                                    ID = line2[0]
                                    if (flag & 0x10):  # reverse complement
                                        line2[1] = str(flag - 0x1 - 0x2 - 0x40)  # modify read1 flag to be revcomp (remove read2 assoc flags)
                                    else:  # forward complements
                                        line2[1] = str(flag - 0x1 - 0x2 - 0x20 - 0x40)  # modify read1 flasg to be forcomp (remove read2 assoc flags)
                                    if ID in PE1:
                                        mapped_pairs_count += 1
                                        self.run_out["mappedsam"].write('\t'.join(PE1[ID]) + '\n')
                                        del PE1[ID]
                                    else:
                                        PE2[ID] = line2
                            else:
                                mapped_pairs_lowqual += 1
                        else:  # an 'unmapped' pair (both pairs unmapped, one of pair unmapped, or both mapped but discordant)
                            unmapped_pairs_count += 1

            print "Records processed: %s, PE in ref: %s, SE in ref: %s" % (i, mapped_pairs_count, mapped_singles_count)

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
            self.run_out.close()
            pass
        except:
            pass

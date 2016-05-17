# sequenceReads.py
# sequenceReads.py stores and processes individual DNA sequence reads.

import sys
from TnAmplicons import misc

try:
    from TnAmplicons import editdist
    editdist_loaded = True
except ImportError:
    sys.stderr.write("Warning: editdist library not loaded, Insertion/Deletion detection in barcodes and primers will not be performed\n")
    editdist_loaded = False

try:
    from TnAmplicons import trim
    trim_loaded = True
except ImportError:
    sys.stderr.write("Warning: trim library not loaded, trimming using python\n")
    trim_loaded = False


def primerDist(primer_l, read, dedup_float, max_diff, end_match):
    # ------------------- calculate distance between primer sequence and first part of read ------------------------------
    """
    gets edit distance between primer and read sequence
    """
    pr = None
    prMismatch = max_diff+1
    prStartPosition = 0
    prEndPosition = 0

    if editdist_loaded:
        pr_i, prMismatch, prStartPosition, prEndPosition = editdist.bounded_distance_list(primer_l, read, dedup_float, max_diff, end_match)

    else:
        for i, key in enumerate(primer_l):
            read = read[0:len(key)]
            if end_match == 0:
                dist = sum(map(lambda x: x[0] != x[1], zip(read, key)))
            else:
                dist = sum(map(lambda x: x[0] != x[1], zip(read[:-end_match], key[:-end_match])))
                for i in range(len(key)-end_match, len(key)):
                    if read[i] != key[i]:
                        dist += 100
            if dist < prMismatch:
                pr_i = i
                prMismatch = dist
                prEndPosition = len(key)

    if pr_i < 0 or prMismatch > max_diff:
        pr = None
        prStartPosition = 0
        prEndPosition = 0
    else:
        pr = primer_l[pr_i]

    return (pr, prMismatch, prStartPosition, prEndPosition)


# ---------------- Class for 2 read sequence data processed with TnAmplicons preprocess ----------------
class TwoSequenceReadSet:
    """
    Class to hold one Illumina two read set.
    Class processes a read by defining sample and project ids. Finally class returns a paired read set for output.
    """
    def __init__(self, name_1, read_1, qual_1, name_2, read_2, qual_2):
        """
        Initialize a TwoSequenceReadSet with names, two read sequences and cooresponding quality sequence.
        Barcode and primer sequence is inferred by their placement in the read names.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            # Casava processed, check for barcode in the second string
            self.primer_string1 = None
            self.primer_string2 = None
            self.primer = None

            self.read_1 = read_1
            self.qual_1 = qual_1
            self.read_2 = read_2
            self.qual_2 = qual_2
            self.trim_left = len(read_1)
            self.trim_right = len(read_2)
            self.goodRead = False
        except IndexError:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] Read names are not formatted in the expected manner, CASAVA 1.8 fastq output\n')
            raise
        except:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] Unknown error occured initiating read\n')
            raise

    def checkLinker(self, prTable, dedup_float, max_diff, endmatch):

        # Primer One Matching
        linker, linkerMismatch, linkerStartPosition, linkerEndPosition = primerDist(prTable.getLinkersequences(), self.read_2, 0, max_diff, endmatch)

        self.linker = [prTable.getLinkerID(linker), linkerMismatch, linkerStartPosition, linkerEndPosition]
        self.goodRead = self.linker[0] is not None
        if self.goodRead:
            self.read_2 = self.read_2[linkerEndPosition:]
            self.qual_2 = self.qual_2[linkerEndPosition:]
            return 1
        else:
            return 0

    def assignPrimer(self, prTable, dedup_float, max_diff, endmatch):
        """
        Given a primerTable object, the maximum number of allowed difference (mismatch, insertion, deletions) and the
        required number of end match bases (final endmatch bases must match) assign a primer pair ID from the read
        sequences.
        """
        # Primer One Matching
        pr1, pr1Mismatch, pr1StartPosition, pr1EndPosition = primerDist(prTable.getP5sequences(), self.read_1, 0, max_diff, endmatch)

        # Primer Two Matching
        pr2, pr2Mismatch, pr2StartPosition, pr2EndPosition = primerDist(prTable.getP7sequences(), self.read_2, 9, max_diff, endmatch)

        self.MI2 = self.read_2[0:pr2StartPosition]

        # Primer Pair Matching
        combined_pr = prTable.getMatch(pr1, pr2)
        self.primer = [combined_pr[0], combined_pr[1], pr1Mismatch, pr1StartPosition, pr1EndPosition, combined_pr[2], pr2Mismatch, pr2StartPosition, pr2EndPosition]
        self.goodRead = self.goodRead and self.primer[0] is not None
        if self.goodRead:
            return 1
        else:
            return 0

    def getPrimer(self):
        """
        Return the reads primer pair ID
        """
        return self.primer[0]

    def trimRead(self, minQ, minL):
        """
        Trim the read by a minQ score
        """
        if (trim_loaded):
            trim_points = trim.trim(self.qual_1, self.qual_2, minQ)
            self.trim_left = trim_points["left_trim"]
            self.trim_right = trim_points["right_trim"]
            if (self.trim_left < minL or self.trim_right < minL):
                self.goodRead = False

    def getFastq(self):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer[0] is not None:
            read1_name = "%s*%s|Linker*%s*%s|%s|%s*%s*%s*%s*%s 1:N:0:%s" % (self.name, self.barcode, self.linker[0], self.linker[1], self.MI2, self.primer[0], self.primer[1], self.primer[2], self.primer[5], self.primer[6], self.barcode)
            read2_name = "%s*%s|Linker*%s*%s|%s|%s*%s*%s*%s*%s 2:N:0:%s" % (self.name, self.barcode, self.linker[0], self.linker[1], self.MI2, self.primer[0], self.primer[1], self.primer[2], self.primer[5], self.primer[6], self.barcode)
            r1 = '\n'.join([read1_name, self.read_1[(self.primer[4] - 2):self.trim_left], '+', self.qual_1[(self.primer[4] - 2):self.trim_left]])
            r2 = '\n'.join([read2_name, self.read_2[self.primer[8]:self.trim_right], '+', self.qual_2[self.primer[8]:self.trim_right]])
        else:
            read1_name = "%s 1:N:0:%s" % (self.name, self.barcode)
            read2_name = "%s 2:N:0:%s" % (self.name, self.barcode)
            r1 = '\n'.join([read1_name, self.read_1, '+', self.qual_1])
            r2 = '\n'.join([read2_name, self.read_2, '+', self.qual_2])
        return [r1, r2]


# ---------------- Class for 2 read sequence data processed with TnAmplicons preprocess ----------------
class OneSequenceReadSet:
    """
    Class to hold a one Illumina read set, assumes the paired reads produced by TnAmplicons preprocess have been merged
    """
    def __init__(self, name_1, read_1, qual_1):
        """
        Initialize a OneSequenceReadSet with name, one read sequences and cooresponding quality sequence.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        self.goodRead = False
        self.primer = None
        # parse the read name for the sample and primer ids
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            if (len(split_name) == 4):
                self.primer = split_name[1].split(":")[4]
        except IndexError:
            sys.stderr.write('ERROR:[OneSequenceReadSet] Read names are not formatted in the expected manner\n')
            raise
        except:
            sys.stderr.write('ERROR:[OneSequenceReadSet] Unknown error occured initiating read\n')
            raise
        self.read_1 = read_1
        self.qual_1 = qual_1

    def getFastqSRA(self):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        read1_name = "%s 1:N:0:%s" % (self.name, self.barcode)
        r1 = '\n'.join([read1_name, self.read_1, '+', self.qual_1])
        return [r1]

    def getFastq(self):
        """
        Create four line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        if self.primer is not None:
            read1_name = "%s 1:N:0:%s:%s" % (self.name, self.barcode, self.primer)
        else:
            read1_name = "%s 1:N:0:%s" % (self.name, self.barcode)
        r1 = '\n'.join([read1_name, self.read_1, '+', self.qual_1])
        return [r1]

    def getFasta(self):
        """
        Create two line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer is not None:
            read1_name = "%s|%s:%s:%i" % (name, self.barcode, self.primer, len(self.read_1))
        else:
            read1_name = "%s|%s:%i" % (name, self.barcode, len(self.read_1))
        r1 = '\n'.join([read1_name, self.read_1])
        return [r1]

# primers.py
#
# primer lookup file should look like, where Read is P5 or R1 or READ1 and P7 or R2 or READ2,
# the '#' character represents a comments and will be ignored

# #Read    Primer_ID   Sequence
# P5  27F_YM1 GTAGAGTTTGATCCTGGCTCAG
# P5  27F_YM2 CGTAGAGTTTGATCATGGCTCAG
# P5  27F_YM3 ACGTAGAGTTTGATTCTGGCTCAG
# P5  27F_YM4 TACGTAGAGTTTGATTATGGCTCAG
# P5  27F_Bif GTACGTAGGGTTCGATTCTGGCTCAG
# P5  27F_Bor CGTACGTAGAGTTTGATCCTGGCTTAG

"""
primer.py parses and stores primer information associated with a double barcoded illumina amplicon project
"""
import sys
from TnAmplicons import misc


# TODO: This class needs more work
class primerTable:
    # ---------------- primer class ----------------
    """
    Class to read in and hold amplicon pcr primer table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, primerfile):
        """
        Initialize a new primerTable object with the file primer table, parses and stores the primer information
        """
        self.P5sequences = []
        self.P5id = {}
        self.P7sequences = []
        self.P7id = {}
        self.Linkersequences = []
        self.Linkerid = {}
        # TODO: add in check for presense of both P5 and P7 in pair
        try:
            prfile = open(primerfile, 'r')
        except IOError:
            sys.stderr.write('ERROR:[Primers] Error cannot open\n', primerfile)
            raise
        f = prfile.readlines()
        line = 0
        for row in f:
            line += 1
            if row[0] == "#" or row[0] == "\n":  # comment or blank line
                continue
            row = row.rstrip()
            try:
                READ, ID, SEQ = row.split('\t')[0:4]
            except ValueError as e:
                sys.stderr.write("ERROR:[Primers] Error reading line %s of primer file: %s\n" % (str(line), str(e)))
                raise
            except:
                sys.stderr.write("ERROR:[Primers] Unexpected error on line %s of the primers file: %s\n" % (line, sys.exc_info()[0]))
                raise
            pseqs = misc.expand_iupac(SEQ.upper())
            if READ in ["P5", "R1", "READ1", "F", "FORWARD"]:
                for pseq in pseqs:
                    if pseq in self.P5sequences:
                        sys.stderr.write("ERROR:[Primers] sequence %s has already been seen within the file\n" % (pseq))
                        raise
                    else:
                        self.P5sequences.extend([pseq])
                        self.P5id[pseq] = ID
            if READ in ["P7", "R2", "READ2", "R", "REVERSE"]:
                for pseq in pseqs:
                    if pseq in self.P7sequences:
                        sys.stderr.write("ERROR:[Primers] sequence %s has already been seen within the file\n" % (pseq))
                        raise
                    else:
                        self.P7sequences.extend([pseq])
                        self.P7id[pseq] = ID
            if READ in ["Linker", "Linkers"]:
                for pseq in pseqs:
                    if pseq in self.Linkersequences:
                        sys.stderr.write("ERROR:[Primers] sequence %s has already been seen within the file\n" % (pseq))
                        raise
                    else:
                        self.Linkersequences.extend([pseq])
                        self.Linkerid[pseq] = ID

        prfile.close()

    def getP5sequences(self):
        """
        Return the list of P5 sequence
        """
        return self.P5sequences

    def getP7sequences(self):
        """
        Return the list of P7 sequences
        """
        return self.P7sequences

    def getMatch(self, seq1, seq2):
        if seq1 in self.P5id.keys():
            id1 = self.P5id[seq1]
        else:
            id1 = None
        if seq2 in self.P7id.keys():
            id2 = self.P7id[seq2]
        else:
            id2 = None
        # at least one primer not id, TODO: allow for the situation of only 1 primer
        if id1 is None or id2 is None:
            return [None, id1, id2]
        # Catch all
        return ["Primer", id1, id2]

    def getLinkersequences(self):
        """
        Return the list of P5 sequence
        """
        return self.Linkersequences

    def getLinkerID(self, seq1):
        if seq1 in self.Linkerid.keys():
            id1 = self.Linkerid[seq1]
        else:
            id1 = None
        return id1

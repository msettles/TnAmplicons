# preprocess_app.py
#
import sys
import traceback
import time
from TnAmplicons import primerTable

from TnAmplicons import TwoReadIlluminaRun
from TnAmplicons import IlluminaTwoReadOutput
from TnAmplicons import validateApp


class preprocessApp:
    """
    Preprocess generic Illumina amplicon data
    """

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, output_prefix, primerFile, primerMaxDiff=4, primerEndMatch=4, insertionSite='TA', batchsize=10000, uncompressed=False, output_unidentified=False, minQ=None, minL=0, verbose=True, debug=False, test=False):
        """
        Start preprocessing double barcoded Illumina sequencing run, perform
        """
        self.verbose = verbose
        try:
            # read in primer sequences if present
            prTable = primerTable(primerFile)
            if verbose:
                sys.stdout.write("primer table length P5 Primer Sequences:%s, P7 Primer Sequences:%s, Linker sequences:%s\n" % (len(prTable.getP5sequences()), len(prTable.getP7sequences()), len(prTable.getLinkersequences())))
            # validate assumptions
            # v = validateApp()
            # if v.validatePrimer(prTable, debug) != 0:
            #    sys.stderr.write("Failed validation\n")
            #    self.clean()
            #    return 1

            # setup output files
            identified_count = 0
            unidentified_count = 0
            self.run_out = {}
            self.run_out["Identified"] = IlluminaTwoReadOutput(output_prefix, uncompressed)
            if output_unidentified:
                self.run_out["Unidentified"] = IlluminaTwoReadOutput(output_prefix+"_Unidentified", uncompressed)

            # establish and open the Illumina run
            self.run = TwoReadIlluminaRun(fastq_file1, fastq_file2)
            self.run.open()
            lasttime = time.time()
            while 1:
                # get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                # process individual reads
                for read in reads:
                    read.checkLinker(prTable, 0, primerMaxDiff, primerEndMatch)
                    read.assignPrimer(prTable, 0, primerMaxDiff, primerEndMatch)
                    if minQ is not None:
                        read.trimRead(minQ, minL)
                    if read.goodRead is True:
                        identified_count += 1
                        self.run_out["Identified"].addRead(read.getFastq())
                    else:
                        unidentified_count += 1
                        if output_unidentified:
                            self.run_out["Unidentified"].addRead(read.getFastq())
                # Write out reads
                for key in self.run_out:
                    self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("processed %s total reads, %s Reads/second, %s identified reads(%s%%), %s unidentified reads\n" % (self.run.count(), round(self.run.count()/(time.time() - lasttime), 0), identified_count, round((float(identified_count)/float(self.run.count()))*100, 1), unidentified_count))
                if test:  # exit after the first batch to test the inputs
                    break
            if self.verbose:
                    sys.stdout.write("%s reads processed in %s minutes, %s (%s%%) identified\n\n" % (self.run.count(), round((time.time()-lasttime)/(60), 2), identified_count, round((float(identified_count)/float(self.run.count()))*100, 1)))
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
            self.run.close()
            for key in self.run_out:
                self.run_out[key].close()
        except:
            pass

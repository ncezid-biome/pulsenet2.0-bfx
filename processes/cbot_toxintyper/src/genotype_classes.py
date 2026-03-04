
from collections import namedtuple
from os.path import exists


class GenotypeHit(object):

    def __init__(self):
        self.reference_id = ''
        self.reference_len = ''
        self.query_id = ''
        self._query_start = 0
        self._query_stop = 0
        self._reference_start = 0
        self._reference_stop = 0
        self.forward = True
        self.identity = 0.
        self.absolute_len = 0
        self.relative_len = 0
        self.query_len = 0
        self.num_mismatches = 0
        self.num_gap_opens = 0
        self.query_seq = None
        self.reference_seq = None
        self.bitscore = 0
        self.evalue = 0
        self.full_match = False

    # ************* IMPORTANT *************
    # 
    # The following methods are used to cope
    # with using 1 indexed BLAST indices.
    # If you are going to use the setters
    # make sure you set it with the 0 indexed
    # value, which is reasonable since you are getting
    # the 0 indexed values from the below
    # property methods to begin with.
    # We did this so I can have good feels about this
    # datastructure (thanks).
    # 

    @property
    def query_start(self):
        return self._query_start - 1

    @query_start.setter
    def query_start(self, val):

        self._query_start = val + 1

    @property
    def query_stop(self):
        return self._query_stop - 1

    @query_stop.setter
    def query_stop(self, val):
        self._query_stop = val + 1

    @property
    def reference_start(self):
        return self._reference_start - 1

    @reference_start.setter
    def reference_start(self, val):
        self._reference_start = val + 1

    @property
    def reference_stop(self):
        return self._reference_stop - 1

    @reference_stop.setter
    def reference_stop(self, val):
        self._reference_stop = val + 1

    @staticmethod
    def from_blast(line):
        # query_id  ref_id|len  iden    alignmentlen    mismatches
        #   gapopens    qstart  qstop   refstart    refstop evalue
        #   bitscore    qlen    queryseq    refseq 

        # Create the object
        hit = GenotypeHit()

        # Split the line
        parts = line.split('\t')

        # Parse the results
        hit.query_id = parts[0]
        hit.reference_id = parts[1]
        hit.reference_len = int((parts[1]).split("|")[1])
        hit.identity = float(parts[2])
        hit.absolute_len = int(parts[3])
        hit.num_mismatches = int(parts[4])
        hit.num_gap_opens = int(parts[5])

        # NOTE:
        # The alignments come in 1-indexed
        hit._query_start = int(parts[6])
        hit._query_stop = int(parts[7])
        hit._reference_start = int(parts[8])
        hit._reference_stop = int(parts[9])

        hit.evalue = float(parts[10])
        hit.bitscore = float(parts[11])

        hit.query_len = int(parts[12])
        hit.coverage = float((hit._query_stop-hit._query_start+1)/(hit.reference_len))
        
        # Not all of the genotyping needs to 
        # have the sequences returned
        if len(parts) == 15:
            hit.query_seq = parts[13]
            hit.reference_seq = parts[14]

        # Check to see if the sequence has been reversed
        # If for whatever reason the hit is one nucleotide
        # long, this will be false
        hit.forward = hit._reference_start < hit._reference_stop

        if not hit.forward:
            hit._reference_start, hit._reference_stop = \
                hit._reference_stop, hit._reference_start

        # Change the percentages into decimals,
        # because why not
        hit.identity = hit.identity / 100.

        # Make sure to avoid a divide by zero error
        if hit.reference_len:
            hit.relative_len = float(hit.absolute_len) / float(hit.reference_len)

        if not hit.num_gap_opens and hit.identity == 1.0 and hit.relative_len == 1.0:
            hit.full_match = True

        return hit

class GenotypeResults(object):

    def __init__(self, filename):
        self._hits = []
        self.load_hits(filename)

    def load_hits(self, filename):
        if isinstance(filename, str):
            if exists(filename):
                with open(filename, 'r') as f:
                    for line in f.readlines():
                        self._hits.append(GenotypeHit.from_blast(line))

            else:

                raise RuntimeError('Provide file path is not a real'
                    ' path for alignment results parsing')

        else:
            # We were given a file handle
            for line in filename.readlines():

                self._hits.append(GenotypeHit.from_blast(line))

    @property
    def hits(self):
        return self._hits

SequenceInfo = namedtuple('SequenceInfo', ['locus', 'allele', 'accession', 'sequence', 'other'])
GenotypeRegion = namedtuple('GenotypeRegion', ['coverage', 'identity', 'locations'])
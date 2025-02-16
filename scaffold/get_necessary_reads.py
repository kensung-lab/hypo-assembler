import pysam
import sys

written_list = set()

print("Reading alignments from %s" % sys.argv[1])
aligns = pysam.AlignmentFile(sys.argv[1], "rb")

contig_lens = {}
for i in range(len(aligns.references)):
    contig_lens[aligns.references[i]] = aligns.lengths[i]

print("Writing reads to %s" % sys.argv[3])
output = open(sys.argv[3], "w")
for align in aligns:
    if not align.is_secondary and align.reference_name is not None:
        if align.reference_start <= 10000 or (align.reference_end is not None and align.reference_end >= contig_lens[align.reference_name] - 10000):
            if align.query_name not in written_list:
                print(">%s\n%s" % (align.query_name, align.query_sequence), file=output)
                
print("Reading alignments from %s" % sys.argv[2])
aligns = pysam.AlignmentFile(sys.argv[2], "rb")

contig_lens = {}
for i in range(len(aligns.references)):
    contig_lens[aligns.references[i]] = aligns.lengths[i]
                
for align in aligns:
    if not align.is_secondary and align.reference_name is not None:
        if align.reference_start <= 10000 or (align.reference_end is not None and align.reference_end >= contig_lens[align.reference_name] - 10000):
            if align.query_name not in written_list:
                print(">%s\n%s" % (align.query_name, align.query_sequence), file=output)

output.close()

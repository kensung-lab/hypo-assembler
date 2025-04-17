import sys

from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re
import pysam

IDENTITY_THRESHOLD = 0.9

contig_filename = sys.argv[1]
overlap_filename = sys.argv[2]
output_path = sys.argv[3]

print("Loading contigs")

contigs = {}
for record in SeqIO.parse(contig_filename, "fasta"):
    contigs[record.id] = record

aligns = pysam.AlignmentFile(overlap_filename, "rb")

fully_contained_contigs = set()
for align in aligns:
    if not align.is_secondary and align.reference_name is not None and align.cigartuples is not None:
        qn = align.query_name
        rn = align.reference_name
        if qn != rn:
            total_clips = 0
            for a, b in align.cigartuples:
                if a == 4 or a == 5:
                    total_clips += b
            
            if total_clips == 0 or total_clips <= 0.01 * align.infer_read_length():
                fully_contained_contigs.add(qn)

write_contigs = []
for qn in contigs:
    if qn not in fully_contained_contigs:
        write_contigs.append(contigs[qn])

print("Removed ", len(contigs) - len(write_contigs), " duplicated contigs")

SeqIO.write(write_contigs, output_path, "fasta")

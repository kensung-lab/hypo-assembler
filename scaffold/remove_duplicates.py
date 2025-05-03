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

remove_initial_dupe = True

if len(sys.argv) > 4:
    if sys.argv[4] == '1':
        remove_initial_dupe = False
    

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



sorted_contigs = []
for qn in contigs:
    if qn not in fully_contained_contigs:
        sorted_contigs.append((qn, str(contigs[qn].seq)))
sorted_contigs.sort()

final_contigs = []
previous_qn = None
current_contig = ""
hap_code = "_1"

if remove_initial_dupe:
    count_removed = 0
    for qn, get_contig in sorted_contigs:
        if previous_qn is None:
            previous_qn = qn
            current_contig = get_contig
        else:
            a = "_".join(qn.split("_")[:-2])
            b = "_".join(previous_qn.split("_")[:-2])
            
            hap_code = qn.split("_")[-1]
            if a == b:
                current_contig += get_contig
                count_removed += 1
            else:
                final_contigs.append(current_contig)
                previous_qn = qn
                current_contig = get_contig
    if previous_qn is not None:
        if len(current_contig) >= 10000:
            final_contigs.append(current_contig)
else:
    for qn, get_contig in sorted_contigs:
        if len(get_contig) >= 10000:
            final_contigs.append(get_contig)

cid = 0
write_contigs = []
for contig in final_contigs:
    write_contigs.append(SeqRecord(Seq(contig), id="contig_%d_%s" % (cid, hap_code), description=""))
    cid += 1

print("Removed ", len(contigs) - len(write_contigs), " duplicated contigs")
    
    
SeqIO.write(write_contigs, output_path, "fasta")

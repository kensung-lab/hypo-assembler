import sys

if len(sys.argv) < 4:
    print("Usage:", sys.argv[0], "<pickled file>", "<long read alignments>", "<output prefix>")
    exit(1)

import pickle

with open(sys.argv[1], "rb") as f:
    overlap_positions, processed_contigs, contigs, duplicates = pickle.load(f)

import pysam
from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

aligns = pysam.AlignmentFile(sys.argv[2], "rb")


final_sequences = []
final_aux = []

unmerged_contigs = set()

count = 0

for cn in overlap_positions:
    ctg1, ctg2, lu, rs, re1, re2, c1s, c1e, c2s, c2e, seq1, seq2 = overlap_positions[cn]
    
    count += 1
    
    cn1 = cn + "_1"
    cn2 = cn + "_2"
    
    if lu == 1: # left uses 1, means support_1 is from re1
        support_1 = []
        for align in aligns.fetch(cn1, re1-1, re1+1):
            len20 = align.infer_read_length() // 5
            if len20 >= 400 and align.reference_start <= re1 - len20 and align.reference_end >= re1 + len20:
                support_1.append(align.query_name)
        
        support_2 = []
        for align in aligns.fetch(cn2, rs-1, rs+1):
            len20 = align.infer_read_length() // 5
            if len20 >= 400 and align.reference_start <= rs - len20 and align.reference_end >= rs + len20:
                support_2.append(align.query_name)
                
    else:
        support_1 = []
        for align in aligns.fetch(cn1, rs-1, rs+1):
            len20 = align.infer_read_length() // 5
            if len20 >= 400 and align.reference_start <= rs - len20 and align.reference_end >= rs + len20:
                support_1.append(align.query_name)
            
        support_2 = []
        for align in aligns.fetch(cn2, re2-1, re2+1):
            len20 = align.infer_read_length() // 5
            if len20 >= 400 and align.reference_start <= re2 - len20 and align.reference_end >= re2 + len20:
                support_2.append(align.query_name)
    
            
    if len(support_1) >= 5:
        if len(support_1) > len(support_2) or len(support_2) < 5:
            # take 1
            final_sequences.append(SeqRecord(Seq(seq1), id=cn))
            if lu == 1:
                final_aux.append((cn, ctg1, 0, c1e, ctg2, c2e, len(contigs[ctg2]), support_1))
            else:
                final_aux.append((cn, ctg2, 0, c2s, ctg1, c1s, len(contigs[ctg1]), support_1))
        else:
            # take 2
            final_sequences.append(SeqRecord(Seq(seq2), id=cn))
            if lu == 1:
                final_aux.append((cn, ctg1, 0, c1s, ctg2, c2s, len(contigs[ctg2]), support_2))
            else:
                final_aux.append((cn, ctg2, 0, c2e, ctg1, c1e, len(contigs[ctg1]), support_2))
    elif len(support_2) >= 5:
        # take 2
        final_sequences.append(SeqRecord(Seq(seq2), id=cn))
        if lu == 1:
            final_aux.append((cn, ctg1, 0, c1s, ctg2, c2s, len(contigs[ctg2]), support_2))
        else:
            final_aux.append((cn, ctg2, 0, c2e, ctg1, c1e, len(contigs[ctg1]), support_2))
    else:
        # not taken
        
        unmerged_contigs.add(ctg1)
        unmerged_contigs.add(ctg2)

for ctg in contigs:
    if ctg not in processed_contigs or ctg in unmerged_contigs:
        final_sequences.append(SeqRecord(Seq(contigs[ctg]), id=ctg))

prefix = sys.argv[3]
SeqIO.write(final_sequences, prefix + ".fa", "fasta")

f = open(prefix + ".aux", "w")
for cn, c1, c1s, c1e, c2, c2s, c2e, support in final_aux:
    print(cn, c1, c1s, c1e, c2, c2s, c2e, len(support), sep='\t', file=f)
f.close()

f = open(prefix + ".dup", "w")
for c1, c2 in duplicates:
    print(c1, c2, sep='\t', file=f)
f.close()

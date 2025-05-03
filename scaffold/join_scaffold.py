import sys

if len(sys.argv) < 6:
    print("Usage:", sys.argv[0], "<draft>", "<overlaps>", "<output>", "<pickle temp out for filter>", "<temp directory>", "<draft2>", "<overlaps2>", "<output2>", "<pickle temp2 out for filter>", "<debug_out>")
    exit(1)

from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re

def find_identity(seq1, seq2):
    path1 = os.path.join(temp_directory, "1.fa")
    path2 = os.path.join(temp_directory, "2.fa")
    path3 = os.path.join(temp_directory, "3.paf")
    SeqIO.write([SeqRecord(Seq(seq1), id=c1)], path1, "fasta")
    SeqIO.write([SeqRecord(Seq(seq2), id=c2)], path2, "fasta")
    os.system("minimap2 -cx asm5 " + path1 + " " + path2 + " > " + path3 + " 2> /dev/null")
    f = open(path3)
    total_match = 0
    for line in f:
        a = line.strip().split("\t")
        is_secondary = False
        for part in a[12:]:
            if part.startswith("tp"):
                if part[-1] == 'S':
                    is_secondary = True
        
        if not is_secondary:
            for part in a[12:]:
                if part.startswith("cg:Z:"):
                    cigar = a[-1][5:]
                    matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
                    for a, b in matches:
                        if b == 'M':
                            total_match += int(a)
    
    if min(len(seq1), len(seq2)) > 0:
        return total_match / min(len(seq1), len(seq2))
    else:
        return 0
        
IDENTITY_THRESHOLD = 0.9

contig_filename = sys.argv[1]
overlap_filename = sys.argv[2]
output_path = sys.argv[3]
temp_directory = sys.argv[5]

is_debug = False
if len(sys.argv) > 10:
    is_debug = True
    debug_path = sys.argv[10]

NEW_SEGMENT_THRESHOLD = 1000

print("Loading contigs")

contigs = {}
for record in SeqIO.parse(contig_filename, "fasta"):
    contigs[record.id] = str(record.seq)


print("Loading overlaps and resolving ordering")

contig_pairs = {}
per_contig = {}
for line in open(overlap_filename):
    a = line.strip().split("\t")
    if a[0] == "R":
        c1 = a[1]
        c2 = a[3]
        orientation = int(a[2]) ^ int(a[4])
        score = int(a[5])
        
        if c1 not in per_contig:
            per_contig[c1] = []
        per_contig[c1].append((c2, orientation, score))

for c1 in per_contig:
    per_contig[c1].sort(key = lambda x: -x[2])
    
    c2, orientation, score = per_contig[c1][0]
        
    if c1 != c2 and (c1, c2) not in contig_pairs:
        contig_pairs[(c1, c2)] = []
    contig_pairs[(c1, c2)].append((c1, c2, score, orientation))
    
pair_to_process = []
for c1, c2 in contig_pairs:
    s = 0
    for a in contig_pairs[(c1, c2)]:
        s += a[2]
    pair_to_process.append((s, c1, c2))
pair_to_process.sort()

print("Joining")

if is_debug:
    debug_file = open(debug_path, "w")

processed_contigs = set()

get_sequences = []
get_segments = []

duplicates = set()
duplicated_contigs = set()

overlap_positions = {}

construct_results = []

for p_idx in range(len(pair_to_process)):
    s, c1, c2 = pair_to_process[len(pair_to_process)-p_idx-1]
    
    if c1 != "NA" and c2 != "NA" and c1 not in processed_contigs and c2 not in processed_contigs:
        
        new_segments = []
        
        # now join the contigs
        
        temp_locs = []
        
        # check orientation
        consistent_orientation = True
        current_orientation = None
        for c1, c2, score, orientation in contig_pairs[(c1, c2)]:
            if current_orientation is not None and orientation != current_orientation:
                consistent_orientation = False
                break
            current_orientation = orientation
            
        if consistent_orientation:
            c1_seq = contigs[c1]
            if current_orientation == -1:
                c2_seq = str(Seq(contigs[c2]).reverse_complement())
            else:
                c2_seq = contigs[c2]
                
            construct_result = c1_seq + c2_seq
            
            record_name = c1 + "_" + c2 + "_" + str(orientation)
            construct_results.append(SeqRecord(Seq(construct_result), id=record_name))
            
            processed_contigs.add(c1)
            processed_contigs.add(c2)

for ctg in contigs:
    if ctg not in processed_contigs:
        construct_results.append(SeqRecord(Seq(contigs[ctg]), id=ctg))

old_construct_results = construct_results
construct_results = []
old_construct_results.sort(key=lambda x: len(x))
for i in range(len(old_construct_results) // 2):
    r1 = old_construct_results[i]
    r2 = old_construct_results[len(old_construct_results) - i - 1]
    
    record_name = r1.id + "_" + r2.id
    new_seq = r1.seq + r2.seq
    
    construct_results.append(SeqRecord(Seq(new_seq), id=record_name))

print(len(old_construct_results), len(construct_results))

SeqIO.write(construct_results, output_path, "fasta")

import pickle

with open(sys.argv[4], "wb") as f:
    pickle.dump([construct_results, processed_contigs, contigs, duplicates], f)

contig_filename = sys.argv[6]
overlap_filename = sys.argv[7]
output_path = sys.argv[8]
temp_directory = sys.argv[5]

print("Loading contigs")

contigs = {}
for record in SeqIO.parse(contig_filename, "fasta"):
    contigs[record.id] = str(record.seq)
    
print("Loading overlaps and resolving ordering")

contig_pairs = {}
per_contig = {}
for line in open(overlap_filename):
    a = line.strip().split("\t")
    if a[0] == "R":
        c1 = a[1]
        c2 = a[3]
        orientation = int(a[2]) ^ int(a[4])
        score = int(a[5])
        
        if c1 not in per_contig:
            per_contig[c1] = []
        per_contig[c1].append((c2, orientation, score))

for c1 in per_contig:
    per_contig[c1].sort(key = lambda x: -x[2])
    
    c2, orientation, score = per_contig[c1][0]
        
    if c1 != c2 and (c1, c2) not in contig_pairs:
        contig_pairs[(c1, c2)] = []
    contig_pairs[(c1, c2)].append((c1, c2, score, orientation))
    
pair_to_process = []
for c1, c2 in contig_pairs:
    s = 0
    for a in contig_pairs[(c1, c2)]:
        s += a[2]
    pair_to_process.append((s, c1, c2))
pair_to_process.sort()

print("Joining")

if is_debug:
    debug_file = open(debug_path, "w")

processed_contigs = set()

get_sequences = []
get_segments = []

duplicates = set()
duplicated_contigs = set()

overlap_positions = {}

construct_results = []

for p_idx in range(len(pair_to_process)):
    s, c1, c2 = pair_to_process[len(pair_to_process)-p_idx-1]
    
    if c1 != "NA" and c2 != "NA" and c1 not in processed_contigs and c2 not in processed_contigs:
        
        new_segments = []
        
        # now join the contigs
        
        temp_locs = []
        
        # check orientation
        consistent_orientation = True
        current_orientation = None
        for c1, c2, score, orientation in contig_pairs[(c1, c2)]:
            if current_orientation is not None and orientation != current_orientation:
                consistent_orientation = False
                break
            current_orientation = orientation
            
        if consistent_orientation:
            c1_seq = contigs[c1]
            if current_orientation == -1:
                c2_seq = str(Seq(contigs[c2]).reverse_complement())
            else:
                c2_seq = contigs[c2]
                
            construct_result = c1_seq + c2_seq
            
            record_name = c1 + "_" + c2 + "_" + str(orientation)
            construct_results.append(SeqRecord(Seq(construct_result), id=record_name))
            
            processed_contigs.add(c1)
            processed_contigs.add(c2)

for ctg in contigs:
    if ctg not in processed_contigs:
        construct_results.append(SeqRecord(Seq(contigs[ctg]), id=ctg))

old_construct_results = construct_results
construct_results = []
old_construct_results.sort(key=lambda x: len(x))
for i in range(len(old_construct_results) // 2):
    r1 = old_construct_results[i]
    r2 = old_construct_results[len(old_construct_results) - i - 1]
    
    record_name = r1.id + "_" + r2.id
    new_seq = r1.seq + r2.seq
    
    construct_results.append(SeqRecord(Seq(new_seq), id=record_name))

SeqIO.write(construct_results, output_path, "fasta")

with open(sys.argv[9], "wb") as f:
    pickle.dump([construct_results, processed_contigs, contigs, duplicates], f)

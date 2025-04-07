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
for line in open(overlap_filename):
    a = line.strip().split("\t")
    c1 = a[0]
    c2 = a[3]
    c1_s = int(a[1])
    c1_e = int(a[2])
    c2_s = int(a[4])
    c2_e = int(a[5])
    score = int(a[6])
    orientation = int(a[7])
    
    if (c1, c2) not in contig_pairs:
        contig_pairs[(c1, c2)] = []
    contig_pairs[(c1, c2)].append((c1, c1_s, c1_e, c2, c2_s, c2_e, score, orientation))
    
pair_to_process = []
for c1, c2 in contig_pairs:
    s = 0
    for a in contig_pairs[(c1, c2)]:
        s += a[6]
    pair_to_process.append((s, c1, c2))
pair_to_process.sort()

print("Joining")

output_file = open(output_path, "w")

if is_debug:
    debug_file = open(debug_path, "w")

processed_contigs = set()

get_sequences = []
get_segments = []

duplicates = set()
duplicated_contigs = set()

overlap_positions = {}

for p_idx in range(len(pair_to_process)):
    s, c1, c2 = pair_to_process[len(pair_to_process)-p_idx-1]
    
    if c1 != "NA" and c2 != "NA": #and c1 not in processed_contigs and c2 not in processed_contigs:
        
        new_segments = []
        
        # now join the contigs
        
        temp_locs = []
        
        # check orientation
        consistent_orientation = True
        current_orientation = None
        for c1, c1_s, c1_e, c2, c2_s, c2_e, score, orientation in contig_pairs[(c1, c2)]:
            if current_orientation is not None and orientation != current_orientation:
                consistent_orientation = False
                break
            current_orientation = orientation
            
            if current_orientation == -1:
                temp_locs.append((c1_s, c1_e, c2_e, c2_s))
            else:
                temp_locs.append((c1_s, c1_e, c2_s, c2_e))
        
        temp_locs.sort()
        consistent_ordering = True
        for i in range(1, len(temp_locs)):
            if temp_locs[i][0] < temp_locs[i-1][1] or temp_locs[i][2] < temp_locs[i-1][3]:
                consistent_ordering = False
                break
        
        if consistent_orientation and consistent_ordering:
            #print("Joining", c1, "and", c2)
            
            c1_locs = []
            c2_locs = []
            c1_seq = None
            c2_seq = None
            for a, b, c, d in temp_locs:
                c1_locs.append(a)
                c1_locs.append(b)
                c2_locs.append(c)
                c2_locs.append(d)
                
            c1_seq = contigs[c1]
            if current_orientation == -1:
                c2_seq = str(Seq(contigs[c2]).reverse_complement())
            else:
                c2_seq = contigs[c2]
                
            c1_overlap_start = min(c1_locs)
            c1_overlap_end = max(c1_locs)
            
            c2_overlap_start = min(c2_locs)
            c2_overlap_end = max(c2_locs)
            
            construct_result = ""
            
            c1_left_remainder = c1_seq[:c1_overlap_start]
            c2_left_remainder = c2_seq[:c2_overlap_start]
            
            left_used = None
            
            if len(c1_left_remainder) < len(c2_left_remainder):
                left_used = 2
                construct_result += c2_left_remainder
            else:
                left_used = 1
                construct_result += c1_left_remainder
                    
            ov_start = len(construct_result)
            
            construct_result_1 = construct_result + c1_seq[c1_overlap_start:c1_overlap_end]
            construct_result_2 = construct_result + c2_seq[c2_overlap_start:c2_overlap_end]
            
            ov_end1 = len(construct_result_1)
            ov_end2 = len(construct_result_2)
            
            # rightmost, after last index
            c1_right_remainder = c1_seq[c1_overlap_end:]
            c2_right_remainder = c2_seq[c2_overlap_end:]
            
            right_used = None
            
            # right
            if len(c1_right_remainder) > len(c2_right_remainder):
                right_used = 1
                construct_result_1 += c1_right_remainder
                construct_result_2 += c1_right_remainder
            else:
                right_used = 2
                construct_result_1 += c2_right_remainder
                construct_result_2 += c2_right_remainder
            
            c1_part = c1_seq[c1_overlap_start:c1_overlap_end]
            c2_part = c2_seq[c2_overlap_start:c2_overlap_end]
            
            new_segments.append((c1_part, c2_part))
            
            if right_used != left_used:
                
                if c1 not in processed_contigs and c2 not in processed_contigs:
                    # check for similarity
                    identity = find_identity(c1_part, c2_part)
                    
                    if identity >= IDENTITY_THRESHOLD:
                    
                        # new contigs are done, now write 

                        processed_contigs.add(c1)
                        processed_contigs.add(c2)

                        print("Joining", c1, c2, "with identity", identity)

                        record_name = c1 + "_" + c2 + "_" + str(orientation)
                        SeqIO.write(SeqRecord(Seq(construct_result_1), id=record_name + "_1"), output_file, "fasta")
                        SeqIO.write(SeqRecord(Seq(construct_result_2), id=record_name + "_2"), output_file, "fasta")

                        overlap_positions[record_name] = (c1, c2, left_used, ov_start, ov_end1, ov_end2, c1_overlap_start, c1_overlap_end, c2_overlap_start, c2_overlap_end, construct_result_1, construct_result_2)
                    else:
                        print(c1, "and", c2, "considered overlap but identity is", identity)
            else:
                if right_used == 1: # for duplicates, do not add to processed_contigs the actual dupes
                    if len(c2_part) >= 0.8 * len(c2_seq):
                        identity = find_identity(c1_part, c2_part)
                        if identity >= IDENTITY_THRESHOLD:
                            processed_contigs.add(c2)
                            duplicates.add((c2, c1))
                            print(c2, "is duplicate of", c1, "with identity", identity)
                        else:
                            print(c2, "is considered duplicate of", c1, "but identity is", identity)
                    else:
                        print(c2, "is considered duplicate of", c1, "but overlap is too short")
                else:
                    if len(c1_part) >= 0.8 * len(c1_seq):
                        identity = find_identity(c1_part, c2_part)
                        if identity >= IDENTITY_THRESHOLD:
                            processed_contigs.add(c1)
                            duplicates.add((c1, c2))
                            print(c1, "is duplicate of", c2, "with identity", identity)
                        else:
                            print(c1, "is considered duplicate of", c2, "but identity is", identity)
                    else:
                        print(c1, "is considered duplicate of", c2, "but overlap is too short")
            if is_debug:
                SeqIO.write(new_segments, debug_file, "fasta")

for ctg in contigs:
    if ctg not in processed_contigs:
        SeqIO.write(SeqRecord(Seq(contigs[ctg]), id=ctg), output_file, "fasta")
output_file.close()

import pickle

with open(sys.argv[4], "wb") as f:
    pickle.dump([overlap_positions, processed_contigs, contigs, duplicates], f)

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
for line in open(overlap_filename):
    a = line.strip().split("\t")
    c1 = a[0]
    c2 = a[3]
    c1_s = int(a[1])
    c1_e = int(a[2])
    c2_s = int(a[4])
    c2_e = int(a[5])
    score = int(a[6])
    orientation = int(a[7])
    
    if (c1, c2) not in contig_pairs:
        contig_pairs[(c1, c2)] = []
    contig_pairs[(c1, c2)].append((c1, c1_s, c1_e, c2, c2_s, c2_e, score, orientation))
    
pair_to_process = []
for c1, c2 in contig_pairs:
    s = 0
    for a in contig_pairs[(c1, c2)]:
        s += a[6]
    pair_to_process.append((s, c1, c2))
pair_to_process.sort()

print("Joining")

output_file = open(output_path, "w")

processed_contigs = set()

get_sequences = []
get_segments = []

duplicates = set()
duplicated_contigs = set()

overlap_positions = {}

for p_idx in range(len(pair_to_process)):
    s, c1, c2 = pair_to_process[len(pair_to_process)-p_idx-1]
    
    if c1 != "NA" and c2 != "NA": #and c1 not in processed_contigs and c2 not in processed_contigs:
        
        new_segments = []
        
        # now join the contigs
        
        temp_locs = []
        
        # check orientation
        consistent_orientation = True
        current_orientation = None
        for c1, c1_s, c1_e, c2, c2_s, c2_e, score, orientation in contig_pairs[(c1, c2)]:
            if current_orientation is not None and orientation != current_orientation:
                consistent_orientation = False
                break
            current_orientation = orientation
            
            if current_orientation == -1:
                temp_locs.append((c1_s, c1_e, c2_e, c2_s))
            else:
                temp_locs.append((c1_s, c1_e, c2_s, c2_e))
        
        temp_locs.sort()
        consistent_ordering = True
        for i in range(1, len(temp_locs)):
            if temp_locs[i][0] < temp_locs[i-1][1] or temp_locs[i][2] < temp_locs[i-1][3]:
                consistent_ordering = False
                break
        
        if consistent_orientation and consistent_ordering:
            #print("Joining", c1, "and", c2)
            
            c1_locs = []
            c2_locs = []
            c1_seq = None
            c2_seq = None
            for a, b, c, d in temp_locs:
                c1_locs.append(a)
                c1_locs.append(b)
                c2_locs.append(c)
                c2_locs.append(d)
                
            c1_seq = contigs[c1]
            if current_orientation == -1:
                c2_seq = str(Seq(contigs[c2]).reverse_complement())
            else:
                c2_seq = contigs[c2]
                
            c1_overlap_start = min(c1_locs)
            c1_overlap_end = max(c1_locs)
            
            c2_overlap_start = min(c2_locs)
            c2_overlap_end = max(c2_locs)
            
            construct_result = ""
            
            c1_left_remainder = c1_seq[:c1_overlap_start]
            c2_left_remainder = c2_seq[:c2_overlap_start]
            
            left_used = None
            
            if len(c1_left_remainder) < len(c2_left_remainder):
                left_used = 2
                construct_result += c2_left_remainder
            else:
                left_used = 1
                construct_result += c1_left_remainder
                    
            ov_start = len(construct_result)
            
            construct_result_1 = construct_result + c1_seq[c1_overlap_start:c1_overlap_end]
            construct_result_2 = construct_result + c2_seq[c2_overlap_start:c2_overlap_end]
            
            ov_end1 = len(construct_result_1)
            ov_end2 = len(construct_result_2)
            
            # rightmost, after last index
            c1_right_remainder = c1_seq[c1_overlap_end:]
            c2_right_remainder = c2_seq[c2_overlap_end:]
            
            right_used = None
            
            # right
            if len(c1_right_remainder) > len(c2_right_remainder):
                right_used = 1
                construct_result_1 += c1_right_remainder
                construct_result_2 += c1_right_remainder
            else:
                right_used = 2
                construct_result_1 += c2_right_remainder
                construct_result_2 += c2_right_remainder
            
            c1_part = c1_seq[c1_overlap_start:c1_overlap_end]
            c2_part = c2_seq[c2_overlap_start:c2_overlap_end]
            
            new_segments.append((c1_part, c2_part))
            
            if right_used != left_used:
                
                if c1 not in processed_contigs and c2 not in processed_contigs:
                    # check for similarity
                    identity = find_identity(c1_part, c2_part)
                    
                    if identity >= IDENTITY_THRESHOLD:
                    
                        # new contigs are done, now write 

                        processed_contigs.add(c1)
                        processed_contigs.add(c2)

                        print("Joining", c1, c2, "with identity", identity)

                        record_name = c1 + "_" + c2 + "_" + str(orientation)
                        SeqIO.write(SeqRecord(Seq(construct_result_1), id=record_name + "_1"), output_file, "fasta")
                        SeqIO.write(SeqRecord(Seq(construct_result_2), id=record_name + "_2"), output_file, "fasta")

                        overlap_positions[record_name] = (c1, c2, left_used, ov_start, ov_end1, ov_end2, c1_overlap_start, c1_overlap_end, c2_overlap_start, c2_overlap_end, construct_result_1, construct_result_2)
                    else:
                        print(c1, "and", c2, "considered overlap but identity is", identity)
            else:
                if right_used == 1: # for duplicates, do not add to processed_contigs the actual dupes
                    if len(c2_part) >= 0.8 * len(c2_seq):
                        identity = find_identity(c1_part, c2_part)
                        if identity >= IDENTITY_THRESHOLD:
                            processed_contigs.add(c2)
                            duplicates.add((c2, c1))
                            print(c2, "is duplicate of", c1, "with identity", identity)
                        else:
                            print(c2, "is considered duplicate of", c1, "but identity is", identity)
                    else:
                        print(c2, "is considered duplicate of", c1, "but overlap is too short")
                else:
                    if len(c1_part) >= 0.8 * len(c1_seq):
                        identity = find_identity(c1_part, c2_part)
                        if identity >= IDENTITY_THRESHOLD:
                            processed_contigs.add(c1)
                            duplicates.add((c1, c2))
                            print(c1, "is duplicate of", c2, "with identity", identity)
                        else:
                            print(c1, "is considered duplicate of", c2, "but identity is", identity)
                    else:
                        print(c1, "is considered duplicate of", c2, "but overlap is too short")
            if is_debug:
                SeqIO.write(new_segments, debug_file, "fasta")

for ctg in contigs:
    if ctg not in processed_contigs:
        SeqIO.write(SeqRecord(Seq(contigs[ctg]), id=ctg), output_file, "fasta")
output_file.close()

if is_debug:
    debug_file.close()

import pickle

with open(sys.argv[9], "wb") as f:
    pickle.dump([overlap_positions, processed_contigs, contigs, duplicates], f)

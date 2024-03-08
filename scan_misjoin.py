import pysam
import sys

from Bio import SeqIO
from Bio.Seq import Seq                                                 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if len(sys.argv) < 4:
    print("Usage: python scan_misjoin.py <initial assembly> <reads alignment file> <output file>")
    exit 1

samfile = pysam.AlignmentFile(sys.argv[2], "rb")

count_softclips_left = {}
count_softclips_right = {}
count_coverage_left = {}
count_coverage_right = {}
count_softclips = {}

for i in range(len(samfile.references)):
    count_softclips_left[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 1)]
    count_softclips_right[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 1)]
    count_coverage_left[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 1)]
    count_coverage_right[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 1)]
    count_softclips[samfile.references[i]] = [0 for _ in range(samfile.lengths[i] // 1000 + 1)]

current_rn = "a"
current_rs = 0

samfile = pysam.AlignmentFile(sys.argv[2], "rb")

for align in samfile:
    if not align.is_secondary:
        get_read_length = align.infer_read_length()
        
        # get left-most position
        
        rs = align.reference_start
        re = align.reference_end
        rn = align.reference_name
        
        if current_rn != rn or current_rs - rs >= 1000000:
            current_rn = rn
            current_rs = rs
            print("Processing", rn, rs)
            
        
        if re is not None:
            
            # check for clips
            
            get_cigar = align.cigartuples
            
            if get_cigar is not None:
                # left clips
                left_clip = 0
                start = 0
                while get_cigar[start][0] == 4 or get_cigar[start][0] == 5:
                    left_clip += get_cigar[start][1]
                    start += 1
                    
                    if start >= len(get_cigar):
                        break
                
                # right clips
                right_clip = 0
                start = len(get_cigar) - 1
                while get_cigar[start][0] == 4 or get_cigar[start][0] == 5:
                    right_clip += get_cigar[start][1]
                    start -= 1
                    
                    if start < 0:
                        break
                
                if left_clip >= 0.1 * get_read_length or left_clip >= 500:
                    count_softclips_left[rn][rs // 1000] += 1
                    count_softclips[rn][rs // 1000] += 1
                count_coverage_left[rn][rs // 1000] += 1
                
                if right_clip >= 0.1 * get_read_length or right_clip >= 500:
                    count_softclips_right[rn][re // 1000] += 1
                    count_softclips[rn][rs // 1000] += 1
                count_coverage_right[rn][re // 1000] += 1
                
                

contigs = {}
records = SeqIO.parse(sys.argv[1], "fasta")
for record in records:
    contigs[record.id] = str(record.seq)

new_contig_records = []
count = 0
for ctg_name in contigs:
    
    new_contigs = []
    current_contig = ""
    
    for i in range(len(count_softclips[ctg_name])):
        if i >= 2 and count_coverage_left[ctg_name][i] >= 20 and count_softclips_left[ctg_name][i] >= 0.5 * count_coverage_left[ctg_name][i]:
            if len(current_contig) > 0:
                new_contigs.append(current_contig)
            current_contig = ""
            count += 1
        elif i < len(count_softclips[ctg_name]) - 1 and count_coverage_right[ctg_name][i] >= 20 and count_softclips_right[ctg_name][i] >= 0.5 * count_coverage_right[ctg_name][i]:
            if len(current_contig) > 0:
                new_contigs.append(current_contig)
            current_contig = ""
            count += 1
        else:
            current_contig += contigs[ctg_name][i*1000:(i+1)*1000]
    
    new_contigs.append(current_contig)
    
    print(ctg_name, len(new_contigs))
    
    i = 1
    for ctg in new_contigs:
        new_contig_records.append(SeqRecord(
            Seq(ctg),
            id="%s_%d" % (ctg_name, i),
            name="",
            description="",
        ))
        i += 1

SeqIO.write(new_contig_records, sys.argv[3], "fasta")

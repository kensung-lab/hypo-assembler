from Bio import SeqIO
import sys

if len(sys.argv) < 9:
    print("Usage:", sys.argv[0], "<contig file 1>", "<contig file 2>", "<quast contig 1 to paternal>", "<quast contig 1 to maternal>", "<quast contig 2 to paternal>", "<quast contig 2 to maternal>", "<output 1>" , "<output 2>")
    exit(1)


contigs = sys.argv[1]

contig_names = set()
print("Reading contigs from", contigs)
for record in SeqIO.parse(contigs, "fasta"):
    contig_names.add(record.id)
    
contigs2 = sys.argv[2]

contig_names2 = set()
print("Reading contigs from", contigs2)
for record in SeqIO.parse(contigs2, "fasta"):
    contig_names2.add(record.id)


quast1 = sys.argv[3]
print("Reading contig 1 misassemblies-paternal:", quast1)
f = open(quast1)

current_contig = None
mis_by_contigs = {}
for contig_name in contig_names:
    mis_by_contigs[contig_name] = []

for line in f:
    if line.strip() in contig_names:
        current_contig = line.strip()
    else:
        a = line.strip().split()
        if a[0] == "Extensive":
            s1, e1, s2, e2 = int(a[-5]), int(a[-4]), int(a[-2]), int(a[-1])
            
            a = max(s1, e1)
            b = min(s2, e2)
            
            mis_by_contigs[contig_name].append((a, b))
            
quast2 = sys.argv[4]
print("Reading contig 1 misassemblies-maternal:", quast2)
f = open(quast2)

current_contig = None
mis_by_contigs2 = {}
for contig_name in contig_names:
    mis_by_contigs2[contig_name] = []
    
for line in f:
    if line.strip() in contig_names:
        current_contig = line.strip()
    else:
        a = line.strip().split()
        if a[0] == "Extensive":
            s1, e1, s2, e2 = int(a[-5]), int(a[-4]), int(a[-2]), int(a[-1])
            
            a = max(s1, e1)
            b = min(s2, e2)
            
            mis_by_contigs2[contig_name].append((a, b))
            
quast1 = sys.argv[5]
print("Reading contig 2 misassemblies-paternal:", quast1)
f = open(quast1)

current_contig = None
mis_by_contigs_2 = {}
for contig_name in contig_names:
    mis_by_contigs_2[contig_name] = []

for line in f:
    if line.strip() in contig_names:
        current_contig = line.strip()
    else:
        a = line.strip().split()
        if a[0] == "Extensive":
            s1, e1, s2, e2 = int(a[-5]), int(a[-4]), int(a[-2]), int(a[-1])
            
            a = max(s1, e1)
            b = min(s2, e2)
            
            mis_by_contigs_2[contig_name].append((a, b))
            
quast2 = sys.argv[6]
print("Reading contig 2 misassemblies-maternal:", quast2)
f = open(quast2)

current_contig = None
mis_by_contigs2_2 = {}
for contig_name in contig_names:
    mis_by_contigs2_2[contig_name] = []
    
for line in f:
    if line.strip() in contig_names:
        current_contig = line.strip()
    else:
        a = line.strip().split()
        if a[0] == "Extensive":
            s1, e1, s2, e2 = int(a[-5]), int(a[-4]), int(a[-2]), int(a[-1])
            
            a = max(s1, e1)
            b = min(s2, e2)
            
            mis_by_contigs2_2[contig_name].append((a, b))

print("Resolving misassemblies")
mis_by_contigs_filtered  = {}

for contig_name in mis_by_contigs:
    coords = []
    
    for s1, e1 in mis_by_contigs[contig_name]:
        found = False
        for s2, e2 in mis_by_contigs2[contig_name]:
            if s1 <= e2 and s2 <= e1:
                found = True
        if found:
            coords.append((s1, e1))
    
    for s2, e2 in mis_by_contigs2[contig_name]:
        found = False
        for s1, e1 in mis_by_contigs[contig_name]:
            if s1 <= e2 and s2 <= e1:
                found = True
        if found:
            coords.append((s2, e2))
    
    coords.sort()
    
    mis_by_contigs_filtered[contig_name] = []
    
    for s, e in coords:
        if len(mis_by_contigs_filtered[contig_name]) == 0:
            mis_by_contigs_filtered[contig_name].append((s, e))
        else:
            ps, pe = mis_by_contigs_filtered[contig_name][-1]
            
            if s <= pe:
                mis_by_contigs_filtered[contig_name][-1] = (ps, max(e, pe))
            else:
                mis_by_contigs_filtered[contig_name].append((s, e))
                
mis_by_contigs_filtered_2  = {}

for contig_name in mis_by_contigs:
    coords = []
    
    for s1, e1 in mis_by_contigs_2[contig_name]:
        found = False
        for s2, e2 in mis_by_contigs2_2[contig_name]:
            if s1 <= e2 and s2 <= e1:
                found = True
        if found:
            coords.append((s1, e1))
    
    for s2, e2 in mis_by_contigs2_2[contig_name]:
        found = False
        for s1, e1 in mis_by_contigs_2[contig_name]:
            if s1 <= e2 and s2 <= e1:
                found = True
        if found:
            coords.append((s2, e2))
    
    coords.sort()
    
    mis_by_contigs_filtered_2[contig_name] = []
    
    for s, e in coords:
        if len(mis_by_contigs_filtered_2[contig_name]) == 0:
            mis_by_contigs_filtered_2[contig_name].append((s, e))
        else:
            ps, pe = mis_by_contigs_filtered_2[contig_name][-1]
            
            if s <= pe:
                mis_by_contigs_filtered_2[contig_name][-1] = (ps, max(e, pe))
            else:
                mis_by_contigs_filtered_2[contig_name].append((s, e))

count_mis = 0
print("Writing filtered misassemblies (contig 1) to", sys.argv[7])
output_mis = open(sys.argv[7], "w")
for contig_name in mis_by_contigs:
    for s, e in mis_by_contigs_filtered[contig_name]:
        count_mis += 1
        print("%s\t%d\t%d" % (contig_name, s, e), file=output_mis)
output_mis.close()

count_mis_2 = 0
print("Writing filtered misassemblies (contig 2) to", sys.argv[8])
output_mis = open(sys.argv[8], "w")
for contig_name in mis_by_contigs:
    for s, e in mis_by_contigs_filtered_2[contig_name]:
        count_mis_2 += 1
        print("%s\t%d\t%d" % (contig_name, s, e), file=output_mis)
output_mis.close()

print("Number of misassemblies (contig 1): %d" % count_mis)
print("Number of misassemblies (contig 2): %d" % count_mis_2)

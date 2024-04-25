import subprocess
import argparse
import os
import re
from Bio import SeqIO

def separate_haplotypes(contig_path_1, contig_path_2, paf_path_p_1, paf_path_m_1, paf_path_p_2, paf_path_m_2, out_path_1, out_path_2):
    records = SeqIO.parse(contig_path_1, "fasta")
    pat = {}
    for record in records:
        pat[record.id] = record
        
    records = SeqIO.parse(contig_path_2, "fasta")
    mat = {}
    for record in records:
        mat[record.id] = record
        
    f = open(paf_path_p_1)
    contig_to_score_p1 = {}
    for line in f:
        a = line.strip().split("\t")
        current_nm = None
        current_len = None
        is_primary = False
        for i in range(len(a)):
            z = a[i].split(":")
            if z[0] == "tp" and z[2] == "P":
                is_primary = True
            if z[0] == "NM":
                current_nm = int(z[2])
                current_len = int(a[3]) - int(a[2])
            if z[0] == "cg":
                get_cigar = z[2]

        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score_p1:
                contig_to_score_p1[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score_p1[a[0]]
                contig_to_score_p1[a[0]] = (cnm + current_nm, cnm + current_len)
                
    f = open(paf_path_p_2)
    contig_to_score_p2 = {}
    for line in f:
        a = line.strip().split("\t")
        current_nm = None
        current_len = None
        is_primary = False
        for i in range(len(a)):
            z = a[i].split(":")
            if z[0] == "tp" and z[2] == "P":
                is_primary = True
            if z[0] == "NM":
                current_nm = int(z[2])
                current_len = int(a[3]) - int(a[2])
            if z[0] == "cg":
                get_cigar = z[2]

        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score_p2:
                contig_to_score_p2[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score_p2[a[0]]
                contig_to_score_p2[a[0]] = (cnm + current_nm, cnm + current_len)
                
    f = open(paf_path_m_1)
    contig_to_score_m1 = {}
    for line in f:
        a = line.strip().split("\t")
        current_nm = None
        current_len = None
        is_primary = False
        for i in range(len(a)):
            z = a[i].split(":")
            if z[0] == "tp" and z[2] == "P":
                is_primary = True
            if z[0] == "NM":
                current_nm = int(z[2])
                current_len = int(a[3]) - int(a[2])
            if z[0] == "cg":
                get_cigar = z[2]

        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score_m1:
                contig_to_score_m1[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score_m1[a[0]]
                contig_to_score_m1[a[0]] = (cnm + current_nm, cnm + current_len)
                
    f = open(paf_path_m_2)
    contig_to_score_m2 = {}
    for line in f:
        a = line.strip().split("\t")
        current_nm = None
        current_len = None
        is_primary = False
        for i in range(len(a)):
            z = a[i].split(":")
            if z[0] == "tp" and z[2] == "P":
                is_primary = True
            if z[0] == "NM":
                current_nm = int(z[2])
                current_len = int(a[3]) - int(a[2])
            if z[0] == "cg":
                get_cigar = z[2]

        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score_m2:
                contig_to_score_m2[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score_m2[a[0]]
                contig_to_score_m2[a[0]] = (cnm + current_nm, cnm + current_len)
                
    s1 = []
    s2 = []

    for cid in pat:
        p1len = 0
        p1score = 0
        if cid in contig_to_score_p1:
            p1len = contig_to_score_p1[cid][1]
            p1score = contig_to_score_p1[cid][0]
        p2len = 0
        p2score = 0
        if cid in contig_to_score_p2:
            p2len = contig_to_score_p2[cid][1]
            p2score = contig_to_score_p2[cid][0]
        m1len = 0
        m1score = 0
        if cid in contig_to_score_m1:
            m1len = contig_to_score_m1[cid][1]
            m1score = contig_to_score_m1[cid][0]
        m2len = 0
        m2score = 0
        if cid in contig_to_score_m2:
            m2len = contig_to_score_m2[cid][1]
            m2score = contig_to_score_m2[cid][0]
        
        tpm = (p1len - p1score) + (m2len - m2score)
        tmp = (m1len - m1score) + (p2len - p2score)
        
        if tpm != tmp:
            print(cid, tpm, tmp)
        
        if tpm >= tmp:
            s1.append(pat[cid])
            s2.append(mat[cid])
        else:
            s1.append(mat[cid])
            s2.append(pat[cid])
            
    SeqIO.write(s1, out_path_1, "fasta")
    SeqIO.write(s2, out_path_2, "fasta")
    

def separate(args):
    if args.assembly1 is None:
        print("-1 or --assembly1 is needed for diploid error evaluation.")
        exit(1)
    if args.assembly2 is None:
        print("-2 or --assembly2 is needed for diploid error evaluation.")
        exit(1)
    if args.paternal is None:
        print("-p or --paternal is needed for diploid error evaluation.")
        exit(1)
    if args.maternal is None:
        print("-m or --maternal is needed for diploid error evaluation.")
        exit(1)
    
    if not os.path.exists(args.workdir):
        os.makedirs(args.workdir)
    
    error_workdir = os.path.join(args.workdir, "pafs")
    if not os.path.exists(error_workdir):
        os.makedirs(error_workdir)
    
    # Run minimap2 for evaluation
    path = os.path.join(error_workdir, "asm1_p.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly1, args.paternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "--eqx", "-t", str(args.threads), args.paternal, args.assembly1], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm1_m.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly1, args.maternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "--eqx", "-t", str(args.threads), args.maternal, args.assembly1], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm2_p.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly2, args.paternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "--eqx", "-t", str(args.threads), args.paternal, args.assembly2], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm2_m.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly2, args.maternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "--eqx", "-t", str(args.threads), args.maternal, args.assembly2], check=True, stdout=f)
        f.close()
        
    separate_haplotypes(args.assembly1, args.assembly2, os.path.join(error_workdir, "asm1_p.paf"), os.path.join(error_workdir, "asm1_m.paf"), os.path.join(error_workdir, "asm2_p.paf"), os.path.join(error_workdir, "asm2_m.paf"), args.output + "_1.fa", args.output + "_2.fa")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--assembly1", type=str, help="First haplotype input")
    parser.add_argument("-2", "--assembly2", type=str, help="Second haplotype input")
    parser.add_argument("-p", "--paternal", type=str, help="First haplotype reference")
    parser.add_argument("-m", "--maternal", type=str, help="Second haplotype reference")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to run", default=1)
    parser.add_argument("-w", "--workdir", type=str, help="Number of threads to run", default="hypo_eval_wd")
    parser.add_argument("-o", "--output", type=str, help="Output prefix", default="H")
    args = parser.parse_args()
    
    separate(args)
  
  
if __name__ == '__main__':
    main()

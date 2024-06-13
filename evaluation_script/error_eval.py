import subprocess
import argparse
import os
import re
from Bio import SeqIO

def evaluate_error(contig_path, paf_path_1, paf_path_2):
    contigs = contig_path
    print("Evaluating contigs from", contig_path, flush=True)
    contigs1_p = {}
    contigs1_m = {}
    deletions_p = {}
    deletions_m = {}
    for record in SeqIO.parse(contigs, "fasta"):
        contigs1_p[record.id] = [0 for _ in range(len(record) + 1)]
        contigs1_m[record.id] = [0 for _ in range(len(record) + 1)]
        deletions_p[record.id] = set()
        deletions_m[record.id] = set()
        
    f = open(paf_path_1)
    contig_to_score = {}
    
    deletions = []
    print("Processing", paf_path_1, flush=True)
    counter = 0
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
        
        if is_primary:
            cigar_ops = re.split("([MIDNSHPX=])", get_cigar)
            if a[4] == "+":
                current_pos = int(a[2])
                
                for i in range(1, len(cigar_ops), 2):
                    get_count = int(cigar_ops[i-1])
                    get_op = cigar_ops[i]
                    
                    if get_op == "=":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 1
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "S":
                        current_pos += 1
                    elif get_op == "X":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 2
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "I":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 3
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "D":
                        deletions_p[a[0]].add(current_pos)
            else:
                current_pos = int(a[3])
                for i in range(1, len(cigar_ops), 2):
                    get_count = int(cigar_ops[i-1])
                    get_op = cigar_ops[i]
                    
                    if get_op == "=":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 1
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "S":
                        current_pos -= 1
                    elif get_op == "X":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 2
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "I":
                        for j in range(get_count):
                            if contigs1_p[a[0]][current_pos] == 0:
                                contigs1_p[a[0]][current_pos] = 3
                            else:
                                contigs1_p[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "D":
                        deletions_p[a[0]].add(current_pos - 1)
                    
        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score:
                contig_to_score[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score[a[0]]
                contig_to_score[a[0]] = (cnm + current_nm, cnm + current_len)
        
        counter += 1
        if counter % 100 == 0:
            print("Evaluated %d alignments." % counter)
    
    f = open(paf_path_2)
    contig_to_score_m = {}
    print("Processing", paf_path_2, flush=True)
    counter = 0
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
        
        if is_primary:
            cigar_ops = re.split("([MIDNSHPX=])", get_cigar)
            if a[4] == "+":
                current_pos = int(a[2])
                
                for i in range(1, len(cigar_ops), 2):
                    get_count = int(cigar_ops[i-1])
                    get_op = cigar_ops[i]
                    
                    if get_op == "=":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 1
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "S":
                        current_pos += 1
                    elif get_op == "X":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 2
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "I":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 3
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos += 1
                    elif get_op == "D":
                        deletions_m[a[0]].add(current_pos)
            else:
                current_pos = int(a[3])
                for i in range(1, len(cigar_ops), 2):
                    get_count = int(cigar_ops[i-1])
                    get_op = cigar_ops[i]
                    
                    if get_op == "=":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 1
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "S":
                        current_pos += 1
                    elif get_op == "X":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 2
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "I":
                        for j in range(get_count):
                            if contigs1_m[a[0]][current_pos] == 0:
                                contigs1_m[a[0]][current_pos] = 3
                            else:
                                contigs1_m[a[0]][current_pos] = 5
                            current_pos -= 1
                    elif get_op == "D":
                        deletions_m[a[0]].add(current_pos - 1)
        if is_primary and current_nm is not None:
            if a[0] not in contig_to_score_m:
                contig_to_score_m[a[0]] = (current_nm, current_len)
            else:
                cnm, clen = contig_to_score_m[a[0]]
                contig_to_score_m[a[0]] = (cnm + current_nm, cnm + current_len)
        
        
        counter += 1
        if counter % 100 == 0:
            print("Evaluated %d alignments." % counter)
    
    sum_nm = 0
    sum_len = 0
    for ctg in contig_to_score:
        if ctg in contig_to_score_m:
            if contig_to_score[ctg][0] / contig_to_score[ctg][1] < contig_to_score_m[ctg][0] / contig_to_score_m[ctg][1]:
                sum_nm += contig_to_score[ctg][0]
                sum_len += contig_to_score[ctg][1]
            else:
                sum_nm += contig_to_score_m[ctg][0]
                sum_len += contig_to_score_m[ctg][1]
    
    print("Assembly accuracy is: %.5f%%" % (100*(sum_len - sum_nm) / sum_len))
    print("Assembly errors are: %.5f per 1kbp" % (1000 * sum_nm / sum_len))
    
    print("Counting non-switching errors", flush=True)
    
    count_unmapped = 0
    count_both = 0
    count_paternal = 0
    count_maternal = 0
    count_multialign = 0
    count_mismatch = 0
    count_insertion = 0
    count_different = 0
    count_deletion = 0
    count_bases = 0
    count_switch = 0
    for ctg in contigs1_m:
        count_bases += len(contigs1_m[ctg])
        state = 0
        pat = contigs1_p[ctg]
        mat = contigs1_m[ctg]
        for i in range(len(contigs1_m[ctg])):
            if pat[i] == 0 and mat[i] == 0:
                count_unmapped += 1
                state = 0
            elif pat[i] == 1 and mat[i] == 1:
                count_both += 1
                state = 0
            elif pat[i] == 1 and mat[i] != 1:
                count_paternal += 1
                if state == 2:
                    count_switch += 1
                state = 1
            elif pat[i] != 1 and mat[i] == 1:
                count_maternal += 1
                if state == 1:
                    count_switch += 1
                state = 2
            elif pat[i] == 5 or mat[i] == 5:
                count_multialign += 1
                state = 0
            elif pat[i] == 2 and mat[i] == 2:
                count_mismatch += 1
                state = 0
            elif pat[i] == 3 and mat[i] == 3:
                if state != 3:
                    count_insertion += 1
                state = 3
            else:
                count_different += 1
                state = 0
    
    for ctg in deletions_p:
        z = deletions_p[ctg].intersection(deletions_m[ctg])
        count_deletion += len(z)
    
    print("Number of bases                      : %d" % count_bases)
    print("Number of unmapped bases             : %d" % count_unmapped)
    print("Number of multiply aligned bases     : %d" % count_multialign)
    print("Number of correct bases (both)       : %d" % count_both)
    print("Number of correct bases (paternal)   : %d" % count_paternal)
    print("Number of correct bases (maternal)   : %d" % count_maternal)
    print("Number of mismatches                 : %d" % count_mismatch)
    print("Number of insertions                 : %d" % count_insertion)
    print("Number of deletions                  : %d" % count_deletion)
    print("Number of ambiguous error            : %d" % count_different)
    print("Number of switches                   : %d" % count_switch)
    
    print("Mismatch rate : %.5f per 10kbp" % (10000 * count_mismatch / (count_bases - count_unmapped - count_multialign)))
    print("Insertion rate: %.5f per 10kbp" % (10000 * count_insertion / (count_bases - count_unmapped - count_multialign)))
    print("Deletion rate : %.5f per 10kbp" % (10000 * count_deletion / (count_bases - count_unmapped - count_multialign)))
    print("Error rate    : %.5f per 10kbp" % (10000 * (count_mismatch + count_insertion + count_deletion) / (count_bases - count_unmapped - count_multialign)))
    print("Switch rate   : %.5f per 10kbp" % (10000 * count_switch / (count_bases - count_unmapped - count_multialign)), flush=True)


def error_assessment(args):
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
    
    error_workdir = os.path.join(args.workdir, "error_eval")
    if not os.path.exists(error_workdir):
        os.makedirs(error_workdir)
    
    # Run minimap2 for evaluation
    path = os.path.join(error_workdir, "asm1_p.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly1, args.paternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm20", "--eqx", "-t", str(args.threads), args.paternal, args.assembly1], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm1_m.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly1, args.maternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm20", "--eqx", "-t", str(args.threads), args.maternal, args.assembly1], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm2_p.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly2, args.paternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm20", "--eqx", "-t", str(args.threads), args.paternal, args.assembly2], check=True, stdout=f)
        f.close()
    
    path = os.path.join(error_workdir, "asm2_m.paf")
    if os.path.isfile(path):
        print("Opening file from", path, flush=True)
    else:
        f = open(path, "w")
        print("Running minimap2: %s to %s" % (args.assembly2, args.maternal))
        minimap2_run = subprocess.run(["minimap2", "-cx", "asm20", "--eqx", "-t", str(args.threads), args.maternal, args.assembly2], check=True, stdout=f)
        f.close()
        
    evaluate_error(args.assembly1, os.path.join(error_workdir, "asm1_p.paf"), os.path.join(error_workdir, "asm1_m.paf"))
    evaluate_error(args.assembly2, os.path.join(error_workdir, "asm2_p.paf"), os.path.join(error_workdir, "asm2_m.paf"))
    
def haplotype_assessment():
    pass
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", type=str, help="the evaluation script to execute")
    parser.add_argument("-1", "--assembly1", type=str, help="First haplotype input")
    parser.add_argument("-2", "--assembly2", type=str, help="Second haplotype input")
    parser.add_argument("-p", "--paternal", type=str, help="First haplotype reference")
    parser.add_argument("-m", "--maternal", type=str, help="Second haplotype reference")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to run", default=1)
    parser.add_argument("-w", "--workdir", type=str, help="Number of threads to run", default="hypo_eval_wd")
    args = parser.parse_args()
    
    if args.command == "diploid_error":
        error_assessment(args)
    else:
        print("Command %s is unknown" % args.command)
        exit(1)
  
  
if __name__ == '__main__':
    main()

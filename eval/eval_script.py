import subprocess
import argparse
import os

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
        
    f = open(os.path.join(error_workdir, "asm1_p.paf"), "w")
    minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "-t", str(args.threads), args.paternal, args.assembly1], check=True, stdout=f)
    f.close()
    
    f = open(os.path.join(error_workdir, "asm1_m.paf"), "w")
    minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "-t", str(args.threads), args.maternal, args.assembly1], check=True, stdout=f)
    f.close()
    
    f = open(os.path.join(error_workdir, "asm2_p.paf"), "w")
    minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "-t", str(args.threads), args.paternal, args.assembly2], check=True, stdout=f)
    f.close()
    
    f = open(os.path.join(error_workdir, "asm2_m.paf"), "w")
    minimap2_run = subprocess.run(["minimap2", "-cx", "asm5", "-t", str(args.threads), args.maternal, args.assembly2], check=True, stdout=f)
    f.close()
    
    # Scan files for evaluation
    
    f = open(os.path.join(error_workdir, "asm1_p.paf"))
    contig_to_score = {}
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
        if is_primary and current_nm is not None:
            contig_to_score[a[0]] = (current_nm, current_len)
    
    f = open(os.path.join(error_workdir, "asm1_m.paf"))
    contig_to_score_m = {}
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
        if is_primary and current_nm is not None:
            contig_to_score_m[a[0]] = (current_nm, current_len)
    
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
    
    print("Assembly 1 accuracy is: %.5f%%" % (100*(sum_len - sum_nm) / sum_len))
    
    f = open(os.path.join(error_workdir, "asm2_p.paf"))
    contig_to_score = {}
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
        if is_primary and current_nm is not None:
            contig_to_score[a[0]] = (current_nm, current_len)
    
    f = open(os.path.join(error_workdir, "asm2_m.paf"))
    contig_to_score_m = {}
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
        if is_primary and current_nm is not None:
            contig_to_score_m[a[0]] = (current_nm, current_len)
    
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
    
    print("Assembly 2 accuracy is: %.5f%%" % (100*(sum_len - sum_nm) / sum_len))

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

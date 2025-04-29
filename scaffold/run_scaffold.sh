#!/bin/bash

set -e -o pipefail

usage (){
echo "  -k <solids>    the list of solid kmers in bitvector format                                     [          required ]"
echo "  -i <contigs>   the contigs to join in fasta format                                         [          required ]"
echo "  -I <contigs>   the 2nd haplotype contigs to join in fasta format                                         [          required ]"
echo "  -l <reads>     the long reads to map as verification                                       [          required ]"
echo "  -t <threads>   the number of threads to use.                                               [ default:        1 ]"
echo "  -o <prefix>    prefix for outputs                                                          [ default:   joined ]"
echo "  -T <tempdir>   directory to store intermediate files                                       [ default:    temp/ ]"
echo "  -f             toggle removing solid kmers with > 2 occurences                             [ default: disabled ]"
echo "  -p             toggle assuming reads as pacbio instead of ONT for mapping                  [ default: disabled ]"
echo "  -m <sortmem>   memory to use in each thread of samtools sort                               [ default:       1G ]"
echo "  -@ <sortthreads>    the threads used to sort alignment files                               [ default:       10 ]"
echo "  -s <stage>     the stage to start on. Assumed temporary files exists.                      [ default:        0 ]"
echo "  -D             debug mode                                                                  [ default: disabled ]"
echo "  -h             display this help and exit"
exit 1
}

solids=""
contigs=""
threads="1"
filter="0"
tempdir="temp"
prefix="scaffold"
longreads=""
kmerlen="17"
readtype="ont"
sortmem="1G"
sortthreads="10"
debugmode=""
stage="0"
while getopts ":k:i:I:t:fT:o:l:pm:@:s:Dh" opt; do
  case $opt in
    k)
        solids="$OPTARG"
        ;;
    i)
        contigs="$OPTARG"
        ;;
    I)
        contigs2="$OPTARG"
        ;;
    t)
        threads="$OPTARG"
        ;;
    f)
        filter="1"
        ;;
    T)
        tempdir="$OPTARG"
        ;;
    o)
        prefix="$OPTARG"
        ;;
    l)
        longreads="$OPTARG"
        ;;
    p)
        readtype="pb"
        ;;
    m)
        sortmem="$OPTARG"
        ;;
    @)
        sortthreads="$OPTARG"
        ;;
    s)
        stage="$OPTARG"
        ;;
    D)
        debugmode="1"
        ;;
    h)
        usage
        ;;
    \?) 
        echo "Invalid option -$OPTARG" >&2
        ;;
  esac
done

if [ "$solids" == "" ]; then
    echo "Option -k <solids> needed."
    exit 1
else
    echo "Reading solid kmers from $solids"
fi
if [ "$contigs" == "" ]; then
    echo "Option -i <contigs> needed."
    exit 1
else
    echo "Reading contigs from $contigs"
fi
if [ "$contigs2" == "" ]; then
    echo "Option -I <contigs> needed."
    exit 1
else
    echo "Reading contigs from $contigs"
fi
if [ "$longreads" == "" ]; then
    echo "Option -l <long reads> needed."
    exit 1
else
    echo "Reading long reads from $longreads"
fi
echo "Using read type $readtype for mappings."
echo "Using $threads threads"
if [ "$filter" == "0" ]; then
    echo "Not filtering > 2 occurences solid kmers."
else
    echo "Filtering out > 2 occurences solid kmers."
fi
echo "Using $tempdir as temporary directory."
mkdir -p $tempdir
echo "Using $sortmem memory on samtools sort."
echo "Output: $prefix.fa"

if (( $(echo "$stage <= 0" | bc -l) )); then
    echo "[SCAFFOLD: STEP 0] Filtering reads"
    minimap2 -I 64G -ax map-$readtype -t $threads $contigs $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map_initial.sorted.bam
    minimap2 -I 64G -ax map-$readtype -t $threads $contigs2 $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map_initial2.sorted.bam

    python get_necessary_reads.py $tempdir/map_initial.sorted.bam $tempdir/map_initial2.sorted.bam $tempdir/filtered_reads.fa
fi

if (( $(echo "$stage <= 1" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 1] Finding scaffolds"
        echo "./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter"
        ./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/filtered_reads.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter
    else
        touch $tempdir/runscaffold.log
        echo "[SCAFFOLD: STEP 1] Finding scaffolds <DEBUG>"
        echo "./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt"
        ./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/filtered_reads.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt 2>&1 | tee -a $tempdir/runscaffold.log
    fi
fi

if (( $(echo "$stage <= 2" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 2] Joining scaffolds"
        echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
        python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt
    else
        echo "[SCAFFOLD: STEP 2] Joining scaffolds <DEBUG>"
        echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
        python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl $tempdir/debug_2.txt > $tempdir/identity.txt
    fi
fi

if (( $(echo "$stage <= 3.1" | bc -l) )); then
    echo "[SCAFFOLD: STEP 3.1] Mapping long reads"
    echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam
    echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam
    echo "samtools index -@ $threads $tempdir/map.sorted.bam"
    samtools index -@ $threads $tempdir/map.sorted.bam
    echo "samtools index -@ $threads $tempdir/map2.sorted.bam"
    samtools index -@ $threads $tempdir/map2.sorted.bam
fi

if (( $(echo "$stage <= 3.2" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 3.2] Filtering scaffolds"
        echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2"
        python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2
    else
        echo "[SCAFFOLD: STEP 3.2] Filtering scaffolds"
        echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2 $tempdir/debug_3.txt"
        python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2 $tempdir/debug_3.txt 2>&1 | tee -a $tempdir/runscaffold.log
    fi
fi

if (( $(echo "$stage <= 3.5" | bc -l) )); then
    echo "[SCAFFOLD: STEP 3.5] Removing duplicates"
    minimap2 -ax asm5 -t $threads $tempdir/iter1_1.fa $tempdir/iter1_1.fa | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/all_vs_all.bam
    minimap2 -ax asm5 -t $threads $tempdir/iter1_2.fa $tempdir/iter1_2.fa | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/all_vs_all2.bam
    mv $tempdir/iter1_1.fa $tempdir/iter1_1_before_duplicate.fa
    mv $tempdir/iter1_2.fa $tempdir/iter1_2_before_duplicate.fa
    python remove_duplicates.py $tempdir/iter1_1_before_duplicate.fa $tempdir/all_vs_all.bam $tempdir/iter1_1.fa
    python remove_duplicates.py $tempdir/iter1_2_before_duplicate.fa $tempdir/all_vs_all.bam $tempdir/iter1_2.fa
fi

mkdir -p $tempdir/iter2
oldtemp=$tempdir
tempdir=$tempdir/iter2

if (( $(echo "$stage <= 4" | bc -l) )); then
    echo "[SCAFFOLD: STEP 4] Creating new contigs"
    echo "./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter"
fi

if (( $(echo "$stage <= 5" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 5] Finding scaffolds"
        echo "./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter"
        ./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $longreads $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter
    else
        echo "[SCAFFOLD: STEP 5] Finding scaffolds <DEBUG>"
        echo "./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt"
        ./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $longreads $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt 2>&1 | tee -a $tempdir/runscaffold.log
    fi
fi

if (( $(echo "$stage <= 6" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 6] Joining scaffolds"
        echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
        python join_scaffold.py $oldtemp/iter1_1.fa $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $oldtemp/iter1_2.fa $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt
    else
        echo "[SCAFFOLD: STEP 6] Joining scaffolds <DEBUG>"
        echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
        python join_scaffold.py $oldtemp/iter1_1.fa $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $oldtemp/iter1_2.fa $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl $tempdir/debug_2.txt> $tempdir/identity.txt
    fi
fi

if (( $(echo "$stage <= 7.1" | bc -l) )); then
    echo "[SCAFFOLD: STEP 7.1] Mapping long reads"
    echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam
    echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam
    echo "samtools index -@ $threads $tempdir/map.sorted.bam"
    samtools index -@ $threads $tempdir/map.sorted.bam
    echo "samtools index -@ $threads $tempdir/map2.sorted.bam"
    samtools index -@ $threads $tempdir/map2.sorted.bam
fi

if (( $(echo "$stage <= 7.2" | bc -l) )); then
    if [ "$debugmode" == "" ]; then
        echo "[SCAFFOLD: STEP 7.2] Filtering scaffolds"
        echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2"
        python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter2_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter2_2
    else
        echo "[SCAFFOLD: STEP 7.2] Filtering scaffolds"
        echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2 debug_3.txt"
        python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter2_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter2_2 $tempdir/debug_3.txt 2>&1 | tee -a $tempdir/runscaffold.log
    fi
fi

if (( $(echo "$stage <= 7.5" | bc -l) )); then
    echo "[SCAFFOLD: STEP 7.5] Removing duplicates"
    minimap2 -ax asm5 -t $threads $tempdir/iter2_1.fa $tempdir/iter2_1.fa | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/all_vs_all.bam
    minimap2 -ax asm5 -t $threads $tempdir/iter2_2.fa $tempdir/iter2_2.fa | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/all_vs_all2.bam
    mv $tempdir/iter2_1.fa $tempdir/iter2_1_before_duplicate.fa
    mv $tempdir/iter2_2.fa $tempdir/iter2_2_before_duplicate.fa
    python remove_duplicates.py $tempdir/iter2_1_before_duplicate.fa $tempdir/all_vs_all.bam ${prefix}_1.fa 1
    python remove_duplicates.py $tempdir/iter2_2_before_duplicate.fa $tempdir/all_vs_all.bam ${prefix}_2.fa 1
fi

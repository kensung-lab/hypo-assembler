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
while getopts ":k:i:I:t:fT:o:l:pm:@:Dh" opt; do
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

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 1] Finding scaffolds"
    echo "./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter"
    ./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter
else
    touch $tempdir/runscaffold.log
    echo "[SCAFFOLD: STEP 1] Finding scaffolds <DEBUG>"
    echo "./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt"
    ./find_scaffold $kmerlen $solids $contigs $contigs2 $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt 2>&1 | tee -a $tempdir/runscaffold.log
fi

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 2] Joining scaffolds"
    echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
    python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt
else
    echo "[SCAFFOLD: STEP 2] Joining scaffolds <DEBUG>"
    echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
    python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl $tempdir/debug_2.txt > $tempdir/identity.txt
fi

echo "[SCAFFOLD: STEP 3] Mapping long reads"
echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam"
minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam
echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam"
minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam
echo "samtools index -@ $threads $tempdir/map.sorted.bam"
samtools index -@ $threads $tempdir/map.sorted.bam
echo "samtools index -@ $threads $tempdir/map2.sorted.bam"
samtools index -@ $threads $tempdir/map2.sorted.bam

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 4] Finalization"
    echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2"
    python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2
else
    echo "[SCAFFOLD: STEP 4] Finalization"
    echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2 $tempdir/debug_3.txt"
    python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam $tempdir/iter1_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam $tempdir/iter1_2 $tempdir/debug_3.txt 2>&1 | tee -a $tempdir/runscaffold.log
fi

mkdir -p $tempdir/iter2
oldtemp=$tempdir
tempdir=$tempdir/iter2

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 1-2] Finding scaffolds"
    echo "./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter"
    ./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter
else
    echo "[SCAFFOLD: STEP 1-2] Finding scaffolds <DEBUG>"
    echo "./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt"
    ./find_scaffold $kmerlen $solids $oldtemp/iter1_1.fa $oldtemp/iter1_2.fa $tempdir/scaffold.txt $tempdir/scaffold2.txt $threads $filter $tempdir/debug.txt 2>&1 | tee -a $tempdir/runscaffold.log
fi

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 2-2] Joining scaffolds"
    echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
    python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt
else
    echo "[SCAFFOLD: STEP 2-2] Joining scaffolds <DEBUG>"
    echo "python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl > $tempdir/identity.txt"
    python join_scaffold.py $contigs $tempdir/scaffold.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir $contigs2 $tempdir/scaffold2.txt $tempdir/intermediate2.fa $tempdir/obj2.pkl $tempdir/debug_2.txt> $tempdir/identity.txt
fi

echo "[SCAFFOLD: STEP 3-2] Mapping long reads"
echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam"
minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map.sorted.bam
echo "minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam"
minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate2.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/map2.sorted.bam
echo "samtools index -@ $threads $tempdir/map.sorted.bam"
samtools index -@ $threads $tempdir/map.sorted.bam
echo "samtools index -@ $threads $tempdir/map2.sorted.bam"
samtools index -@ $threads $tempdir/map2.sorted.bam

if [ "$debugmode" == "" ]; then
    echo "[SCAFFOLD: STEP 4-2] Finalization"
    echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2"
    python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2
else
    echo "[SCAFFOLD: STEP 4-2] Finalization"
    echo "python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2 debug_3.txt"
    python filter_scaffold.py $tempdir/obj.pkl $tempdir/map.sorted.bam ${prefix}_1 $tempdir/obj2.pkl $tempdir/map2.sorted.bam ${prefix}_2 $tempdir/debug_3.txt 2>&1 | tee -a $tempdir/runscaffold.log
fi

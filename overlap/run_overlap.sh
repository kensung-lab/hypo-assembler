#!/bin/bash

usage (){
echo "  -k <solids>    the list of solid kmers in bitvector format                                     [          required ]"
echo "  -i <contigs>   the contigs to join in fasta format                                         [          required ]"
echo "  -l <reads>     the long reads to map as verification                                       [          required ]"
echo "  -b <aligns>    long reads alignment to contigs, to filter reads                            [ default:     none ]"
echo "  -t <threads>   the number of threads to use.                                               [ default:        1 ]"
echo "  -o <prefix>    prefix for outputs                                                          [ default:   overlap ]"
echo "  -T <tempdir>   directory to store intermediate files                                       [ default:    temp/ ]"
echo "  -f             toggle removing solid kmers with > 2 occurences                             [ default: disabled ]"
echo "  -m <sortmem>   memory to use in each thread of samtools sort                               [ default:       1G ]"
echo "  -h             display this help and exit"
exit 1
}

solids=""
contigs=""
threads="1"
filter="0"
tempdir="temp"
prefix="overlap"
longreads=""
readtype="ont"
sortmem="1G"
kmerlen="17"
aligns=""
while getopts ":k:i:t:fb:T:o:l:pm:h" opt; do
  case $opt in
    k)
        solids="$OPTARG"
        ;;
    i)
        contigs="$OPTARG"
        ;;
    t)
        threads="$OPTARG"
        ;;
    f)
        filter="1"
        ;;
    b)
        aligns="$OPTARG"
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

echo "[OVERLAP: STEP 1] Finding overlaps"
echo "./find_overlap $solids $contigs $threads $filter > $tempdir/overlap.txt"
./find_overlap $kmerlen $solids $contigs $threads $filter > $tempdir/overlap.txt

echo "[OVERLAP: STEP 2] Joining overlaps"
echo "python join.py $contigs $tempdir/overlap.txt $tempdir/intermediate.fa $tempdir/obj.pkl"
python join_overlap.py $contigs $tempdir/overlap.txt $tempdir/intermediate.fa $tempdir/obj.pkl $tempdir > $tempdir/identity.txt

if [ "$aligns" == "" ]; then
    echo "Alignments not provided; Skipping Step 3: Filtering long reads to remap"
    echo "[OVERLAP: STEP 4] Mapping long reads"
    echo "minimap2 -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS > $tempdir/map.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $longreads | samtools view -bS > $tempdir/map.bam
    echo "samtools sort -@ $threads -o $tempdir/map.sorted.bam $tempdir/map.bam"
    samtools sort -@ $threads -m $sortmem -o $tempdir/map.sorted.bam $tempdir/map.bam
    echo "samtools index -@ $threads $tempdir/map.sorted.bam"
    samtools index -@ $threads $tempdir/map.sorted.bam
else
    echo "[OVERLAP: STEP 3] Filtering long reads by previous alignment"
    python filter_alignments.py $aligns $tempdir/filtered_reads.fa
    echo "[OVERLAP: STEP 4] Mapping long reads"
    echo "minimap2 -ax map-$readtype -t $threads $tempdir/intermediate.fa $tempdir/filtered_reads.fa | samtools view -bS > $tempdir/map.bam"
    minimap2 -I 64G -ax map-$readtype -t $threads $tempdir/intermediate.fa $tempdir/filtered_reads.fa | samtools view -bS > $tempdir/map.bam
    echo "samtools sort -@ $threads -o $tempdir/map.sorted.bam $tempdir/map.bam"
    samtools sort -@ $threads -m $sortmem -o $tempdir/map.sorted.bam $tempdir/map.bam
    echo "samtools index -@ $threads $tempdir/map.sorted.bam"
    samtools index -@ $threads $tempdir/map.sorted.bam
fi

echo "[OVERLAP: STEP 5] Finalization"
echo "python filter.py $tempdir/obj.pkl $tempdir/map.sorted.bam $prefix"
python filter_overlap.py $tempdir/obj.pkl $tempdir/map.sorted.bam $prefix

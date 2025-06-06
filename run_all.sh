#!/bin/bash

set -e -o pipefail

usage (){
echo "  -1 <reads_1>        short read pair 1                               [          required ]"
echo "  -2 <reads_2>        short read pair 2                               [          required ]"
echo "  -l <long_reads>     long reads file                                 [          required ]"
echo "  -d <draft>          draft assembly                                  [ default:      none]"
echo "  -o <output_prefix>  prefix of the output files                      [ default:      hypo]"  
echo "  -B <long_read_map>  mapping of long reads to draft                  [ default:      none]"
echo "  -t <threads>        the number of threads to use.                   [ default:        1 ]"
echo "  -m <sortmem>        the memory used to sort alignment files         [ default:       1G ]"
echo "  -@ <sortthreads>    the threads used to sort alignment files        [ default:       10 ]"
echo "  -k <kmer length>    the length of the solid kmer used               [ default:       17 ]"
echo "  -s <estimated_size> estimated genome size (suffix K/M/G accepted)   [ default:       3G ]"
echo "  -M <kmcmem>         max memory to be used by KMC (in GB)            [ default:       12 ]"
echo "  -T <tempdir>   directory to store intermediate files                [ default:    temp/ ]"
echo "  -D             debug mode                                           [ default: disabled ]"
echo "  -h             display this help and exit"
exit 1
}

reads1=""
reads2=""
longreads=""
longbam=""
draft=""
threads="1"
sortmem="1G"
kmerlen="17"
tempdir="temp/"
genomesize="3G"
outputpref="hypo"
sortthreads="10"
kmcmem="12"
debugmode=""
while getopts "1:2:l:d:B:t:T:hs:o:m:k:@:M:D" opt; do
  case $opt in
    1)
        reads1="$OPTARG"
        ;;
    2)
        reads2="$OPTARG"
        ;;
    l)
        longreads="$OPTARG"
        ;;
    d)
        draft="$OPTARG"
        ;;
    B)
        longbam="$OPTARG"
        ;;
    t)
        threads="$OPTARG"
        ;;
    T)
        tempdir="$OPTARG"
        ;;
    o)
        outputpref="$OPTARG"
        ;;
    s)
        genomesize="$OPTARG"
        ;;
    m)
        sortmem="$OPTARG"
        ;;
    M)
        kmcmem="$OPTARG"
        ;;
    k)
        kmerlen="$OPTARG"
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

if [ "$reads1" == "" ]; then
    echo "Option -1 <reads_1> needed."
    exit 1
else
    echo "Short read file 1: $reads1"
fi

if [ "$reads2" == "" ]; then
    echo "Option -2 <reads_1> needed."
    exit 1
else
    echo "Short read file 2: $reads2"
fi

if [ "$longreads" == "" ]; then
    echo "Option -l <long reads> needed."
    exit 1
else
    echo "Long reads: $longreads"
fi

echo "Using $threads threads"
echo "Using $tempdir as temporary directory."
mkdir -p $tempdir

if ! [ -x "$(command -v minimap2)" ]; then
    echo 'Error: minimap2 is not installed.'
    exit 1
fi

if ! [ -x "$(command -v samtools)" ]; then
    echo 'Error: samtools is not installed.'
    exit 1
fi

if ! [ -x "$(command -v kmc)" ]; then
    echo 'Error: KMC is not installed.'
    exit 1
fi

touch $tempdir/run.log

if [ "$draft" == "" ]; then
    echo "Initial draft is not provided. Running FylE." | tee -a $tempdir/run.log
    
    if ! [ -x "$(command -v flye)" ]; then
        echo 'Error: flye is not installed.' | tee -a $tempdir/run.log
        exit 1
    else
        flye --nano-raw $longreads --genome-size $genomesize --threads $threads --out-dir $tempdir/flye/ 2>&1 | tee $tempdir/flye.log
        draft=$tempdir/flye/assembly.fasta
    fi
else
    echo "Using initial draft: $draft"
fi

if [ "$longbam" == "" ]; then
    echo "Mapping long reads to draft" | tee -a $tempdir/run.log
    minimap2 -ax map-ont -t $threads $draft $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/long_align.bam
    longbam=$tempdir/long_align.bam
else
    echo "Long read mapping: $longbam" | tee -a $tempdir/run.log
fi

echo "[STEP 1] Getting solid kmers" | tee -a $tempdir/run.log
echo $reads1 > $tempdir/shorts.txt
echo $reads2 >> $tempdir/shorts.txt
./suk -k $kmerlen -i @"$tempdir"/shorts.txt -t $threads -m $kmcmem -e -w $tempdir/suk_kmc -o $tempdir/SUK 2>&1 | tee $tempdir/suk.log
# mv SUK_k17.bv $tempdir/SUK_k17.bv

echo "[STEP 2] Scanning misjoin" | tee -a $tempdir/run.log
python scan_misjoin.py $draft $longbam $tempdir/misjoin.fa 2>&1 | tee -a $tempdir/misjoin.log

echo "[STEP 3] Finding overlaps" | tee -a $tempdir/run.log
./run_overlap.sh -k $tempdir/SUK_k$kmerlen.bv -i $tempdir/misjoin.fa -l $longreads -t $threads -o $tempdir/overlap -T $tempdir/overlap_temp 2>&1 | tee -a $tempdir/overlap.log

echo "[STEP 4] Realignment for polishing" | tee -a $tempdir/run.log
minimap2 -I 64G -ax map-ont -t $threads $tempdir/overlap.fa $longreads | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/overlap_long.bam
minimap2 -I 64G -ax sr -t $threads $tempdir/overlap.fa $reads1 $reads2 | samtools view -bS | samtools sort -@ $sortthreads -m $sortmem -o $tempdir/overlap_short.bam

echo "[STEP 5] Polishing" | tee -a $tempdir/run.log
if [ "$debugmode" == "" ]; then
    ./hypo -d $tempdir/overlap.fa -s $genomesize -B $tempdir/overlap_long.bam -C 60 -b $tempdir/overlap_short.bam -r @"$tempdir"/shorts.txt -c 100 -L $kmcmem -t $threads -p 1 -o $tempdir/polished 2>&1 | tee $tempdir/polish.log
else
    ./hypo -d $tempdir/overlap.fa -s $genomesize -B $tempdir/overlap_long.bam -C 60 -b $tempdir/overlap_short.bam -r @"$tempdir"/shorts.txt -c 100 -L $kmcmem -t $threads -p 1 -o $tempdir/polished -i 2>&1 | tee $tempdir/polish.log
fi

echo "[STEP 6] Scaffolding" | tee -a $tempdir/run.log
if [ "$debugmode" == "" ]; then
    ./run_scaffold.sh -k $tempdir/SUK_k$kmerlen.bv -i $tempdir/polished_1.fa -I $tempdir/polished_2.fa -l $longreads -t $threads -o $tempdir/scaffold -T $tempdir 2>&1 | tee $tempdir/scaffold.log
else
    ./run_scaffold.sh -k $tempdir/SUK_k$kmerlen.bv -i $tempdir/polished_1.fa -I $tempdir/polished_2.fa -l $longreads -t $threads -o $tempdir/scaffold -T $tempdir -D 2>&1 | tee $tempdir/scaffold.log
fi
cp $tempdir/scaffold_1.fa ${outputpref}_1.fa
cp $tempdir/scaffold_2.fa ${outputpref}_2.fa

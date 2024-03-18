# Hypo-Assembler: A diploid genome polisher and assembler
=======================================================================

Hypo is a package of both polisher, capable of correcting draft genomes, and an assembler to assemble a diploid genome. It exploits unique genomic kmers to both selectively polish only segments of contigs and to increase mapping accuracy.

## Installation
Hypo is only available for Unix-like platforms (Linux and MAC OS). Currently, we provide the option of installation from source with CMake.
CmakeLists is provided in the project root folder. 

### Pre-requisites
The following requirements are assumed to be installed (with path to their binaries available in $PATH).
- Zlib
- OpenMP
- GCC (>=7.3)
- [KMC3](https://github.com/refresh-bio/KMC)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [FlyE](https://github.com/fenderglass/Flye) (optional, if draft assembly is not provided)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [BioPython](https://biopython.org/wiki/Download)

### Building

First, you need to clone the repository recusrively (for the submodules)
```
git clone --recursive https://github.com/kensung-lab/hypo-assembler
```

After cloning, the easiest way to run the pipeline is to run build_all.sh. This will write a run_all directory, from where the run can be easily executed.
```
./build_all.sh
cd run_all
```

The following optional options are available for build_all.sh
```console
 Usage: ./build_all.sh <args>
 
** Optional parameters
    -t <thread count>
    The number of threads to run make
    [Default] 10
    
    -n
    If present, will run optimized for native (i.e. -march=native)
    
    -o <build directory>
    Directory to put executables and scripts.
    [Default] run_all
```


### Usage

Usage with run_all.sh is as follows:

```console
 Usage: ./run_all.sh <args>

** Mandatory arguments
    -1 <short read file 1>
    The file containing the first of the short read input
    
    -2 <short read file 2>
    The file containing the pair of the short read input
    
    -l <long reads file>
    The file containing long reads

** Optional parameters
    -t <thread count>
    The number of threads to run
    [Default] 1
    
    -o <output prefix>
    The prefix for the output file names. The outputs will be <prefix>_1.fa and <prefix>_2.fa
    [Default] hypo
    
    -d <assembly draft>
    The file containing the original draft. If not provided, FlyE will be run from the given long reads
    [Default] None
    
    -B <long reads mapping to draft>
    Mappings of the long reads to the original draft as provided in -d. If not provided, minimap2 will be run to align the long reads
    [Default] None
    
    -s <estimated size>
    Estimated size of the assembled genome. K/M/G suffix accepted
    [Default] 3G
    
    -T <tempdir>
    Directory to write temporary files
    [Default] temp
    
    -m <sortmem>
    The memory used for each samtools sorting thread. Samtools will use number of threads * sortmem memory.
    [Default] 1G
    
    -@ <sortthreads>
    The threads used for each samtools sorting. Separate from the normal number of threads.
    [Default] 10
    
    -k <kmer length>
    The length of solid kmers used in the method. Recommended is 17 for human genome.
    [Default] 17
```

### Demo

We have included a small demo data to test the installation in demo/

For the puposes of this demo, the following external library is assumed:
- [minimap2 2.26](https://github.com/lh3/minimap2/releases/tag/v2.26)
- [FlyE 2.92](https://github.com/fenderglass/Flye/releases/tag/2.9.2)
- [KMC 3.2.2](https://github.com/refresh-bio/KMC/releases/tag/v3.2.2)

Assuming build_all.sh is already run:
```
cd run_all
./run_all.sh -1 ../demo/il1.fq -2 ../demo/il2.fq -l ../demo/ont.fq.gz -t 40
```

This will run hypo-assembler on 40 threads.

The outputs will be hypo_1.fa and hypo_2.fa in the run_all directory.

To make sure the program runs properly, you can compare the output with the results in demo/expected_result/hypo_1.fa and demo/expected_result/hypo_2.fa, i.e.
```
diff hypo_1.fa ../demo/expected_result/hypo_1.fa
diff hypo_2.fa ../demo/expected_result/hypo_2.fa
```

If both commands execute without any output, then the program works properly.

### Running hypo-polisher

We recommend running the full package of hypo-assembler with -d option to improve a quality of a draft genome. However, we also provide the option of running the polisher separate from the assembler.

To run them, simply run the binary "hypo" built from the polisher directory.
```
Usage: hypo <args>

 ** Mandatory args:
	-r, --reads-short <str>
	Input file name containing reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.

	-d, --draft <str>
	Input file name containing the draft contigs (in fasta/fastq format; can be compressed). 

	-b, --bam-sr <str>
	Input file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). 

	-c, --coverage-short <int>
	Approximate mean coverage of the short reads. 

	-s, --size-ref <str>
	Approximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). 

** Optional args:
	-B, --bam-lr <str>
	Input file name containing the alignments of long reads against the draft (in bam/sam format; must have CIGAR information). 
	[Only Short reads polishing will be performed if this argument is not given]

	-o, --output <str>
	Output file name. 
	[Default] hypo_<draft_file_name>.fasta in the working directory.

 	-t, --threads <int>
	Number of threads. 
	[Default] 1.

 	-p, --processing-size <int>
	Number of contigs to be processed in one batch. Lower value means less memory usage but slower speed. 
	[Default] All the contigs in the draft.

 	-m, --match-sr <int>
	Score for matching bases for short reads. 
	[Default] 5.

 	-x, --mismatch-sr <int>
	Score for mismatching bases for short reads. 
	[Default] -4.

 	-g, --gap-sr <int>
	Gap penalty for short reads (must be negative). 
	[Default] -8.

 	-M, --match-lr <int>
	Score for matching bases for long reads. 
	[Default] 3.

 	-X, --mismatch-lr <int>
	Score for mismatching bases for long reads. 
	[Default] -5.

 	-G, --gap-lr <int>
	Gap penalty for long reads (must be negative). 
	[Default] -4.

 	-n, --ned-th <int>
	Threshold for Normalised Edit Distance of long arms allowed in a window (in %). Higher number means more arms allowed which may slow down the execution.
	[Default] 20.

 	-q, --qual-map-th <int>
	Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. 
	[Default] 2.

 	-i, --intermed
	Store or use (if already exist) the intermediate files. 
	[Currently, only Solid kmers are stored as an intermediate file.].

 	-h, --help
	Print the usage. 
```

### Contact

Other than raising issues in github, you can contact joshuac@comp.nus.edu.sg for specific issues.

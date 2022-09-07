# Hypo-Assembler: A diploid genome polisher and assembler
=======================================================================

Hypo is a package of both polisher, capable of correcting draft genomes, and an assembler to assemble a diploid genome. It exploits unique genomic kmers to both selectively polish only segments of contigs and to increase mapping accuracy.

Hypo provided three different modules:
+ Hypo-short: A paired-end short-read (Illumina) based polisher. Its output is a single corrected genome from a given draft genome given Illumina reads only.
+ Hypo-hybrid: A paired-end short-read (Illumina) and ONT-reads based polisher. Its output, when prompted, is two genomes representing different haplotypes.
+ Hypo-assembler: An assembler based on paired-end short read (Illumina) and ONT reads. When prompted, it outputs two genomes, one for each haplotype.

## Installation
Hypo is only available for Unix-like platforms (Linux and MAC OS). Currently, we provide only the option of installation from source with CMake.

CmakeLists is provided in the project root folder. 

### Pre-requisites
For installing from the source, the following requirements are assumed to be installed already (with path to their binaries available in $PATH).
- Zlib
- OpenMP
- GCC (>=7.3)
  * Following are the commands to update GCC (say to GCC 8) on an Ubuntu machine (from say GCC 5):
  ```console
    sudo apt-get update; sudo apt-get install build-essential software-properties-common -y;
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y; sudo apt update; 
    sudo apt install gcc-snapshot -y; sudo apt update
    sudo apt install gcc-8 g++-8 -y; 
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 --slave /usr/bin/g++ g++ /usr/bin/g++-8
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5;
  ```
- [HTSLIB](https://github.com/samtools/htslib) (version >=1.10)
  + If htslib version 1.10 or higher is not installed, we recommend using `install_htslib.sh` in the project folder to install it locally.

- [KMC3](https://github.com/refresh-bio/KMC)
  + Can also be installed using conda as follows: 
  ```console
  conda install -c bioconda kmc
  ``` 

### Building Executable
Run the following commands to build a binary (executable) `hypo` in `build/bin` :
If htslib version 1.10 or higher is installed:
```console
  git clone --recursive https://github.com/kensung-lab/hypo-assembler hypo-assembler
  cd hypo-assembler
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make -j 8
```
If htslib is not installed or the version is smaller than 1.10:
```console
  git clone --recursive https://github.com/kensung-lab/hypo-assembler hypo-assembler
  cd hypo
  chmod +x install_htslib.sh
  ./install_htslib.sh
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make -j 8
```
**Notes:** 
* If `--recursive` was omitted from `git clone`, please run `git submodule init` and `git submodule update` before compiling.
* If target machine is different from the one on which Hypo is being compiled, exclude the flag `-Doptimise_for_native=ON`.


## Usage of the tool: 
```console
 Usage: hypo short-polish <args>

 ** Mandatory args:
	-r, --reads-short <str>
	Input file name containing illumina reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.

	-d, --draft <str>
	Input file name containing the draft contigs (in fasta/fastq format; can be compressed). 

	-c, --coverage-short <int>
	Approximate mean coverage of the short reads. 

	-s, --size-ref <str>
	Approximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). 

** Optional args:
    -b, --bam-sr <str>
	Input file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). 
    [Default] An alignment will be done using minimap2.

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

 	-q, --qual-map-th <int>
	Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. 
	[Default] 2.

 	-i, --intermed
	Store or use (if already exist) the intermediate files. 

 	-h, --help
	Print the usage. 
```

```console
 Usage: hypo hybrid-polish <args>

 ** Mandatory args:
	-r, --reads-short <str>
	Input file name containing illumina reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.
    
    -l, --reads-long <str>
	Input file name containing long reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.

	-d, --draft <str>
	Input file name containing the draft contigs (in fasta/fastq format; can be compressed). 

	-c, --coverage-short <int>
	Approximate mean coverage of the short reads. 

	-s, --size-ref <str>
	Approximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). 

** Optional args:
    -b, --bam-sr <str>
	Input file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). 
    [Default] An alignment will be done using minimap2.
    
    -2, --diploid
    Whether to produce diploid genome or not.
    [Default] Produce one haploid genome.

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

 	-q, --qual-map-th <int>
	Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. 
	[Default] 2.

 	-i, --intermed
	Store or use (if already exist) the intermediate files. 
	[Currently, only Solid kmers are stored as an intermediate file.].

 	-h, --help
	Print the usage. 


```

```console
 Usage: hypo assemble <args>

 ** Mandatory args:
	-r, --reads-short <str>
	Input file name containing illumina reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.
    
    -l, --reads-long <str>
	Input file name containing long reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.

	-c, --coverage-short <int>
	Approximate mean coverage of the short reads. 

	-s, --size-ref <str>
	Approximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). 

** Optional args:
    -2, --diploid
    Whether to produce diploid genome or not.
    [Default] Produce one haploid genome.

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

 	-q, --qual-map-th <int>
	Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. 
	[Default] 2.

 	-i, --intermed
	Store or use (if already exist) the intermediate files. 
	[Currently, only Solid kmers are stored as an intermediate file.].

 	-h, --help
	Print the usage. 


```

### Output File
If no `--output` (or `-o`) is provided, `hypo_X.fasta` will be created in the working directory where <X> is the name of the draft file. 

### Intermediate Files
If `--intermed` (or `-i`) is used, hypo will store the intermediate files corresponding to the solid kmers in a folder named `aux` in the first run. Another file indicating the progress of the run will also be created in that folder. In the next run (with `-i`), those intermediate files will be used instead. Remove `-i` or delete `aux` folder to start Hypo from the beginning.

## External Libraries

 * [sdsl-lite](https://github.com/simongog/sdsl-lite) has been used for rank-select and bit-vectors data-structures.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.
 * [suk](https://github.com/Ritu-Kundu/suk) has been used as the module to compute the solid (unique) kmers.
 * An adapted version of [FALCON-Sense](https://github.com/marbl/canu) was used for consensus.
 * [FlyE](https://github.com/fenderglass/Flye) is used to construct intial draft for hypo-assembler.

## Contact
Other than raising issues on Github, you can contact Joshua Casey (joshuac@comp.nus.edu.sg) for getting help in installation/usage or any other related query.

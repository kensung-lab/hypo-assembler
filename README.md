# Hypo-Assembler: A diploid genome polisher and assembler
=======================================================================

Hypo is a package of both polisher, capable of correcting draft genomes, and an assembler to assemble a diploid genome. It exploits unique genomic kmers to both selectively polish only segments of contigs and to increase mapping accuracy.

## Installation
Hypo is only available for Unix-like platforms (Linux and MAC OS). Currently, we provide the option of installation from source with CMake and a Docker image.

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
- [HTSLIB](https://github.com/samtools/htslib) (version >=1.16)
  + If htslib version 1.16 or higher is not installed, we recommend using `install_htslib.sh` in the project folder to install it locally.

- [KMC3](https://github.com/refresh-bio/KMC)
- [minimap2](https://github.com/lh3/minimap2)
- [FlyE](https://github.com/fenderglass/Flye)
- [samtools](https://github.com/samtools/samtools)


To install and compile, use cmake to build i.e.

```
mkdir build
cd build
cmake ..
make
```

If there are issues on building concerning htslib, run install_htslib.sh to fix out the htslib dependency, or put a fully built htslib on external/htslib.

Docker image is available on dockerhub: https://hub.docker.com/r/jcsyd/hypo-assembler
To run it, simply pull:
```
docker pull jcsyd/hypo-assembler:v0.9
```

### Usage

```console
 Usage: hypo <args>

** Mandatory arguments
    -1 or --short-read-1
    The file containing the first of the short read input
    
    -2 or --short-read-2
    The file containing the pair of the short read input
    
    -l or --long-read
    The file containing long reads

** Optional parameters
    -3 or --hic-read-1
    The file containing the first of the hi-c read input
    [Default] None
    
    -4 or --hic-read-2
    The file containing the second of the hi-c read input
    [Default] None
    
    -t or --threads
    The number of threads to run (including FlyE and minimap2)
    [Default] 1
    
    -w or --workdir
    The working directory to put temporary files and results
    [Default] $PWD/hypo_wd
    
    -F or --flye-path
    The path to flye executable
    [Default] flye (expected to find on $PATH)
    
    -M or --minimap2-path
    The path to minimap2 executable
    [Default] minimap2 (expected to find on $PATH)
    
    -S or --samtools-path
    The path to samtools executable
    [Default] samtools (expected to find on $PATH)
    
    -s or --size-ref
    Estimated genome size
    [Default] 3000000000
    
    -c or --coverage-short
    Estimated coverage of the short reads
    [Default] Estimated from reads
    
    -C or --coverage-long
    Estimated coverage of the long reads
    [Default] Estimated from reads
    
    -@ or --samtools-thread
    The number of thread to use for samtools (will be used in parallel with the initial number of threads, (samtools -@)
    [Default] 1
    
    -Z or --samtools-memory
    The memory used by samtools (samtools -m)
    [Default] 768M
    
    -X or --samtools-temp
    The path for temporary files by samtools (samtools -T)
    [Default] hypo_wd/short_read_initial and hypo_wd/long_read_initial
    
    -I or --debug
    Turns on debug mode to print out temporary files. Not recommended due to high I/O requirements, unless there is an issue encountered.
    [Default] Off
```

### Demo

We have included a small demo data to test the installation in demo/

To run the demo, after compilation, run the following command:

```
hypo -1 demo/il1.fq -2 demo/il2.fq -l demo/nanopore.fq.gz -3 demo/hic1.fq -4 demo/hic2.fq -t 40 -s 1500000
```

This will run hypo-assembler on 40 threads. It is estimated to finish within 20 minutes.
The expected result is hypo_wd/final_1.fa and hypo_wd/final_2.fa.

To verify the result, refer to [this page](https://github.com/kensung-lab/hypo-assembler/tree/main/eval).

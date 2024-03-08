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


The easiest way to run the pipeline is to run build_all.sh. This will write a run_all directory, from where the run can be easily executed.

```
./build_all.sh
cd run_all
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

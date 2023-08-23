# Hypo-Assembler Evaluation Toolkit
=======================================================================

The evaluation toolkit is designed to accurately assess an assembly of a diploid genome.

## Pre-requisites

The toolkit is written on python3. The following python packages are needed:
- pysam (https://github.com/pysam-developers/pysam)
- BioPython (https://biopython.org/)

The following packages are also needed for the different evaluations:
- minimap2 (https://github.com/lh3/minimap2)
- YAK (https://github.com/lh3/yak)
- QUAST (https://github.com/ablab/quast)

## Usage

Evaluation is divided into several evaluation modules. Depending on the available dataset, some modules can be run.

### Error evaluation for diploid genome

The requirements for diploid genome error evaluation is a diploid reference genome, paternal and maternal.

```console
 Usage: python eval_script.py diploid_error <args>

** Mandatory arguments
    -1 or --assembly1
    The file containing the first assembled haplotype
    
    -2 or --assembly2
    The file containing the second assembled haplotype
    
    -p or --paternal
    The file containing the paternal reference
    [Note] For the evaluation purposes, maternal and paternal are interchangable and only reflected on the report
    
    -m or --maternal
    The file containing the maternal reference
    [Note] For the evaluation purposes, maternal and paternal are interchangable and only reflected on the report

** Optional parameters
    -t or --threads
    The number of threads to run 
    [Default] 1
    
    -w or --workdir
    The working directory to put temporary files and results
    [Default] $PWD/hypo_eval_wd
```

### Assembly contiguity evaluation for diploid genome

The requirements for diploid genome error evaluation is a diploid reference genome, paternal and maternal.

```console
 Usage: python eval_script.py diploid_contiguity <args>

** Mandatory arguments
    -1 or --assembly1
    The file containing the first assembled haplotype
    
    -2 or --assembly2
    The file containing the second assembled haplotype
    
    -p or --paternal
    The file containing the paternal reference
    [Note] For the evaluation purposes, maternal and paternal are interchangable and only reflected on the report
    
    -m or --maternal
    The file containing the maternal reference
    [Note] For the evaluation purposes, maternal and paternal are interchangable and only reflected on the report

** Optional parameters
    -t or --threads
    The number of threads to run 
    [Default] 1
    
    -w or --workdir
    The working directory to put temporary files and results
    [Default] $PWD/hypo_eval_wd
    
    -r or --relaxed
    Set a relaxed condition on haplotype evaluation across contigs, i.e. ignore haplotype switches across contigs
```

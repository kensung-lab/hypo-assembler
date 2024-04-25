# Hypo-Assembler Evaluations

In this page we show the evaluations of hypo-assembler results compared to other assemblers, in order to reproduce the results as shown in our paper. The relevant scripts are included in this directory.

## Contiguity Evaluations

The evaluations are done using [QUAST](https://github.com/ablab/quast). 

The reference genomes used in the evaluations are [paternal](https://www.ncbi.nlm.nih.gov/assembly/GCA_021951015.1/) and [maternal](https://www.ncbi.nlm.nih.gov/assembly/GCA_021950905.1).

The script [separate_hap.py](https://github.com/kensung-lab/hypo-assembler/tree/main/evaluation_script/evaluate_hap.py) is first used to separate the haplotypes based on the reference genome, since we will be assessing the haplotype switching error separately. The better result (before or after the haplotype separation script) are shown.


| Assembler      | Total Size    | NGA50 (Mbp) | # Contigs | Coverage | # Misassemblies |
| -------------- | ------------- | ----------- | --------- | -------- | --------------- |
| Hypo-assembler | 2.932 / 2.918 | 50.9 / 50.9 | 259 / 259 | 97.8     |  12 / 12        |
| HifiAsm        | 3.033 / 3.015 | 54.5 / 40.3 | 445 / 402 | 96.9     |  41 / 45        |
| HifiAsm (Hi-C) | 3.075 / 2.909 | 44.9 / 45.6 | 405 / 399 | 97.1     | 46 / 50         |
| FALCON-Phase   | 3.027 / 3.027 | 29.5 / 29.3 | 503 / 512 | 95.6     | 115 / 129       |
| HifiAsm (Trio) | 2.936 / 3.033 | 55.7 / 55.6 | 351 / 307 | 97.9     | 31 / 35         | 
| Merfin Trio    | 2.985 / 2.913 | 48.2 / 47.4 | 319 / 288 | 97.7     | 78 / 81         |


## Error Evaluations

Error evaluations are done to the reference genome using the script [error_eval.py](https://github.com/kensung-lab/hypo-assembler/blob/main/evaluation_script/error_eval.py) The script measures errors that are not haplotype errors, measuring switching errors separately.

Assembly        | Haplotype | Mismatch | Insertion | Deletion | Switch  |
--------------- | --------- | -------- | --------- | -------- | ------  |
Hypo-assembler  | 1         | 0.60041  | 0.06504   | 0.01903  | 0.02943 |
Hypo-assembler  | 2         | 0.61082  | 0.05894   | 0.02889  | 0.02914 |
Hifiasm         | 1         | 0.66837  | 0.04915   | 0.03335  | 0.02543 |
Hifiasm         | 2         | 0.64583  | 0.04997   | 0.03755  | 0.02241 |
Hifiasm (Hi-C)  | 1         | 0.66837  | 0.04905   | 0.03345  | 0.02343 |
Hifiasm (Hi-C)  | 2         | 0.66533  | 0.04907   | 0.03655  | 0.02441 |
FALCON-Phase    | 1         | 3.79093  | 1.29039   | 0.69904  | 0.15772 |
FALCON-Phase    | 2         | 3.92378  | 1.23909   | 0.67125  | 0.16810 |
Hifiasm (Trio)  | 1         | 0.66082  | 0.05894   | 0.03889  | 0.02122 |
Hifiasm (Trio)  | 2         | 0.65763  | 0.06707   | 0.04614  | 0.02445 |
Merfin Trio     | 1         | 1.39405  | 0.64931   | 0.49502  | 0.09748 |
Merfin Trio     | 2         | 1.49850  | 0.69039   | 0.41294  | 0.09755 |

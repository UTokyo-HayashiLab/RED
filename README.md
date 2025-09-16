## RED: A Recombination-Aware Extension of Edit Distance

This repository contains the implementation of **RED**, a dynamic programming–based algorithm for detecting recombination breakpoints and parental origins in viral genomes.  
We also provide benchmark datasets, evaluation scripts, and comparison baselines (e.g., Bootscan, RED-Band).  

## Overview
Recombination is a key evolutionary process in many viruses, such as HIV-1 and hepatitis B virus (HBV). Accurate detection of recombination breakpoints is essential for phylogenetic analysis and molecular epidemiology.  

**RED** extends the classical edit distance to explicitly model *switches* between parental lineages, enabling precise inference of recombination breakpoints.  
Key advantages:
- No machine learning or training data required  
- Transparent and interpretable dynamic programming formulation  
- Higher breakpoint precision compared to Bootscan on both synthetic and real HBV benchmarks  

## Contents
- `red.c` — RED implementation 
- `red_band.c` — Band-limited RED for faster runtime  
- `bootscan.c` — Bootscan baseline  
- `eval.py` — Evaluation script for accuracy, precision, recall, F1, MAE  

## Usage
1. Compile Bootscan, RED and RED-Band:

   ```bash
   # Compile Bootscan
   gcc -O3 -o bootscan bootscan.c

   # Compile RED
   gcc -O3 -o red red.c

   # Compile RED-Band
   gcc -O3 -o red_band red_band.c

2. Running examples
   ```bash
   unzip 1200_synthetic_eval_varied.zip
   cd 1200_synthetic_eval_varied/sample_1
    
   # Run RED on sample_1
   ./../../red \
    --parents-aln parents.aligned.fasta \
    --child       child.fasta \
    --out         sample_1_red \
    --mismatch    1 \
    --gap         1 \
    --switch      20 \
    --win         21 \
    --kappa       1 \
    --stick       1 \
    --mu-inf      2 \
    --mu-eq       1 \
    --gap-swell   1 \
    --gap-near    10 \
    --reward      0 \

   # Run Bootscan on sample_1
   ./../../bootscan \
    --parents-aln parents.aligned.fasta \
    --child       child.fasta \
    --out         sample_1_bs \
    --window  21 \
    --step        1

   # Run RED-Band on sample_1
   ./../../red_band \
    --parents-aln parents.aligned.fasta \
    --child       child.fasta \
    --out         sample_1_red_band \
    --mismatch    1 \
    --gap          1 \
    --switch      20 \
    --win         21 \
    --kappa       1 \
    --stick       1 \
    --mu-inf      2 \
    --mu-eq       1 \
    --gap-swell   1 \
    --gap-near    10 \
    --reward      0

3. Evaluation
   ```bash
   python eval.py \
   --red-origins 1200_synthetic_eval_varied/sample_1/sample_1_red.red.origins.tsv \
   --truth-breaks 1200_synthetic_eval_varied/sample_1/truth.breakpoints.json \
   --truth-origin 1200_synthetic_eval_varied/sample_1/truth.origin.tsv 

## Citation

Katsuhiko Hayashi, RED: a recombination-aware edit distance for parental origin tracing in mosaic genomes, [in preparation].


## License

This project is licensed under the Creative Commons CC0 1.0 Universal (Public Domain Dedication).
CC0 1.0 Universal Summary
You can copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission.

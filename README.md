# ATLANTIDES_pipeline

This document gathers all my BASH and R scripts that I created and used during my master 2 project: Investigate the genetic effect of parallel domestication of Atlantic Salmon in Europe and North America: a genome-wide scan for selection signatures

Comments are written for problems I encountered or interesting things to know about the bioinformatics analyses.

Two mains objectives : 
- define the genetic structure of each population
- highligh genes under selection in farmed population by 3 methods : Tajima's D (site-frequency based method), FST (population differenciation), XP-EHH (haplotype differenciation based method)

Following steps 

1. Variants/genotypes calling 
2. SNPs extraction (only bi-allelic variants)
3. Hard filtering ("PASS" filtering)
4. Data quality statistics
5. Tajima calculation
6. Fst calculation
7. Soft filtering (according to step 4.)
8. Structure analyses
9. Phasing
10. Haplotypes analysis


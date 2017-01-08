# Evolution of viruses
Bioinformatics seminar

Group ub201617_beta



#### Calculating edit distance matrix

First we identified genes and their sizes. We found that Ebola virus consists of seven genes, four of which are smaller (named VP24, VP30, VP35 and VP40), two medium sized (named GP and NP) and one large gene (named L). Smaller ones are of length ~1500 bp, medium sized ~2500-3000 bp and larger one ~7000 bp. 

Next we randomly pick two samples and perform global alignment on each gene using dynamic programming. Since we are only interested in edit distance, we do not perform traceback. We found that calculation time is in the order of couple of minutes (for a single pair of genes) and memory requirements are significant. It was not possible to calculate edit distance of the largest gene L on a machine with 8GB of memory.

We ignore gene L due to memory limitations and we use multiprocessing to speed up calculations. Calculation of whole edit distance matrix took ~20 hours.
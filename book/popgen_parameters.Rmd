# Population genetics indices

## Computing population genetics parameters in STRAF

<center><img src="img/capture_popgen_parameters_1.png" class="capture"/></center>

## Details on population genetics indices

### Hardy-Weinberg equilibrium

A population is considered at Hardy-Weinberg equilibrium when

Why is it important to check? If a locus presents a significant deviation from
HWE, it means that a process is influencing the distribution of allele and genotype
frequencies in the population.

__Inbreeding__


__Population structure__

Individuals closer from each other are in general more likely to mate with each other.

__Selection__

Locus is under selection. This is very unlikely that STR loci as they are 
supposed to evolve neutrally. However, they could be found near loci under selection

In forensics, we need to assume HWE as indices computed (for example, a match probability)
would be biased if

### Heterozygosities

### F-statistics

$F_{\textrm{IS}}$ 

$F_{\textrm{ST}}$


:::note
__One concept, multiple estimators.__

Several __estimators__ of $F_{\textrm{ST}}$ exist (for example, Weir and Cockerham's, Nei's, 
Hudson's $F_{\textrm{ST}}$). It's like if each population geneticist decided to develop their
own estimator! Why is that? In statistics, what we call an __estimator__
is. It is important to keep in mind that these estimators rely on a specific __model__,
with underlying assumptions. It explains why some estimators are more or less reliable
depending on the case and observed data, and each of them has been developed for 
a different situation.
:::

## Linkage disequilibrium (LD)

### What is linkage disequilibrium?

__Coming soon.__

### How to compute LD in STRAF?

__Coming soon.__

<center><img src="img/capture_popgen_parameters_2.png" class="capture"/></center>

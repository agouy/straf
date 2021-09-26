# Introduction {-}

In this chapter, we will briefly introduce some essential concepts in genetics.

## DNA and genetic variation {-}

Each of our cells contains 23 pairs of __chromosomes__, composed of a long __DNA__ 
(deoxyribonucleic acid) molecule. Under this somewhat barbaric name
is hiding a simple concept. This molecule is the support of the information 
used by the body to function and development. This information is encoded by a 
chain of __nucleotides__ of four types that can be referred to using the letters
A (Adenine), T (Thymine), C (Cytosine) and G (Guanine).

Developments in biotechnologies enabled the characterisation of the DNA of
individuals. These techniques also led to the discovery that this DNA varies
between individuals. This __genetic variation__, also called __polymorphism__, 
can be used to characterise individuals and populations based on their DNA.

## Markers of polymorphism {-}

DNA variation can take different forms: it can for example be a 
__Single Nucleotide Polymorphism__ (**SNP**), when a mutation occurs and changes a
nucleotide at a given position in the genome. In that case, we would observe
different nucleotides in a population at a single position. 

There can also be __insertions__ and __deletions__ (sometimes referred to as
**InDels**), of one or multiple nucleotides.

Finally, other markers of genetic variation are __Copy Number Variants__ (**CNVs**), 
when a sequence is repeated a certain number of times. 
They can contain more or less repetitive units. These units can contain more or 
less nucleotides. __Short Tandem Repeats__ (**STRs**) a type of genetic 
polymorphism consisting in short sequences from 2 to 7 base pairs that are 
repeated. The __number of repeats varies__ among individuals, therefore characterizing 
their length can be useful to identify individuals. STRs are still nowadays one of 
the most commons markers used in forensic genetics.

## Polymorphism and forensics {-}

### DNA profiling and typing {-}

As DNA varies between individuals, DNA typing became a central element of the forensic
scientist toolkit. For example, typical questions forensic genetics aims at
answering include:

* What is the probability that a randomly-picked person in a population would 
match the individual of interest in terms of DNA?

* Which proportion of the population has the same combination of
genetic variants as the sample of interest?

### The role of population genetics {-}

To answer these questions, it is crucial to first get a good characterisation 
of __genetic variants__ (or **alleles**) frequencies in populations of interest,
at different __loci__ across the genome. Indeed, these frequencies can vary 
widely among populations.

In this context, STRAF has been designed to facilitate the analysis of 
__population data__ in forensic genetics.


:::digression
__Digression - Why STRs and not other markers?__

Nowadays, it is easier and cheaper to generate whole-genome sequences. One could
wonder why STRs are still so popular in forensic practice, and have not been replaced
by SNPs that are easily generated from Next-Generation Sequencing (NGS) data.

This is mainly due to the fact that STRs have a high mutation rate, therefore
are more diverse in human populations. This explains their high __power of discrimination__.

Furthermore, they can also be used in deciphering mixture components, a very common
case in forensics.

Finally, they can be combined in multiplex assays, which is convenient when low 
amounts of biological material can be recovered.

For all these reasons, STRs remain the dominant marker used in forensic genetics.
:::

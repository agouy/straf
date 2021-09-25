--- 
title: "The STRAF Book"
author: "Alexandre Gouy and Martin Zieger"
date: "2021-09-25"
site: bookdown::bookdown_site
documentclass: book
output: bookdown::bs4_book
classoption: a5paper
biblio-style: apalike
link-citations: yes
description: "Online version of the STRAF Book."
---

# Preface {-}

<center><img src="img/cover.png" class="cover" alt="book_cover" width="300"/></center>

## What is this book? {-}

This is the online version of __The STRAF Book__, which is currently under
active development. It is dedicated to the STRAF software, a web application
for the analysis of genetic data in forensics practice.

## Forensic and population genetics, lost sisters {-}

Genetics has many faces, and forensic and population genetics are two of them.
If we were to summarise their respective scopes, we could say that the former 
is the application of genetics to legal matters, and the latter aims at 
understanding genetic differences within and between populations, a fundamental matter in 
evolutionary biology.

Forensic genetics and population genetics have always been tightly linked
disciplines. This is likely because quite a number of questions they address
are similar. Even though problems in forensics and population genetics seem 
different, they often are the same question, simply phrased differently.

As an example, DNA profiling, used in criminal investigations or parental testing,
aims at matching different DNA samples and understanding how related are some
samples in terms of DNA. In population genetics, a common goal is to
characterise the genetic diversity of a set of populations, by looking at
how related individuals are within and between populations. Hence you can now
imagine why the two fields are linked: they both want to __understand and quantify__ 
the __relatedness__ of a set of samples.

Software and metrics developed in the population genetics for the study of the 
evolution of species are now used routinely in forensic genetics practice. 
But forensics is not just _applied population genetics_. The legal implications
and unique situations encountered in the forensics world also led to the 
development of relevant statistical tools and metrics with a more specific purpose.

## And then there was STRAF {-}

STRAF was born from the encounter of two scientists: a forensic geneticist and
a population geneticist. In 2017, in Bern, Switzerland, Martin Zieger came to 
visit a population genetics lab, where Alexandre Gouy was pursuing his 
Ph.D. thesis at that time.

This encounter led to a fruitful collaboration when they realised that some tools
used in population genetics could be leveraged by the forensics community. The
most striking example is the computation of forensics parameters, that describe
for example how good are our loci at discriminating samples. These 
parameters were typically computed using a spreadsheet that had been created by 
one of the suppliers of assays used to genotype samples. It is the mythical 
PowerStats v1.2 spreadsheet, allowing to compute forensic statistics and allele 
frequencies in Microsoft Excel. It has been since then removed from the Internet, 
and forensic geneticists started sharing this spreadsheet among each other, circulating 
almost secretly, "under the cloak" as French speakers would say.

As similar operations were done in routine in population genetics, we already had 
some scripts for the analysis of STR data. Then, after we applied them to an existing 
dataset, we decided to put everything into a web application so that the forensics 
community could benefit from it.

A few weeks later, STRAF was born, and after four year, STRAF had become a 
widely used tool by the forensics community, but not only. 
It has been used as a support for teaching population genetics, and has 
been used in evolutionary biology studies. 
The positive reception of the software in the community motivated its 
development over the years until the release of STRAF 2.0 in 2021. 

STRAF's story highlights the importance of communication between fields.

## What will you learn? {-}

By reading this book, our hope is that you will:

* Get an overview of common __concepts__ in forensic and population genetics

* Learn how to use the __STRAF software__ for STR data analysis through __practical applications__

* Be able to __interpret__ common metrics and analyses used in forensics practice

## Outline {-}

The book is organised as follow:

* We'll start by an __Introduction__ to essential forensic and population genetics concepts.

* In __Chapter 1__, we will focus on data, from its generation to its preparation for
downstream analysis in STRAF.

* In __Chapter 2__, we will review __forensic parameters__ that can be computed in STRAF,
and discuss their interpretation.

* In __Chapter 3__, we will review essential population genetics concepts and
describe __population genetics indices__ that can be computed in STRAF.

* In __Chapter 4__, we will focus on __multivariate statistics__ and how they can provide 
insights into population structure, with a particular focus on Principal Component
Analysis (PCA) and Multidimensional Scaling (MDS), two widely used approaches in genetics.

* In __Chapter 5__, we gather recommendations around potential next analysis steps
by presenting STRAF's __file conversion__ capabilities and useful methods implemented in
__other software__.

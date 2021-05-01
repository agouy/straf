STRAF is a browser-based application that allows to perform forensics and population genetics analysis of STR data.


### The STRAF Book

Click on the image below to access our online book with a lot more details about the software!

<div align='center'><a href='https://agouy.github.io/straf/' target='_blank'><img src='cover.png' align='center' width='269' /></a></div>

### Input file format

STRAF accepts tab delimited txt-tables containing genotypes.

* The first column, named __ind__, needs to contain the sample ID 
* The second column, , named __pop__, contains the population ID (this column must exist even if a single population is studied)
* The next columns correspond to genotypes: for haploid samples, one column per locus must be reported; for diploid data, two columns per locus (with the same name)
* Genotypes must be encoded as numbers (STRAF accepts point alleles)
* Missing data (e.g. null alleles) must be indicated with a “0”.

This format is designed to facilitate the input file generation from
a typical Excel file (Save as > Text (Tab-delimited) (*.txt)).

Examples of diploid and haploid input files can be downloaded using the 
following links:

<a id="raw-url" href="exampleSTRAFdiplo.txt" download="exampleSTRAFdiplo.txt" target="_blank">Download STRAF diploid example file.</a>

<a id="raw-url" href="exampleSTRAFhaplo.txt" download="exampleSTRAFhaplo.txt" target="_blank">Download STRAF haploid example file.</a>

### Using STRAF

STRAF computes standard forensics parameters. Some standard population
genetics analysis can be achieved if the samples are assigned to different 
populations. STRAF generates downloadable tables, and plots can be personnalized 
(using the Graphical parameters section on the left panel)
before saving them. Details about the methods can be found in Gouy & Zieger (2017).

Use the left panel to choose and upload your file. If no error appear
once the file is uploaded, three tabs appear on the right page: 
Data, Forensics analysis and Population genetics analysis.

On the Data tab, three checkboxes allow to 1. display the dataset, 
2. plot the distribution of alleles frequencies per locus, and 3. display a 
table of allele frequencies. This table is formatted as in most forensics data
reports (rows = alleles; columns = loci). You can download this table as a 
TSV file readable in Excel by clicking the Download button.

On the Forensics analysis tab, you can compute all the standard
   forensics parameters for your dataset. You can download the table
   and plot the results.
   
On the Population genetics analysis tab, standard population
   genetics statistics can be computed (F-statistics, HWE, LD, ...).
   A PCA can also be performed to study population structure or
   discover outliers.

### Reference

Please cite STRAF if you use it for your project
using the following reference:

> Gouy, A., & Zieger, M. (2017). STRAF - A convenient online tool for STR data evaluation in forensic genetics. Forensic Science International: Genetics, 30, 148-151.


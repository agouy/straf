### Udpates

* 2.0.7 (03/10/2021) – STRAF includes now more robust HW and LD tests. STRAF calls
Genepop to run the tests. Graphical representations are now interactive. Minor
bug fixes (FST matrix disappeared).
* 2.0.0 (01/05/2021) – STRAF has been partially reimplemented as an R package, implying a lot of code refactoring. MDS with reference population improved. The STRAF book has been populated. Genepop-like format is not supported anymore and input parameters have been simplified (only the "ploidy"" input parameter remains). Overall performance has been optimised.
* 1.4.6 (25/04/2021) – It is now possible to add custom reference frequency data to perform an MDS.
* 1.4.5 (25/04/2021) – Genepop conversion now works for haploid data. Common alleles are displayed for MDS.
* 1.4.4 (24/04/2021) – Various UI improvements (select and checkboxes inputs), code optimisation. MDS now excluding alleles with at least 1 missing value instead of imputing with the mean.
* 1.4.3 (18/04/2021) – It is now possible to add uploaded population to the STRidER MDS. Minor UI improvements.
* 1.4.2 (17/04/2021) – MDS on STRidER reference database has been added.
* 1.4.1 (08/03/2021) – STRAF now accepts point alleles for haploid data. In that case, all allele values are multiplied by 10.
* 1.4.0 (09/02/2021) – STRAF can now perform an MDS. All results can be downloaded as Excel files. Minor bug fixes in file conversion.
* 1.3.3 (03/02/2021) – Minor bug fixes in file conversion and PIC computation. Improved graphics.
* 1.3.2 (02/02/2021) – STRAF can convert files to the Arlequin format.
* 1.3.0 (01/02/2021) – STRAF can convert files to the Genepop and Familias formats. A File conversion tab has been added.
* 1.2.2 (09/01/2021) – STRAF has moved to an AWS server (without any changes in the License).
* 1.1.2 (03/10/2020) – Preparation for integration to Tercen; few minor updates
* 1.0.5 (18/09/2018) – a few bug fixes; new page with updates, license, data usage
* 1.0.4 (30/01/2018) - percentages of explained variance appear now on the PCA plot + style improvements
* 1.0.3 (18/12/2017) - a few bugs fixed, and it is now possible to analyse a single locus
* 1.0.2 (20/10/2017) - forensics parameters and population genetics indices can be computed for each population separately; PCA coordinates and eigenvectors can be downloaded
* 1.0.1 (22/09/2017) - allele frequencies can now be computed and downloaded for each population separately
* 1.0.0 (02/06/2017) - STRAF 1.0.0 release

### Future improvements

* Genepop call for LD and HW exact tests
* HGDP database integration to references
* HGPD matching probabilities
* Formal testing procedures
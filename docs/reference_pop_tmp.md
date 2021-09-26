# Reference populations analysis

__Coming soon.__

## Approach

MDS on allele frequencies.

By default, STRIDER frequencies (loci with less than 10 populations have been excluded).

<center><img src="img/capture_ref_pop_1.png" class="capture"></center>

<center><img src="img/capture_ref_pop_2.png" class="capture"></center>

Then it is possible to add samples to the existing MDS.

<center><img src="img/capture_ref_pop_3.png" class="capture"></center>



## Preparing a custom allele frequency database

Custom database:

| D1S1656 |            |         |        |
|--------|-------------|---------|--------|
| Allele | Switzerland | Germany | France |
| 9      | 0.12        | 0.09    | 0      |
| 10     | 0.40        | 0.35    | 0.28   |
| 10.2   | 0.31        | 0.41    | 0.5    |
| 11     | 0.17        | 0.15    | 0.22   |
| __D2S1338__ |        |         |        |
| Allele | Switzerland | Germany | France |
| 19     | 0.40        | 0.38    | 0.42   |
| 20     | 0.42        | 0.26    | 0.28   |
| 21     | 0.31        | 0.36    | 0.3    |


It only takes a few steps to generate an input file in a format that is suitable
for use in STRAF. From Excel, for example, we can start from a spreadsheet looking 
like this:

<center><img src="img/capture_ref_excel_1.png" class="capture_small"></center>

Then, one simply needs to save this table as a tab-delimited text file. This can be
achieved by clicking on `Save As` > `Text (Tab-delimited) (*.txt)`

<center><img src="img/capture_ref_excel_2.png" class="capture"/></center>


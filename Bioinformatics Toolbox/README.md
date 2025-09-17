# Bioinformatics Toolbox
## Code:
- [working_wth_fasta.ipynb](https://github.com/OVER-CODER/OE3E40-CGP/blob/1901eea76b881c83e82005ce9b3f85a75ebb73b2/Bioinformatics%20Toolbox/working_wth_fasta.ipynb)

## Get Databases Functions:
- (getgenpept): GenPept
- (getgenbank): GenBank
- (getembl): European Molecular Bio lab
- (getpir): Protein Sequence Database PIR-PSD
- (getpdb): Protein Data Bank PDB
- (getgeodata): NCBI Gene Expression Omnibus(GEO)


## Important functions:
Mitochondria genbank sequence
```
mitochondria = getgenbank('NC_012920.1','SequenceOnly',true)
```

To determine mono-, di-, and trinucleotide content, and to locate open reading frames.
```
aacount, basecount, codoncount, dimercount, nmercount, ntdensity 
```
Search for specific patterns within a sequence
```
seqshowwords, seqwordcount
```

Search for open reading frames 
```
seqshoworfs
```

Plot monomer densities and combined monomer densities in a graph. In the MATLAB Command window, type 
```
ntdensity(mitochondria) 
```
Count the nucleotides using the function basecount.           
```
basecount(mitochondria)
```
Count the nucleotides in the reverse complement of a sequence using the function seqrcomplement.
```
 basecount(seqrcomplement(mitochondria))
```

Use the function basecount with the chart option to visualize the nucleotide distribution.
```
basecount(mitochondria,'chart','pie’); 
```
Count the dimers in a sequence and display the information in a bar chart.
```
dimercount(mitochondria,'chart','bar’)
```
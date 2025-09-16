# Bioinformatics Toolbox
## Get Databases Functions:
- (getgenpept): GenPept
- (getgenbank): GenBank
- (getembl): European Molecular Bio lab
- (getpir): Protein Sequence Database PIR-PSD
- (getpdb): Protein Data Bank PDB
- (getgeodata): NCBI Gene Expression Omnibus(GEO)


## Important functions:
- to determine mono-, di-, and trinucleotide content, and to locate open reading frames.
```
aacount, basecount, codoncount, dimercount, nmercount, ntdensity 
```
- search for specific patterns within a sequence
```
seqshowwords, seqwordcount
```

- search for open reading frames 
```
seqshoworfs
```

- Plot monomer densities and combined monomer densities in a graph. In the MATLAB Command window, type 
```
ntdensity(mitochondria) 
```

-  Count the nucleotides using the function basecount.           
```
basecount(mitochondria)
```
- Count the nucleotides in the reverse complement of a sequence using the function seqrcomplement.
```
 basecount(seqrcomplement(mitochondria))
```
# virlink  
A package for analyze of VirScan data for cross-reactivity and co-occurrence  

This packages contains functions to calculate pairwise sequence similarity scores and z-score correlations for VirScan profile.  

### To install this package:  
```
{
if(!require(devtools)) install.packages("devtools")
devtools::install_github(repo = "siyangxia419/virlink", 
                         ref = "main",
                         upgrade = "never")
}
```

### The current version of the pakcage (0.2.1) contains 6 functions:  
`get_kmer`: generate kmers of user-selected length from a longer peptide sequence  
`make_kmer_dataframe`: generate a kmer dataset from a peptide dataset using `get_kmer`  
`peptide_cooccurrence`: calculuate a contigency table of two peptide's presence in multiple individuals  
`to_pairwise`: given a data frame, perform pairwise analysis with any function on all pairs of columns  
`peptide_pairwise_alignment`: perform pairwise sequence alignment for a data frame of peptides, using `to_pairwise` and `pairwiseAlignment` in the package "Biostrings"  
`peptide_pairwise_correlation`: perform pairwise correlation or co-occurrence analysis for a data frame of peptides, using `to_pairwise`.  

The package also contains two example datasets: `peptide_df` and `peptide_z`  

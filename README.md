# Kinase-specificity-evolution
A repository of files relating to the *Bradley and Beltrao, 2018* bioRxiv manuscript about the evolution of protein kinase substrate recognition.

## Protein kinase sequence alignment

A FASTA file of protein kinase sequences is given in the following file:
```
kinase_domain_opisthokont.fa
```
This file is essentially a compilation of all kinase domain sequences found in *KinBase* for its 9 animal and fungal species. The MSA of these kinase domain sequences is represented in the following file:
```
opistho_linsi_alcorrect_trim80_filter190_mod.fa
```
The name of each sequence has been modified to include the Group/Family/Subfamily classifications provided by *KinBase*. Pseudokinases were removed from the original sequence file, and then an alignment was generated using the MAFFT-LINSI method. Manual adjustments to the alignment were made in *JalView*. Alignment positions with more than 20% 'gap' characters were filtered from the MSA using *trimAL*. Finally, truncated kinase domain sequences with fewer than 190 amino acids (~75% of the kinase domain) were also filtered. This MSA contains 2094 rows (kinases) and 246 columns (positions). 

## Protein kinase global phylogeny

The trimmed MSA referred to above was used to generate a maximum likelihood (ML) phylogeny of the 9 opisthokont kinomes that were used for this analysis. The ML phylogeny was generated using the RAxML tool; amino acids substitutions were modelled using the LG substitution matrix and a gamma model to account for the heterogeneity of rates between sites. A neighbour joining (NJ) generated with the *ape* package in *R* was used as the starting tree. The resulting phylogeny is found in the following file:
```
RAxML_bestTree.alcorrect_linsi190_NJ
```
We recommend *FigTree* and *iTOL* as particularly good tools for the visualisation of phylogenies.

## Protein kinase ancestral sequence reconstruction

Ancestral sequences were then reconstructed for every node in the global phylogeny. This was achieved using the *CodeML* program of the *PAML* package. As with the phylogenetic reconstruction, an LG substitution matrix and a gamma model were employed. No molecular clock was assumed, and four discrete categories were specificied for the gamma model. The *CodeML* configuration file used for this analysis has the following name:
```
bradley_opistho_alcorrect_NJstarter.tcl
```
As an output, *CodeML* produces a phylogeny with the different ancestral nodes labelled. The name of this file is as follows:
```
opisthokont_alcorrect_midpoint_NJstart.tre
```
It is important to note that this phylogeny does **not** contain the correct branch lengths that were produced by RAxML.

Posterior probabilities for every alignment site and every ancestral node were outputted by *CodeML* in a file called 'rst'. Please follow the link to download the file (**warning: file size is 1.2 GB**):
```
https://drive.google.com/open?id=1CrD1AWOBcxOIoUy_p1bH730Z7Hd0aJD4
```
We provide as an '.csv' file the predicted amino acids for every ancestral *Family* and *Subfamily* node tested . This includes the posterior probabilities for the most likely amino acid at each site:
```
Family_posterior_probabilities_corrected.csv
Subfamily_posterior_probabilities_corrected.csv
```
The predicted sequences for every ancestral node in the phylogeny is also given in the following FASTA file:
```
ancestral_seq.fa
```
We also include here a very simple custom *R* function called *anc_node query()* that can be used to parse the 'rst' file and extract posterior probabilities for any site and node of interest:
```
anc_node_query.R
```
## Calculation of divergence scores and mapping to protein kinase structures

The divergence scores calculated for every *Family* and *Subfamily* tested is present in the following '.csv' files:
```
tested
```
The results of this phylogenetic analysis were then mapped to the 3D structure of the protein kinase (as presented in Figure 2). The correspondence of the kinase MSA positions to kinase domain positions, and the kinase domain positions to the protein kinase A structural numbering, is given in the following file:
```
Bradley_kinase_mapping.csv
```
To make the underlying method wmore accessible, all steps stated above (sequence alignment -> tree generation -> ancestral sequence reconstruction -> divergence score calculation -> 3D mapping) have been repeated in the same way for a simplified example (PLK vs. AUR kinases) in the following markdown file:
```
test
```
The R code used for the analysis of this 'toy' example is provided in the following file:
```
test
```


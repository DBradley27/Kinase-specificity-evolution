# Kinase-specificity-evolution
Files relating to the *Bradley and Beltrao, 2018* bioRxiv manuscript (about the evolution of protein kinase substrate recognition) have been deposited here.

## Protein kinase sequence alignment

A list of protein kinase sequences is given in the following FASTA file:
```
kinase_domain_opisthokont.fa
```
This file is a compilation of all kinase domain sequences found in *KinBase* for its 9 animal and fungal species. The MSA of these kinase domain sequences is given in the following file:
```
opistho_linsi_alcorrect_trim80_filter190_mod.fa
```
The name of each sequence has been modified to include the Group/Family/Subfamily classifications provided by *KinBase*. Pseudokinases were removed from the original sequence file, and then an alignment was generated using the MAFFT-LINSI method. Manual adjustments to the alignment were made in *JalView*. Alignment positions with more than 20% 'gap' characters were filtered from the MSA using *trimAL*. Finally, truncated kinase domain sequences with fewer than 190 amino acids (~75% of the kinase domain) were also filtered. This MSA contains 2094 rows (kinases) and 246 columns (positions). 

## Protein kinase global phylogeny

The trimmed MSA referred to above was used to generate a maximum likelihood (ML) phylogeny of the 9 opisthokont kinomes that were used for this analysis. The ML phylogeny was generated using the RAxML tool; amino acids substitutions were modelled using the LG substitution matrix and a gamma model to account for the heterogeneity of rates between sites. A neighbour joining (NJ) phylogeny generated with the *ape* package in *R* was used as the starting tree. The resulting ML phylogeny is present in the following file:
```
RAxML_bestTree.alcorrect_linsi190_NJ
```
We recommend *FigTree* and *iTOL* as particularly good tools for the visualisation of phylogenies.

## Protein kinase ancestral sequence reconstruction

Ancestral sequences were then reconstructed for every node in the global phylogeny. This was achieved using the *CodeML* program of the *PAML* package. As with the phylogenetic reconstruction, an LG substitution matrix and a gamma model were employed. No molecular clock was assumed, and four discrete categories were specifiied for the gamma model. The *CodeML* configuration file used for this analysis has the following name:
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
We also include here a very simple custom *R* function called *anc_node query()* that can be used to parse the 'rst' file and extract posterior probabilities for any node-site combination of interest:
```
anc_node_query.R
```
## Calculation of divergence scores and mapping to protein kinase structures

The divergence scores calculated for every *Family* and *Subfamily* tested are present in the following '.csv' files:
```
Family_scores.csv
Subfamily_scores.csv
```
The results of this phylogenetic analysis were then mapped to the 3D structure of the protein kinase (as presented in Figure 2). The correspondence of the kinase MSA positions to kinase domain positions, and the kinase domain positions to the protein kinase A structural numbering, is given in the following file:
```
Bradley_kinase_mapping.csv
```
'-' positions in the 'PRKACA_Hs_no.' column represent deletions in the human PRKACA sequence relative to the rest of the MSA. Conversely, '-' positions in the 'pkafm_kin_no' column represent insertions in the human PRKACA sequence relative to the protein kinase domain. Discontinuities in 'PRKACA_Hs_no.' and/or 'pkafm_kin_no' numbering represent PRKACA positions and/or kinase domain positions that were removed from the trimmed alignment. 

To make the underlying method more accessible, all steps stated above (sequence alignment -> tree generation -> ancestral sequence reconstruction -> divergence score calculation -> 3D mapping) have been repeated in the same way for a simplified example (PKA vs. PKG kinases) in the following file:
```
Kinase_analysis_example.pdf
```
The R code used for the analysis of this PKA/PKG example is provided in the following markdown file:
```
Kinase_analysis_example.Rmd
```

## Evolutionry analysis of phosphorylation data

All phosphorylation sites used for the evolutionary analysis of phosphorylation motifs are present in the following files:

```
fg_list_publication.rds
fg_list_publication_prokaryotic.rds
```
Phosphorylation motifs were identified in each species using an R-based implementation of the *motif-x* method. For all motifs passing certain criteria (specified in the manuscript), binomial p-values were used to calculate the extent of the phosphomotif enrichment relative to a background set of randomly shuffled phosphorylation peptides. To clearly illustrate the methods used, this analysis has been repeated for a reduced dataset containing only four species (*Plasmodium falciparum*, *Plasmodium berghei*, *Toxoplasma gondii*, and *Tetrahymena thermophila*). All R code used for this analysis is present in the following two documents:
```
Phosphomotifs_analysis_example.Rmd
Phosphomotifs_analysis_example.pdf
```



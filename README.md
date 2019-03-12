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

# anc_seq_data: 'rst' file produced by CodeML
# node: the ancestral node of interest. A phylogeny with ancestral nodes labelled should have been produced by CodeML
# site: one site of interest from the MSA that was used to generate the ancestral sequence reconstructions.


anc_node_query <- function(anc_seq_data,node_number,site) {

  x <- readLines(anc_seq_data)
  
  # These are the character vectors that we will use to retrieve the posterior probabilities of interst
  
  grepcom <- paste('Prob distribution at node ',node_number,', by site',sep='')
  grepcom2 <- paste('Prob distribution at node ',node_number+1,', by site',sep='')
  
  # Use grep functions to find the lines of interest from the 'rst' file
  
  query_anc_grep <- grep(grepcom,x)
  query_anc_grep_end <- grep(grepcom2,x)
  
  # Retrieve the sequences of interest
  anc_seq_node <- x[query_anc_grep:query_anc_grep_end]
  
  # For the node of interest, find the posterior probabilities for each site and amino acid
  PPs <- rapply(strsplit(anc_seq_node[grep(': A',anc_seq_node)],split=':'),function(x) x[2])
  
  # Convert the posterior probabilities into a numeric vector
  PPs_numeric <- regmatches(PPs,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",PPs))
  PPs_numeric <- lapply(PPs_numeric, function(x) as.numeric(x))
  
  aa_probs <- unlist(strsplit(PPs[site],split=' '))[-1]
  
  # Generate an ordered vector of posterior probabilities for a given node and site
  ordered_ancestrals <- aa_probs[order(PPs_numeric[[site]],decreasing=TRUE)]
  
  return(ordered_ancestrals)
  
}

# As an example, we will find the posterior probabilities of the 20 amino acids at site 124 for ancestral node 3966, which lies
# at the base of the 'GRK' subfamily. Site 124 corresponds to the third residue of the KP_N motif.

example <- anc_node_query('rst',3966,124)

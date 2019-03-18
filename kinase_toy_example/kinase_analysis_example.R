# Set the working directory
setwd("~/Documents/Work/Post-PhD/Manuscript_2_revisions/Kinase_toy_example")

# Load in the required packages
detach("package:seqinr", unload=TRUE)
library('bio3d')
library(ape)

# Read in the 'rst' file containg the ancestral posterior probabilities
x <- readLines('rst_example')

# Read in the kinase alignment file

fas <- read.fasta('PKA_PKG_AKT_named_PRKACA_al_trim_filter.fa')

# Read in the tree with ancestral nodes annotated by codeml

tree <- read.tree('kinase_example_codeml_annotated.tre')

tips <- tree$tip.label
nodes <- tree$node.label

# Specify the family of interest

family <- 'PKA'

# Retrieve the tip labels containing the family of interest
Group <- grep(family,tips)

# Families containing fewer than 5 sequences in this dataset are not considered further
if (length(Group) < 5) {next}

# Retrieve the sequences containing the family of interest
seqs <- (fas[[2]][grep(family,fas$id),])

# Separate into chunks

chunkvec <- 0
count <- 0

# Here I introduce a threshold that allows me to determine whether non-continuities in tip labels are due to the presence of 'contaminant'
# sequences within a clade or actually due to division of the family sequences into separate clades

thresh <- 10

# Here, the 'count' variable is used to track the clade membership of different sequences by recording discontinuities in the tip labels

for (i in 2:length(Group)) {
  
  if (Group[i]-Group[i-1] < thresh) {
    chunkvec <- c(chunkvec,count)
  }
  else {
    count = count+1
    chunkvec <- c(chunkvec,count)
  }
}

# Now we select from our group those sequences belonging to the largest clade

tab <- sort(table(chunkvec))
major_clade <- Group[which(chunkvec == names(tab)[length(tab)])]

# Clades containing fewer than 5 sequences in this dataset are not considered further
if(length(major_clade) < 5) {next}


# Determine the clade purity by calculating the 'frac' variable

globalpath <- nodepath(tree)

Query_anc <- getMRCA(tree,major_clade)
path <- nodepath(tree,from=getMRCA(tree,tips),to=Query_anc)
Query_anc_anc <- path[length(path)-1]

ancestral_path <- globalpath[grep(Query_anc,globalpath)]
frac <- length(intersect(major_clade,unlist(ancestral_path)))/length(ancestral_path)
major_clade_rev <- rev(major_clade)
frac_rev <- length(intersect(major_clade_rev,unlist(ancestral_path)))/length(ancestral_path)
count=1


# Now we need to find the ancestor node of the major clade

tips_new <- rapply(lapply(strsplit(tips[major_clade],split='_'), function(x) x[2:length(x)]), function(x) paste(x,collapse='_'))

if (rapply(strsplit(rownames(seqs), split=''), function(x) x[length(x)])[1] == '_') {
  tips_new <- paste(tips_new,'_',sep='')
} 


# Retrieve all sequences from the major clade

seqs <- seqs[rownames(seqs) %in% tips_new,]

# Find the ancestral node of the major clade

Query_anc <- getMRCA(tree,major_clade)
path <- nodepath(tree,from=getMRCA(tree,tips),to=Query_anc)

# Find the node directly anetecedent to the ancestral node

Query_anc_anc <- path[length(path)-1]
globalpath <- nodepath(tree)
ancestral_path <- globalpath[grep(Query_anc_anc,globalpath)]

# Find the sister clade and its ancestral node

treepos <- grep(Query_anc_anc,ancestral_path[[1]])
Sister_anc <- setdiff(unique(rapply(ancestral_path, function(x) x[treepos+count])),Query_anc)



# Now time to extract the posterior probabilites

aavec <- c('A', 'R', 'N', 'D', 'C', 'Q' , 'E' , 'G' , 'H' , 'I', 'L' , 'K' , 'M',  'F' , 'P',  'S',  'T',  'W',  'Y',  'V') 
library(Biostrings)

# For the query case

grepcom <- paste('Prob distribution at node ',Query_anc,', by site',sep='')
grepcom2 <- paste('Prob distribution at node ',Query_anc+1,', by site',sep='')

query_anc_grep <- grep(grepcom,x)
query_anc_grep_end <- grep(grepcom2,x)
anc_seq_node <- x[query_anc_grep:query_anc_grep_end]

PPs <- rapply(strsplit(anc_seq_node[grep(': A',anc_seq_node)],split=':'),function(x) x[2])
PPs <- regmatches(PPs,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",PPs))
PPs_query <- PPs

# For the sister clade

grepcom3 <- paste('Prob distribution at node ',Sister_anc,', by site',sep='')
grepcom4 <- paste('Prob distribution at node ',Sister_anc+1,', by site',sep='')

query_anc_grep <- grep(grepcom3,x)
query_anc_grep_end <- grep(grepcom4,x)
anc_seq_node <- x[query_anc_grep:query_anc_grep_end]

PPs <- rapply(strsplit(anc_seq_node[grep(': A',anc_seq_node)],split=':'),function(x) x[2])
PPs <- regmatches(PPs,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",PPs))
PPs_sister <- PPs

# Now we must iterate through every one of the kinase domain positions and calculate a divergence score between
# the query clade and the sister clade.

scorevec <- NULL

aavec <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

for (i in 1:length(PPs_query)) {
  
  # Most likely query aa (ancestral) and associated probability
  query_aa <- aavec[which(PPs_query[[i]] == max(PPs_query[[i]]))][1]
  query_aa_prob <- max(PPs_query[[i]])
  
  # Calculate the residue conservation at this site also
  query_aa_cons <- conserv(seqs,method='similarity')[i]
  query_product <- query_aa_cons
  
  # Most likely sister aa (ancestral) and associated probability
  sister_aa <- aavec[which(PPs_sister[[i]] == max(PPs_sister[[i]]))][1]
  sister_aa_prob <- NULL
  
  # The divergence score is calculated differently depending on whether the predicted ancestral amino acid
  # for the sister clade is the same or different from that of the query clade
  
  if (query_aa == sister_aa) {
    #Probability that ancestral amino acid in the sister is the same as that for the query
    multiplier <- 1
    sister_aa_prob <- max(PPs_sister[[i]])
  } else {
    #Probability that ancestral amino acid in the sister is different from that of the query
    multiplier <- -1
    sister_aa_prob <- sum(as.numeric(PPs_sister[[i]][which(aavec != query_aa)]))
  }
  
  sister_product <- multiplier*as.numeric(sister_aa_prob)
  
  score <- query_product - sister_product
  scorevec <- c(scorevec, score)
}


# Now we must finally map the results to structure

detach("package:bio3d", unload=TRUE)
library(seqinr)

noneg <- scorevec

fasta_al <- read.fasta('PKA_PKG_AKT_named_PRKACA_al.fa')

count = 0
almap <- numeric(3)

# Extract the human PRKACA sequence
PDB_al <-fasta_al[[length(fasta_al)]][1:length(fasta_al[[length(fasta_al)]])]

# Iterate through each alignment position, and record whenever there is a 'gap' in the human PRKACA sequence
for (n in 1:length(PDB_al)) {
  
  if (PDB_al[n] == "-") {
    
    count = count
    pair_vec <- c(n,'-','')
    almap <- rbind(almap,pair_vec)
  }
  
  else {
    
    count = count + 1
    pair_vec <- c(n,count,'')
    almap <- rbind(almap,pair_vec)
  }
}
almap <- almap[-1,]

# For mapping, we only need to retain positions that were kept in the trimmed alignment after trimAl processing

colmap <- read.table('colnumbering.txt',col.names=FALSE,stringsAsFactors = FALSE,sep=',')
colmap <- as.numeric(colmap[,1])
colmap <- colmap+1
trim_map <- almap[as.numeric(almap[,1]) %in% colmap,]

# Now the aggregated scores obtained previously can be mapped to the human PRKACA sequence
names(noneg) <- trim_map[,2]

# Extract the top 5 sites
noneg[order(scorevec,decreasing=TRUE)[1:5]]

# Extract the top 5 sites
noneg[order(scorevec,decreasing=TRUE)[1:5]]

# Subtract 1 to map to PDB 1atp

as.numeric(names(noneg[order(scorevec,decreasing=TRUE)[1:5]]))-1



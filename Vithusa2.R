source("http://bioconductor.org/biocLite.R")
biocLite("igraph")

library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)
library(seqinr)  # there is a conflict with a getSequence() in biomaRt package
library(biomaRt)
library(stringr)

# 1. Obtain proteins known to be implicated in cardio-vascular disease
# 2. Create a network of interacting proteins for each disease protein from STRING or similar
# 3. Download the snps associated with each protein in this network
# 4. Compute statistics such as number and type of snp for all proteins (not all snps will cause mutation)
# 5. Use igraph for network statistics and charts to display information

# ---------- last bit---------------
ppi <- file.path('C://R-files//UKCI2016','interactions1.tsv') %>% read.delim(na.strings='')  
ppi <- ppi[,c("X.node1","node2","combined_score")] 
g <- graph.data.frame(ppi)
snplist <-V(g)$name
ad <- get.adjacency(g)
nodesize=igraph::degree(g)# *5 # change the number if you want bigger nodes usually x 5

# normal genes are 'lightnlue' and sized 10; hubs are 'darkgray' color & sized 15"
# normal genes have degree of 1 - 29; hubs 30 - 50 
for (i in 1:length(nodesize))
  if(nodesize[i] <= 30){nodesize[i] <-10} else {nodesize[i]<-15}

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
for (i in 1:length(nodecolor))
  if(nodesize[i] <= 10){nodecolor[i] <- "lightblue"} else {nodecolor[i]<-"darkgray"}

nodelabel<-V(g)$name

g<-as.undirected(g); 
plot(g, edge.color="darkgray", 
     vertex.color=nodecolor,
     vertex.label=nodelabel,
     vertex.label.cex=0.6, 
     vertex.label.font=0.5, 
     vertex.frame.color="darkgreen",
     vertex.size=18,
     vertex.label.color="black", 
     vertex.label.family = "sans", 
     layout=layout.kamada.kawai(g))

gs<-get_gstatistics(g)
head(gs)

tkplot(g,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,
       edge.arrow.size=0, edge.curved=FALSE)

# --------------------------------------- data structures ----------------------------
dstrings1 <- c("CDKN2B","FMN2","MTHFD1L")

#---------------------------------------- Protein Networks --------------------------------
# Each disease gene(protein) was entered into STRING with no more than 50 interactions selected and
# was saved as a text file.

CDKN2B <- file.path('C://R-files//UKCI2016','CDKN2B.txt') %>% read.delim(na.strings='')  # load it in.
CDKN2B <- CDKN2B[,c("X.node1","node2","combined_score")] # just keep PPI and score
cdkn2b <-graph.data.frame(CDKN2B)       # convert to igraph object
cdkn2b.ppi <-V(cdkn2b)$name          # obtain all protein names to get a list snps for.

SMAD3 <- file.path('C://R-files//UKCI2016','SMAD3.txt') %>% read.delim(na.strings='')  # load it in.
SMAD3 <- CDKN2B[,c("X.node1","node2","combined_score")] # just keep PPI and score
smad3 <-graph.data.frame(SMAD3)       # convert to igraph object
smad3.ppi <-V(smad3)$name          # obtain all protein names to get a list snps for.

# --------------------------------------- Connect to BiomaRt ENSEMBLE ----------------------------------
# obtain details of human proteins database and save in "mart", then get lists of proteins(ensembl ID) and save in "ENSG"
mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ENSG <- getBM(mart=mart, attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"), 
              filters="hgnc_symbol", values=cdkn2b.ppi)

#--- in the above code change contents of dstrings1 with handcrafted list of proteins you want OR use list
#--- like cdkn2a.ppi instead. ENSG is then passed on to next stage.

#--------------------------------------- Connect to BiomaRt SNPS --------------------------------------
# open access to the remote ensembl database using biomaRt
snpmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

#------------------------------- Get the SNPS for our list of proteins --------------------------------
# Create a large dataframe with fields for: 
# hgnc_symbol; total_snp_num; entries for each of the 18 types of snp;
snps<-data.frame()  # create an empty dataframe
for (i in 1:length(ENSG[,1])) {
  ens_gene <- ENSG[i,1]
  snp1 <- collate_snps(ens_gene)   # test out new collate_snp() function
  if(is.null(snp1)){cat(paste("No SNPs found for ",ENSG[i,1],"\n"));
    snp1 <- cbind(snp1,ENSG[i,2]);
    snps<-rbind(snps,snp1)}else{
      cat(paste(ENSG[i,1]," has ", snp1$total[1]," SNPs\n"))
      snp1 <- cbind(snp1,ENSG[i,2])  # add the human gene name rather than use ensemble name
      snps<-rbind(snps,snp1)}   # add next snps details to dataframe
}

# http://www.r-bloggers.com/irefr-ppi-data-access-from-r/
# --------------------------- visual notification ------------------------------------------------
# system('CMD /C "ECHO Hi Ken, your R Script or whatever it is has finished running && PAUSE"',invisible=FALSE, wait=FALSE)


# collate_snps() assumes that the snpmart has already been setup and makes the appropriate calls to biomart.
# It's function is to gather the snp data on a single protein and return the details.
# usage: snp_struct <- collate_snps(ensemble_gene_id)
collate_snps <- function(ensggene){
  snptemp <- getBM(mart=snpmart, 
                   attributes=c("refsnp_id","associated_gene","allele","chr_name","chrom_start","chrom_end","consequence_type_tv"), 
                   filters="ensembl_gene", values=ensggene,checkFilters=FALSE)
  if(nrow(snptemp) == 0){null_snps<-no_snps();return(null_snps)}
  
  # Here we have some really inelegant programming but I dont know a quick way of doing! So brute force!
  # create a dataframe holding all the snp possibilities, set that zero
  all_snps <- data.frame(downstream_gene_variant=0,
                         intron_variant=0,
                         upstream_gene_variant=0,
                         non_coding_transcript_variant=0,
                         splice_region_variant=0,
                         missense_variant=0,
                         three_prime_UTR_variant=0, # R does not like variables starting with a number
                         five_prime_UTR_variant=0,    # ditto
                         coding_sequence_variant=0,
                         non_coding_transcript_exon_variant=0,
                         splice_donor_variant=0,
                         synonymous_variant=0,
                         NMD_transcript_variant=0,
                         frameshift_variant =0,
                         stop_gained=0,
                         inframe_deletion=0,
                         inframe_insertion=0,
                         splice_acceptor_variant=0,
                         total=sum(table(snptemp[,7])),
                         stringsAsFactors=FALSE)
  
  temp_table <- as.data.frame(table(snptemp[,7]))   # sntemp[7] contains the consequencetype
  
  for (i in 1:nrow(temp_table)) {
    if(temp_table$Var1[i]=='3_prime_UTR_variant')
    {all_snps$three_prime_UTR_variant=temp_table[i,2]}
    if(temp_table$Var1[i]=='5_prime_UTR_variant')
    {all_snps$five_prime_UTR_variant=temp_table[i,2]}
    
    if(temp_table$Var1[i]=='downstream_gene_variant')
    {all_snps$downstream_gene_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='intron_variant')
    {all_snps$intron_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='upstream_gene_variant')
    {all_snps$upstream_gene_variant=temp_table[i,2]}
    
    if(temp_table$Var1[i]=='non_coding_transcript_variant')
    {all_snps$non_coding_transcript_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='splice_region_variant')
    {all_snps$splice_region_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='missense_variant')
    {all_snps$missense_variant=temp_table[i,2]}    
    
    if(temp_table$Var1[i]=='coding_sequence_variant')
    {all_snps$coding_sequence_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='splice_donor_variant')
    {all_snps$splice_donor_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='synonymous_variant')
    {all_snps$synonymous_variant=temp_table[i,2]}        
    
    if(temp_table$Var1[i]=='NMD_transcript_variant')
    {all_snps$NMD_transcript_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='frameshift_variant')
    {all_snps$frameshift_variant=temp_table[i,2]}  
    if(temp_table$Var1[i]=='stop_gained')
    {all_snps$stop_gained=temp_table[i,2]}        
    
    if(temp_table$Var1[i]=='inframe_deletion')
    {all_snps$inframe_deletion=temp_table[i,2]}  
    if(temp_table$Var1[i]=='inframe_insertion')
    {all_snps$inframe_insertion=temp_table[i,2]}  
    if(temp_table$Var1[i]=='splice_acceptor_variant')
    {all_snps$splice_acceptor_variant=temp_table[i,2]}        
    
  }
  
  return(all_snps)
}

# create a dataframe with all zero entries.
no_snps <- function(){
  null_snp  <- data.frame(downstream_gene_variant=0,
                          intron_variant=0,
                          upstream_gene_variant=0,
                          non_coding_transcript_variant=0,
                          splice_region_variant=0,
                          missense_variant=0,
                          three_prime_UTR_variant=0, # R does not like variables starting with a number
                          five_prime_UTR_variant=0,    # ditto
                          coding_sequence_variant=0,
                          non_coding_transcript_exon_variant=0,
                          splice_donor_variant=0,
                          synonymous_variant=0,
                          NMD_transcript_variant=0,
                          frameshift_variant =0,
                          stop_gained=0,
                          inframe_deletion=0,
                          inframe_insertion=0,
                          splice_acceptor_variant=0,
                          total=0,
                          stringsAsFactors=FALSE)
  return(null_snp)  
}

library(dplyr)
library(tidyr)
library(igraph)

# obviously change the path.

CDKN2B <- file.path('~//Desktop//UKCI2016','CDKN2B.txt') %>% read.delim(na.strings='')  # load it in.
CDKN2B <- CDKN2B[,c("X.node1","node2","combined_score")] # just keep PPI and score
g <-graph.data.frame(CDKN2B)       # convert to igraph object

head(CDKN2B)

ad <- get.adjacency(smad3)
nodesize=igraph::degree(g)# *5 # change the number if you want bigger nodes usually x 5

# normal genes are 'lightnlue' and sized 10; hubs are 'darkgray' color & sized 15"
# normal genes have degree of 1 - 29; hubs 30 - 50 
for (i in 1:length(nodesize))
  if(nodesize[i] <= 30){nodesize[i] <-10} else {nodesize[i]<-15}

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
for (i in 1:length(nodecolor))
  if(nodesize[i] <= 10){nodecolor[i] <- "lightblue"} else {nodecolor[i]<-"darkgray"}

nodelabel<-V(g)$name

g<-as.undirected(g); 
plot(g, edge.color="darkgray", 
     vertex.color=nodecolor,
     vertex.label=nodelabel,
     vertex.label.cex=0.6, 
     vertex.label.font=0.5, 
     vertex.frame.color="darkgreen",
     vertex.size=18,
     vertex.label.color="black", 
     vertex.label.family = "sans", 
     layout=layout.kamada.kawai(g))

gs<-get_gstatistics(g)
head(gs)

# Calculate some statistics about the disease gene network
get_gstatistics <- function(gt) {
  gstats <- data.frame(
    modularity=modularity(gt, membership(cluster_walktrap(gt))),
    avepath=average.path.length(gt),
    nedges=ecount(gt),
    nverts=vcount(gt),
    transit=transitivity(gt),
    avedegree=mean(degree(gt)),
    diameter=diameter(gt,weights=NA),
    connect=is.connected(gt),
    closeness=closeness(gt),
    betweenness=betweenness(gt,directed=FALSE),
    density=graph.density(gt),
    hubness=hub_score(gt)$vector,
    authority=authority.score(gt)$vector)
  #power=bonpow(gt))
  return(gstats)
}











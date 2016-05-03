# vithusa.R
# http://kateto.net/network-visualization

library(dplyr)
library(tidyr)
library(igraph)
library(NCBI2R)
library(ggplot2)
library(biomaRt)
library(scales)
library(grid)
library(RColorBrewer)
library(xtable)

CVDproteins <- c("PSRC1","MIA3", "SMAD3","CDKN2A","CDKN2B","CXCL12","PDGFD","LIPA","PLTP","CETP","FMN2","HSPE1")

# Read in all protein text files based on names in CVDproteins
for (i in 1:length(CVDproteins)){

fname <- paste(CVDproteins[i],'txt',sep=".")
ptemp <- file.path('C://R-files//UKCI2016',fname) %>% read.delim(na.strings='')  # load it in.
ptemp <- ptemp[,c("X.node1","node2","combined_score")]  # just keep P1 and P2 and conf score
if(exists("PPI")){
  PPI <- rbind(PPI,ptemp)} 
else{
  PPI <- ptemp}
}


PPI<- get_interactions(CVDproteins)
head(PPI)

g <-graph.data.frame(ppi)       # convert to igraph object
ad <- get.adjacency(g)
nodesize=igraph::degree(g)# *5 # change the number if you want bigger nodes usually x 5
E(g)$weight=0.6
nodelabel<-V(g)$name

# normal genes are 'yellowgreen' and sized 5; hubs are 'tomato' color & sized 8"
# disease genes can be any colour but are always square in shape
# normal genes have degree of 1 - 19; hubs => 20 
for (i in 1:length(nodesize)) # NODESIZE
  if(nodesize[i] <= 30){nodesize[i] <- 5} else {nodesize[i] <- 8}

nodecolor=character(ncol(ad))  # NODECOLOR create a character for every column in adjaceny matrix
for (i in 1:length(nodecolor))
  if(stats[i,12] <= 0.8){nodecolor[i] <- "yellowgreen"} else {nodecolor[i]<-"tomato"}

nodeshape=character(ncol(ad)) # NODESHAPE create a character for every column in adjaceny matrix
for (i in 1:length(nodeshape))   # Look for our disease gene
  if(nodelabel[i] %in% CVDproteins){nodeshape[i] <- "square"} else {nodeshape[i]<-"circle"}

PPI<-ppi # OUT

g<-as.undirected(g); g<-simplify(g)
plot(g, edge.color="lightgray", vertex.color=nodecolor,vertex.label=nodelabel,vertex.shape=nodeshape,
     vertex.label.cex=0.6, vertex.label.font=1, 	vertex.frame.color="white",
     vertex.size=nodesize, vertex.label.color="black",edge.width=E(g)$weight, #*4
     layout=layout.kamada.kawai(g))   # layout.random; layout.kamada.kawai; 
                                            # layout.fruchterman.reingold; layout.graphopt; layout.spring 

# Now create a table for statistics on disease proteins based on "gs" structure
# hubness, authority, closeness betweenness
gs<-get_gstatistics(g)
head(gs)
gs$proteins <- rownames(gs)

stats<-gs

CVD <- gs[c(CVDproteins),c(12,9,10,13)]
CVD <- CVD[order(CVD$hubness,decreasing=TRUE),]
stats <- stats[order(stats$hubness,decreasing=TRUE),]
head(stats)



shite<-xtable(CVD,caption="Graph statistics measures for known disease proteins")#,digit=c(0,-1,-1,4,-1))
digits(shite) <- c(0,2,-2,3,4)
print(shite, include.rownames=TRUE)
      
shite<-xtable(stats[1:10,c(12,9,10,13)],caption="Graph statistics measures for top ranking proteins based on hubness")#,digit=c(1,-1,2,2,-2))
digits(shite) <- c(0,2,-2,3,4)
print(shite, include.rownames=TRUE)

xtable(stats[1,c(1:8)],caption="Overall graph statistics measures ")

# g is not fully interconnected, so count the associated nodes for each component
comps <- decompose.graph(g)
comps <- sapply(comps,vcount)
percentage <- sort((comps/vcount(g) * 100))
nodes<-sort(comps)
bg <- data.frame(numnodes=as.numeric(nodes),percentage)

# ------------------------------- FUNCTIONS ---------------------------------------

# Calculate some statistics about the disease gene network
get_gstatistics <- function(gt){
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


interactions_from_gene <- function(gene_id){                                                                                                                                                                         
  xmlrec <- entrez_fetch(db="gene", id=gene_id, rettype="xml", parsed=TRUE)                                                                                                                                          
  XML::xpathSApply(xmlrec,                                                                                                                                                                                           
                   "//Gene-commentary[Gene-commentary_heading[./text()='Interactions']]//Other-source[Other-source_src/Dbtag/Dbtag_db[./text()='GeneID']]//Other-source_anchor",
                   XML::xmlValue)                                                                                                                                  
}

# ----------- experimental stuff ---------------

# Detect sources of bias in reporting snps, keep count of published papers for each gene,
# that is to say the more researched that gene is it will have more SNPs. The rentrez package 
# can access ncbi (pubmed) database. The more recently discovered or less exciting genes will have
# fewer papers.

library(rentrez)

query <- "(McGarry K[Author]) AND sunderland[Affiliation]"
paper_search <- entrez_search(db="pubmed", query, retmax=100)
paper_search$count
kenspapers <- entrez_summary(db="pubmed", id=paper_search$ids)

date_and_cite <- extract_from_esummary(kenspapers, c("pubdate", "authors",  "title"))

what <- entrez_search(db="gene", term="(CDKN2A[GENE]) AND (Homo sapiens[ORGN])")
what$ids

# Start collating information for all proteins: 1.freq of snps, 2. no of publications 3. no of interactions
allproteins <-get.data.frame(g, what="vertices")

npubs <- count_articles(allproteins[1:length(allproteins)])
nsnps <- count_snps(allproteins[1:length(allproteins)])

# Problems connecting to ncbi server means we must send list of proteins in small batches when
# obtaining count of interactions. "Error in curl::curl_fetch_memory(url, handle = hand)", I'm not
# sure if its me overloading server or some bug elsewhere. A second problem is that some names dont have
# entrez  counterparts causing this error "Error in tryScan(getURL, sep = sep, quiet =". This
# forces mne to remove offending protein from allproteins list.
# I have Removed: Dvl1, HCVgp1, Dok1, Gsc, Cdc26, tat, EBNA-LP, NP_005216.1, HLA-F, RPSAP15, GCNT6  

# allproteins<-setdiff(allproteins,"HCVgp1") when you want to remove offending protein.

nacts25<- count_interactions(allproteins[1:25,1])
nacts50<- count_interactions(allproteins[26:50,1])
nacts75<- count_interactions(allproteins[51:75,1])
nacts100<- count_interactions(allproteins[76:100])
nacts125<- count_interactions(allproteins[101:125])
nacts150<- count_interactions(allproteins[126:150])
nacts175<- count_interactions(allproteins[151:175])
nacts200<- count_interactions(allproteins[176:200])
nacts225<- count_interactions(allproteins[201:225])
nacts250<- count_interactions(allproteins[226:250])

nacts275<- count_interactions(allproteins[251:275])
nacts300<- count_interactions(allproteins[276:300])
nacts399<- count_interactions(allproteins[301:399])

nacts400<- count_interactions(allproteins[400:499])

nacts500 <-count_interactions(allproteins[500:599])
nacts600 <-count_interactions(allproteins[600:699])
nacts700 <-count_interactions(allproteins[700:length(allproteins)]) 

nacts <- rbind(nacts25,nacts50,nacts75,nacts100,nacts125,nacts150,nacts175,nacts200,nacts225,nacts250,
               nacts275,nacts300,nacts399,nacts400,nacts500,nacts600,nacts700)

colnames(nacts) <- c("protein","nacts")
colnames(npubs) <- c("protein","npubs")
colnames(nsnps) <- c("protein","nsnps")
nacts<-as.data.frame(nacts)

alldata <- merge(nacts,npubs,by="protein")
alldata <- merge(alldata,nsnps,by="protein")

alldata$nacts<-as.numeric(levels(alldata$nacts))[alldata$nacts] 
alldata$npubs<-as.numeric(levels(alldata$npubs))[alldata$npubs]
alldata$nsnps<-as.numeric(levels(alldata$nsnps))[alldata$nsnps]


# now use ggplot2 for Fig 3.a
ggplot(alldata, aes(y = nacts, x = nsnps))  + 
          geom_point(color="blue",size=4,alpha=0.4) +
          theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=20,face="bold")) +
          #scale_size_continuous(range = c(2,4)) +
          scale_x_continuous(limits=c(0,20000)) +
          scale_y_continuous(limits=c(0,1000)) + 
          #coord_fixed(ratio=1)  +
          #geom_jitter(width = 1.0,height = 1.5) +
          labs(x="number of SNPS per protein", y = "number of interactions per protein")

# now use ggplot2 for Fig 3.b
ggplot(alldata, aes(y = nacts, x = npubs))  + 
        geom_point(color="red",size=4,alpha=0.4) + 
        theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=20,face="bold")) +
        scale_x_continuous(limits=c(0,2500)) +
        scale_y_continuous(limits=c(0,850)) + 
        #coord_fixed(ratio=1)  +
        #geom_jitter(width = 0.25,height = 0.5) +
        labs(x="number of publications per protein", y = "number of interactions per protein")

# now use ggplot2 for Fig 3.c
ggplot(alldata, aes(y = nsnps, x = npubs))  + 
  geom_point(color="tomato",size=4,alpha=0.4) + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=15,face="bold")) +
  scale_x_continuous(limits=c(0,2500)) +
  scale_y_continuous(limits=c(0,20000)) + 
  #coord_fixed(ratio=1)  +
  #geom_jitter(width = 0.25,height = 0.5) +
  labs(x="number of publications per protein", y = "number of SNPs per protein")


# ================================================================================
# Start of my functions......

# See how many research articles are written about our proteins. Uses rentrez package.
count_articles <- function (CVDP){
  for (i in 1:length(CVDP)){
    
    print(CVDP[i])
    pname <- paste(CVDP[i],'[GENE]) AND (Homo sapiens[ORGN])',sep="")
    ids<-entrez_search(db="pubmed", term=pname,retmax=40000)
    atemp <- cbind(CVDP[i],length(ids$ids))
    
    if(i!=1){
      articles <- rbind(articles,atemp)} 
    else{
      articles <- atemp}
  }
  return(articles)
}


#xmlrec <- entrez_fetch(db="gene", id="(TP53[GENE])AND (Homo sapiens[ORGN])", rettype="xml", parsed=TRUE) 
crap<-unique(interactions_from_gene(what$ids))

#------------------------------------

# For each disease protein get the proteins they interact with. Uses NCBI2R package.
get_interactions <- function(CVDP){

for (i in 1:length(CVDP)){
  print(CVDP[i])
  pname <- paste(CVDP[i],'[sym]',sep="")
  ids<-GetIDs(pname)
  plist<-GetInteractions(ids)
  plist<-unique(plist[13])
  
  cvd<-rep(CVDP[i],nrow(plist))
  
  ptemp <- cbind(cvd,plist)
  
  if(i!=1){
    ppi <- rbind(ppi,ptemp)} 
  else{
    ppi <- ptemp}
  }
  
  return(ppi)
}

# Count how many interactions each protein has.
count_interactions <- function(CVDP) {
  
      for (i in 1:length(CVDP)){
      
      print(CVDP[i])
      pname <- paste(CVDP[i],'[sym]',sep="")
      ids<-GetIDs(pname)
      plist<-GetInteractions(ids[1])
      plist<-unique(plist[13])
      
      ptemp <- cbind(CVDP[i],nrow(plist))
      
      if(i!=1){
        ppi <- rbind(ppi,ptemp)} 
      else{
        ppi <- ptemp}
    }
    
  return(ppi)
}

# Count hopw many SNPs are associated with each gene
count_snps <- function (CVDP){
  
    for (i in 1:length(CVDP)){
    print(CVDP[i])
    pname <- paste(CVDP[i],'[sym]',sep="")
    ids <- GetIDs(pname)
    snplist <- GetSNPsInGenes(ids[1],MaxRet=200000)
    ptemp <- cbind(CVDP[i],length(snplist))
    
    if(i!=1){
      snps <- rbind(snps,ptemp)} 
    else{
      snps <- ptemp}
    
  }
  return(snps)
}

fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[2]
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
    theme(axis.text.x=element_text(size=7,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=7,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=8,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=8,color=color.axis.title, vjust=1.25)) +
    
    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}


#--- TKPLOT allows you drag nodes around and create a better graph --------
tkplot(g,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,
       edge.arrow.size=0, edge.curved=FALSE)







# kirsty1.R
# 22/04/2016
# Code developed by Kirsty Emery and Ken McGarry to analysis differences between wildtype and mutant
# source("http://bioconductor.org/biocLite.R")
# biocLite("dplyr")
library(dplyr)
library(tidyr)
library(igraph)
library(Peptides)
library(seqinr)  # there is a conflict with a getSequence() in biomaRt package
library(NCBI2R) # is removed from CRAN but archived vesrion still works.
library(biomaRt) 
library(stringr)
library(protr)
library(Biostrings)

CVDproteins <- c("PSRC1","MIA3", "SMAD3", "CDKN2A","CDKN2B","CXCL12","PDGFD","LIPA","PDGFD","PLTP","CETP","FMN2","HSPE1")

mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

ENSG <- getBM(mart=mart, attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position"),
              filters="hgnc_symbol", values=CVDproteins)

snpmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snpSMAD3 <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
                   filters="ensembl_gene", values=ENSG[12,2],checkFilters=FALSE)
patho <- snpSMAD3[grep("pathogenic", snpSMAD3[,6]), ] # get harmful snps
patho <- patho[unique(patho[1,])]
patho <- patho[!duplicated(patho[,1]),]
head(patho)

alle <- getBM(attributes=c("snp","refsnp_id","chrom_start","allele","chr_name","ensembl_transcript_stable_id"),
                    filters=c("snp_filter","downstream_flank","upstream_flank"),
                    values=list(patho[,1],238,238),mart=snpmart, checkFilters=FALSE)





snpCDKN2A <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
                   filters="ensembl_gene", values=ENSG[1,2],checkFilters=FALSE)
coding <- snpCDKN2A[grep("coding_sequence_variant", snpCDKN2A[,6]), ] #  "consequence_type_tv" variable
head(coding)
patho <- snpCDKN2A[grep("pathogenic", snpCDKN2A[,6]), ] #  "consequence_type_tv" variable
head(patho)


alleCDKN2A <- getBM(attributes=c("snp","refsnp_id","chrom_start","allele","chr_name","ensembl_transcript_stable_id"),
              filters=c("snp_filter","downstream_flank","upstream_flank"),
              values=list(snpCDKN2A[1:10,1],238,238),mart=snpmart, checkFilters=FALSE)

process_snippet <- function(thesnippet){
  
  thesnippet <- str_replace_all(thesnippet, "[%/]", "") # remove percents and backslash.
  thesnippet<-seqinr::s2c(thesnippet)
  n<-length(thesnippet)
  
  posW <- n/2    # the wildtype position
  posM <- 1+(n/2)    # the mutation position
  
  downstream <- thesnippet[1:((n/2)-1)]  # to the left of the position
  upstream <- thesnippet[((n/2)+2):n] # to right of the position
  
  wildtype <- c(downstream,thesnippet[posW],upstream)
  mutant <- c(downstream,thesnippet[posM],upstream)
  
  snippet_data <- data.frame(wildtype="EMPTY",mutant="EMPTY")
  snippet_data$wildtype <- c2s(wildtype)
  snippet_data$mutant <- c2s(mutant)
  
  return(snippet_data)
}

CDKN2A1 <- process_snippet(alleCDKN2A[1,1])
CDKN2A2<- process_snippet(alleCDKN2A[2,1])

CDKN2A1w<- (CDKN2A1$wildtype)
CDKN2A1m<-(CDKN2A1$mutant)

CDKN2A2w<- (CDKN2A2$wildtype)
CDKN2A2m<-(CDKN2A2$mutant)

#Convert nucleotide sequence to peptide string
Conversion <- function(NC) {
  NUC<-s2c(NC)
  
  AA<-translate(NUC, frame=0, sens = "F", numcode=1, NAstring="X", ambiguous= "FALSE")
  
  peptide<-toString(AA)
  
  return(peptide)
  
}

#Use of Conversion function
CDKN2A1wAA<-Conversion(CDKN2A1w)
CDKN2A1mAA<-Conversion(CDKN2A1m)

CDKN2A2wAA<-Conversion(CDKN2A2w)
CDKN2A2mAA<-Conversion(CDKN2A2m)


#Calculate physiochemical properties of peptide
AAproperties<-function(peptide) {
  
  aindex<-aindex(peptide)
  
  boman<-boman(peptide)
  
  charge<-charge(peptide)
  
  instaindex<-instaindex(peptide)
  
  lengthpep<-lengthpep(peptide)
  
  p_values<-c(aindex,boman,charge,instaindex,lengthpep)
  return(p_values)
  
}

#Use of AAproperties function
PropertiesCDKN2a1w<-AAproperties(CDKN2A1wAA)
PropertiesCDKN2a1m<-AAproperties(CDKN2A1mAA)

PropertiesCDKN2a2w<-AAproperties(CDKN2A2wAA)
PropertiesCDKN2a2m<-AAproperties(CDKN2A2mAA)

#Collate together to make single list
PropertiesCDKN2A1<-c(PropertiesCDKN2a1w,PropertiesCDKN2a1m)
PropertiesCDKN2A2<-c(PropertiesCDKN2a2w,PropertiesCDKN2a2m)

#Make table comparing two strains
Table<-function(properties) {
  trial<-matrix(c(properties),ncol=2)
  colnames(trial)<-c('Wildtype','Mutant')
  rownames(trial)<-c('aindex','boman','charge','instaindex','lengthpep' )
  trial.table<-as.table(trial)
  
  return(trial.table)
}

#Use table function
physiochemicalCKDN2A1<-Table(PropertiesCDKN2A1)
physiochemicalCKDN2A2<-Table(PropertiesCDKN2A2)

#Compute the amino acid composition of a protein sequence and turn into a data frame
datachange<-function(data) {
  aacomp<-aacomp(data)
  frame<-data.frame(aacomp)
  
  return(frame)
}

#Use datachange function
dataCDKN2A1<-datachange(CDKN2A1wAA)
dataCDKN2A12<-datachange(CDKN2A1mAA)

dataCDKN2A2<-datachange(CDKN2A2wAA)
dataCDKN2A22<-datachange(CDKN2A2mAA)

#So we have names for the rows when we merge the two tables together
rownames<-c("Tiny","Small","Alipathic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")

#Convert row names to a data frame
rows<-data.frame(rownames)

#Collate all 3 tables together
TableCDKN2A1<-c(rows,dataCDKN2A1,dataCDKN2A12)
TableCDKN2A2<-c(rows,dataCDKN2A2,dataCDKN2A22)

#Show Table for physiochemical properties
physiochemicalCKDN2A1
physiochemicalCKDN2A2

#Print table of amino acid composition
data.frame(TableCDKN2A1)
data.frame(TableCDKN2A2)


#For next Protein (CDKN2B) (Shortened, without functions etc):
snpCDKN2B <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
                   filters="ensembl_gene", values=ENSG[2,2],checkFilters=FALSE)

coding1 <- snpCDKN2B[grep("coding_sequence_variant", snpCDKN2B[,6]), ] #  "consequence_type_tv" variable
head(coding1)

patho1 <- snpCDKN2B[grep("pathogenic", snpCDKN2B[,6]), ] #  "consequence_type_tv" variable
head(patho1)

alleCDKN2B <- getBM(attributes=c("snp","refsnp_id","chrom_start","allele","chr_name","ensembl_transcript_stable_id"),
                    filters=c("snp_filter","downstream_flank","upstream_flank"),
                    values=list(snpCDKN2B[1:10,1],238,238),mart=snpmart, checkFilters=FALSE)

CDKN2B1 <- process_snippet(alleCDKN2B[1,1])
CDKN2B2<- process_snippet(alleCDKN2B[2,1])
CDKN2B3<- process_snippet(alleCDKN2B[3,1])

CDKN2B1w<- (CDKN2A1$wildtype)
CDKN2B1m<-(CDKN2A1$mutant)

CDKN2B2w<- (CDKN2A2$wildtype)
CDKN2B2m<-(CDKN2A2$mutant)

CDKN2B3w<-(CDKN2B3$wildtype)
CDKN2B3m<-(CDKN2B3$mutant) 

#Use of Conversion function
CDKN2B1wAA<-Conversion(CDKN2B1w)
CDKN2B1mAA<-Conversion(CDKN2B1m)

CDKN2B2wAA<-Conversion(CDKN2B2w)
CDKN2B2mAA<-Conversion(CDKN2B2m)

CDKN2B3wAA<-Conversion(CDKN2B3w)
CDKN2B3mAA<-Conversion(CDKN2B3m)

#Use of AAproperties function
PropertiesCDKN2B1w<-AAproperties(CDKN2B1wAA)
PropertiesCDKN2B1m<-AAproperties(CDKN2B1mAA)

PropertiesCDKN2B2w<-AAproperties(CDKN2B2wAA)
PropertiesCDKN2B2m<-AAproperties(CDKN2B2mAA)

PropertiesCDKN2B3w<-AAproperties(CDKN2B3wAA)
PropertiesCDKN2B3m<-AAproperties(CDKN2B3mAA)

#Collate together to make single list
PropertiesCDKN2B1<-c(PropertiesCDKN2B1w,PropertiesCDKN2B1m)
PropertiesCDKN2B2<-c(PropertiesCDKN2B2w,PropertiesCDKN2B2m)
PropertiesCDKN2B3<-c(PropertiesCDKN2B3w,PropertiesCDKN2B3m)

#Use table function
physiochemicalCKDN2B1<-Table(PropertiesCDKN2B1)
physiochemicalCKDN2B2<-Table(PropertiesCDKN2B2)
physiochemicalCKDN2B3<-Table(PropertiesCDKN2B3)

#Use datachange function
dataCDKN2B1<-datachange(CDKN2B1wAA)
dataCDKN2B12<-datachange(CDKN2B1mAA)

dataCDKN2B2<-datachange(CDKN2B2wAA)
dataCDKN2B22<-datachange(CDKN2B2mAA)

dataCDKN2B3<-datachange(CDKN2B3wAA)
dataCDKN2B32<-datachange(CDKN2B3mAA)

#Collate all 3 tables together
TableCDKN2B1<-c(rows,dataCDKN2B1,dataCDKN2B12)
TableCDKN2B2<-c(rows,dataCDKN2B2,dataCDKN2B22)
TableCDKN2B3<-c(rows,dataCDKN2B3,dataCDKN2B32)

#Show Table for physiochemical properties
physiochemicalCKDN2B1
physiochemicalCKDN2B2
physiochemicalCKDN2B3

#Print table of amino acid composition
data.frame(TableCDKN2B1)
data.frame(TableCDKN2B2)
data.frame(TableCDKN2B3)

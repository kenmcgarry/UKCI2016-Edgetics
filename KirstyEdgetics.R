# KirstyEdgetics.R
# 1st stage in Edgetics: 
#     1. get protein cDNA sequences (key proteins invovled in the disease)
#     2. get snp alleles and modify sequence
#     3. use peptides package to generate statistics on modified protein sequences
#     4. detect if protein characteristics are modified (many mutations wont)
#     5. express changed protein as a complex network indicating changes to interactions with partners
#               source("http://bioconductor.org/biocLite.R")
#               biocLite("dplyr")

library(dplyr)
library(tidyr)
library(igraph)
library(Peptides)
library(ppiPre)
library(seqinr)  # there is a conflict with a getSequence() in biomaRt package
library(NCBI2R) # is removed from CRAN but archived vesrion still works.
library(biomaRt) 
library(stringr)

CVDproteins <- c("PSRC1","MIA3", "SMAD3","CDKN2A","CDKN2B","CXCL12","PDGFD","LIPA","PLTP","CETP","FMN2","HSPE1")

# obtain details of human proteins database and save in "mart", then get CDKN2A protein ensembl ID and save in "ENSG"
mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ENSG <- getBM(mart=mart, attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position"),
              filters="hgnc_symbol", values=CVDproteins)

# SNP --- setup database calls
snpmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

# Now obtain a list of SNP's attributed to proteins of interest, open up hsapiens_snp database & save to "mart" 
# and download the id's of the SNPs attributed. 
snp <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
             filters="ensembl_gene", values=ENSG[12,2],checkFilters=FALSE)

# search for coding mutations
coding <- snp[grep("coding_sequence_variant", snp[,6]), ] #  "consequence_type_tv" variable
head(coding)

patho <- snp[grep("pathogenic", snp[,6]), ] #  "consequence_type_tv" variable
head(patho)


# alle dataframe contains the snp and other information.
# Just be sure when changing flank sizes to keep both same and both even numbers 10,20,40
alle <- getBM(attributes=c("snp","refsnp_id","chrom_start","allele","chr_name","ensembl_transcript_stable_id"),
          filters=c("snp_filter","downstream_flank","upstream_flank"),
          values=list(snp[1:10,1],10,10),mart=snpmart, checkFilters=FALSE)


somedata <- process_snippet(alle[4,1]) # snp no 4, has AA changed!

dna <-(s2c(as.character(somedata[1]))) # "1" for wild type "2" for mutant
AA <- translate(dna)
aacomp(seq = c2s(AA))

# A function to determine if Amino Acid has changed thanks to SNP



# Function to create two snippets - one is the wild type the other the snp mutant
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




# dna2protein.R
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
#library(protr)
library(seqinr)  # there is a conflict with a getSequence() in biomaRt package
#library(Interpol)
library(NCBI2R) # is removed from CRAN but archived vesrion still works.
library(biomaRt) 
library(Biostrings)
data(BLOSUM50)

CVDproteins <- c("PSRC1","MIA3", "SMAD3", "CDKN2A","CDKN2B","CXCL12","PDGFD","LIPA","PDGFD","PLTP","CETP","FMN2","HSPE1")

# obtain details of human proteins database and save in "mart", then get CDKN2A protein ensembl ID and save in "ENSG"
mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ENSG <- getBM(mart=mart, attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position"),
              filters="hgnc_symbol", values=CVDproteins)

# Now obtain a list of SNP's attributed to proteins of interest, open up hsapiens_snp database & save to "mart" 
# and download the id's of the SNPs attributed. 
snpmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

seqlist <- biomaRt::getSequence(id="SMAD3",type="hgnc_symbol",seqType="gene_exon", mart=mart)
seq1 <- s2c(seqlist[1,1])

# CDKN2A is located from = 21967753 to 21995301 which gives 27548 BP's, 
# downloading gene_exon_intron gives 27549 BP's, we get a slight difference.
# Thats the full coding and noncoding sequence i.e. everything, I need to download only those
# coding sequences, modify with snp allele. Look into frequency of snp's i.e. how many
# can occur for a single person? Check frequency of snp occurences.
# Now renumber snp locations to match 1 - 27548

length(seq1)
x <- seqinr::translate(seq=seq1)

# Just to confirm that our DNA to Amino conversion has worked download the actual AA sequence & compare
testlist <- biomaRt::getSequence(id="SMAD3",type="hgnc_symbol",seqType="peptide", mart=mart)
testseq1 <- s2c(testlist[1,1])


# SNP --- downloads

snp <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
             filters="ensembl_gene", values=ENSG[12,2],checkFilters=FALSE)
# search for coding mutations
coding <- snp[grep("coding_sequence_variant", snp[,8]), ] #  "consequence_type_tv" variable
head(coding)

patho <- snp[grep("pathogenic", snp[,6]), ] #  "consequence_type_tv" variable
head(patho)


#snpseq<-getBM(attributes=c("refsnp_id","snp"),
#              filters=c("snp_filter","downstream_flank","upstream_flank"),
#              values=list(c("rs797015820","rs372649305","rs1800586"),10,10),
#              mart=snpmart, checkFilters=FALSE)

alle <- getBM(attributes=c("snp","refsnp_id","chrom_start","allele","chr_name","ensembl_transcript_stable_id"),
          filters=c("snp_filter","downstream_flank","upstream_flank"),
          values=list(snp[1:10,1],10,10),mart=snpmart, checkFilters=FALSE)




# This part uses NCBI2R package, GetSNPsInGenes() will return list of snps associated with gene of interest
#
singlegene<-GetIDs("SMAD3[sym]") 
snplist<- GetSNPsInGenes(singlegene, batchsize=50, MaxRet = 30000, showurl = FALSE,quiet = TRUE, smt = FALSE, sme = FALSE)
snpinfo <-GetSNPInfo(snplist)

#singlegene<-GetIDs("BCL2L1[sym] human")
#ggi<-GetGeneInfo(singlegene)
#ppi<-GetInteractions(singlegene)

#### LOAD FUNCTIONS #######

process_snippets <- function(allelist){
  N <- nrow(allelist)
  snippet_data <- data.frame()
  z <- s2c(alle[1,1])
  z <- data.frame(sapply(z, function(x) gsub("%", "", x)))

  
  return(snippet_data)
}


# Extract name of 1st ENST from a list and return it with its sequence.
get_ENST <- function(ENST){
  if(is.null(ENST)){  # error trapping for calls without ENST's!
    return(NULL)}
  
  # debug: ENST
  ENST <- strsplit(as.character(ENST), ";")
  
  ENST <- 
    
  enstmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  # debug:transcript <- biomaRt::getSequence(id="ENST00000200676",type="ensembl_transcript_id",seqType="cdna", mart=enstmart)
  
  return(ENST)
}

dotPlot(s2c(alle[1,1]), s2c(transcript[1,1]))
transcript1 = biomaRt::getSequence(id="ENST00000327367",type="ensembl_transcript_id",seqType="coding",mart=enstmart)
transcript2 = biomaRt::getSequence(id="ENST00000540846",type="ensembl_transcript_id",seqType="coding",mart=enstmart)
transcript3 = biomaRt::getSequence(id="ENST00000439724",type="ensembl_transcript_id",seqType="coding",mart=enstmart)

pep = biomaRt::getSequence(id="ENST00000200676",type="ensembl_transcript_id",seqType="peptide",mart=enstmart)

x <- "GCTGCTAGGGAATCCAGATGG"

globalAlign <- pairwiseAlignment((x),
                                 (transcript3[1,1]),
                                substitutionMatrix = BLOSUM50, gapOpening = -2, 
                                gapExtension = -8, scoreOnly = FALSE,type="local")

# gene_exon, transcript_exon, transcript_exon_intron, 
# gene_exon_intron, cdna, coding, coding_transcript_flank, 
# coding_gene_flank, transcript_flank, gene_flank, 
# peptide, 3utr or 5utr.

y<-seqinr::translate(seq=s2c(transcript1[1,1]))  # translate the dna of the transcript into its peptide (amino acids)
                                        # and it should be identical to pep
z<-seqinr::s2c(transcript1[1,1])


# DO A HARD CODED EXAMPLE FIRST AS PROOF OF PRINCIPLE.
# i.e. SMAD gene.
# hgnc  ensembl           Chrom    start        stop       BP
# SMAD3 ENSG00000166949   15       67063763     67195195   131432

testlist <- biomaRt::getSequence(id="SMAD3",type="hgnc_symbol",seqType="peptide", mart=mart)
testseq1 <- s2c(testlist[1,1])
smad3 = biomaRt::getSequence(start=67063763, end=67195195, chromosome=15,type="entrezgene",
                                   seqType="_exon",mart=ensembl)
smad3dna<-seqinr::s2c(smad3[1,1])
smad3pep<-seqinr::translate(seq=smad3dna) 

# Get SNPs for this gene
snp <- getBM(mart=snpmart, attributes=c("refsnp_id","associated_gene","minor_allele_freq","allele","chr_name","clinical_significance","chrom_start","consequence_type_tv"), 
             filters="ensembl_gene", values=ENSG[12,2],checkFilters=FALSE)


# GET THE TRANSCRIPT VERSION AS PROOF OF PRINCIPLE.
enstmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ensembl <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
snpmart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
snpmart <- useMart("snp", dataset="hsapiens_snp")

transcript1 = biomaRt::getSequence(id="ENST00000327367",type="ensembl_transcript_id",seqType="coding",mart=enstmart)
transcript2 = biomaRt::getSequence(id="ENST00000540846",type="ensembl_transcript_id",seqType="coding",mart=enstmart)
transcript3 = biomaRt::getSequence(id="ENST00000439724",type="ensembl_transcript_id",seqType="coding",mart=enstmart)

transinfo <- getBM(attributes=c("ensembl_transcript_id","transcript_start","transcript_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","chromosome_name","transcript_biotype"),
            filters ="ensembl_transcript_id", values="ENST00000327367", mart=ensembl)

transcript1dna<-seqinr::s2c(transcript1[1,1]) # "ENST00000327367"
transcript1pep<-seqinr::translate(seq=transcript1dna) 

# NOW FIND SNPS WITHIN RANGE OF TRANSCRIPT1 ("ENST00000327367")
# start           end
# 67065845        67195195
snpstranscript <- getBM(attributes=c("refsnp_id","allele","chrom_start","chrom_strand"),
                        filters = c("chr_name","start","end"),
                        values = list(15,67065845,67195195), mart = snpmart)

snpstranscript <- getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_strand'),
                        filters = c('chromosomal_region'),
                        values = "15:67065845:67195195",mart = snpmart)



data = data.frame(chr = "chr17", start = 63973115, end = 64437414)
data$query = paste(gsub("chr",'',data$chr),data$start,data$end,sep=":")
out = getBM(attributes = c('refsnp_id','allele','chrom_start'), 
  filters = 'chromosomal_region', values = data$query, mart = snpmart)


regions = list( c("1:661517:661671", "1:787463:787591"))
attribs = (c("refsnp_id","allele","chrom_start","chrom_strand"))


filters <- listFilters(snpmart)
filters[grep("start|end|strand", filters$name),]


result = getBM(filters="chromosomal_region", values=regions, attributes=attribs, mart=snpmart)


      getBM(c('refsnp_id','allele','chrom_start','chrom_strand'), 
          filters = c('chr_name','chrom_start','chrom_end'), 
          values = list(8,148350,148612), mart = snpmart)
    
      
     
result <- getBM(c('refsnp_id','allele','chrom_start','chrom_strand',
                'consequence_type_tv'), filters = c('chr_name','chrom_start', 'chrom_end'), 
              values = list(15,67065845,67195195),mart=snpmart)



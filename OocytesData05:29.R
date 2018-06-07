
#PART 1:Setting up the data.table
library(Biobase)
library(GEOquery)
library(AnnotationDbi)
library(mouse4302.db)
library(gplots)
library(pathVar)
library(KEGGREST)
library(mouse4302.db)
library(limma)
source("https://bioconductor.org/biocLite.R")
library(cogena)


a1 <- getGEO("GSE11667")
a1 <- a1[[1]] # just make expression set easier to access

dat.a1 <- exprs(a1)
pdat.a1 <- pData(a1)

#assigning colnames based on the pdata 
colnames(dat.a1) <- pdat.a1$title

#export the file from R into txt file
write.table(dat.a1,file="GSE11667_epression_data.txt",quote=F,sep="\t",col.names=T,row.names=F)
write.table(pdat.a1,file="GSE11667_phenotype_data.txt",quote=F,sep="\t",col.names=T,row.names=F)

#converting probes to genes  
probe_to_symbol <- select(mouse4302.db,keys=rownames(dat.a1),keytype="PROBEID",columns=c("PROBEID","SYMBOL"))
rownames(probe_to_symbol) <- probe_to_symbol$PROBEID
probe_to_symbol <- probe_to_symbol[!is.na(probe_to_symbol[,2]),]

# create an empty variable. This will ultimately be our new matrix with gene symbols instead of probe IDs
dat.a1_symbol <- NULL
# do a for loop where each iteration, we get data for a different gene symbol
for(i in 1:length(unique(probe_to_symbol$SYMBOL))) {
  uniSymbol <- unique(probe_to_symbol$SYMBOL)[i] # assign uniSymbol to the gene symbol I want
  # get all the probe IDs assigned to our gene symbol
  probeNames <- probe_to_symbol[which(uniSymbol==probe_to_symbol$SYMBOL),"PROBEID"]

  # get the data that goes with the probe IDs
  temp_dat <- dat.a1[probeNames,]
  # if we have more than 1 probe ID for a single gene symbol, we take the average of the rows.
  if( length(probeNames) > 1 ){
    temp_dat <- colMeans(temp_dat)
  }
  # add our newly created row to our final expression matrix
  dat.a1_symbol <- rbind(dat.a1_symbol, temp_dat)
}
# set the rownames to be the gene symbols
rownames(dat.a1_symbol) <- unique(probe_to_symbol$SYMBOL)
View(dat.a1_symbol)

#Project: assessing the effect of maternal age on GV and MII oocyte gene expression

#Part 1:comparing the avg number of transcripts in aged &young mice in GV and MII stages

#1a) subset into separate data frames (young, old, GV, MII)

dat.a1_symbol<-data.matrix(dat.a1_symbol)
#Aged
cn<- colnames(dat.a1_symbol)
agedcol<-cn[grep("aged", cn)]
aged<-dat.a1_symbol[,agedcol]
#Young
youngcol<-cn[grep("young", cn)]
young<-dat.a1_symbol[,youngcol]
#GV Stage
GVcol<-cn[grep("GV", cn)]
GV<-dat.a1_symbol[,GVcol]
#MII Stage
MIIcol<-cn[grep("MII", cn)]
MII<-dat.a1_symbol[,MIIcol]


#Alternatively: can also write a function to do this, as written below

grep_columns<- function(column_name) {
grepped_column_names <- cn[grep(column_name, cn)]
grepped<-dat.a1_symbol[, grepped_column_names]
return(grepped)
}
aged<- grep_columns("aged")
young<- grep_columns("young")
GV<- grep_columns("GV")
MII<- grep_columns("MII")

# 1b)Does the avg expression drastically differ between young & old GV ( oocytes), young & old MII (eggs)?
#find the grand mean of the replicates for each stage, and compare them

agedGV<- GV[,1:4]
agedGVmean<- rowMeans(agedGV)

youngGV<- GV[, 5:8]
youngGVmean<- rowMeans(youngGV)

agedMII<- MII[,1:4]
agedMIImean<- rowMeans(agedMII)

youngMII<- MII[, 5:8]
youngMIImean<- rowMeans(youngMII)

grandmeans<- cbind(youngGVmean, agedGVmean, youngMIImean, agedMIImean)

#log transform all matrices
logdata<-log(dat.a1_symbol+1, 2) 
logyoungGV<- log(youngGV+1, 2)
logagedGV<-log(agedGV+1, 2)
logGV<-log(GV+1,2)
logyoungMII<-log(youngMII+1, 2)
logagedMII<- log(agedMII+1, 2)
logMII<-log(MII+1,2)

#observation: we discover that there are similar numbers of expression values detected in young & old oocytes, and young & old eggs. fewer transcripts seem to be detected in eggs than in oocytes

#PART 2) two-way unsupervised hierachical cluster to see how the replicates independently cluster


colnames(dat.a1_symbol)<-c("Aged GV replicate 1", "Aged GV replicate 2", "Aged GV replicate 3", "Aged GV replicate 4", "young GV replicate 1", "young GV replicate 2", "young GV replicate 3", "young GV replicate 4", "aged MII replicate 1","aged MII replicate 2", "aged MII replicate 3", "aged MII replicate 4", "young MII replicate 1", "young MII replicate 2", "young MII replicate 3", "young MII replicate 4")
corr.dist = function(x) as.dist(1-cor(t(x))) #pearson correlation as our distance metric
clustered<- heatmap(dat.a1_symbol, RowV=TRUE, ColV=TRUE, margins=c(10,1),labRow=FALSE, rowdistfun=corr.dist, scaled="row")

#observation: all replicates independently clustered according to their age, with branch distances fairly minimal between replicates.

###############################################################################################

#PART 3) DIFFERENTIAL GENE ANALYSIS Using T Test and Limma

#3a) T test for young versus old GV mice 

pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(dat.a1_symbol)) { # For each gene : 
  x = logyoungGV[i,] # young GV expression of gene number i
  y = logagedGV[i,] # aged GV expression of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

sum(pvalue<0.05)   
GVdeg <- p.adjust(pvalue, "BH") < 0.05 # adjust for multiple corrections and set p value cutoff 
GVgenes<- row.names(dat.a1_symbol)[which(GVdeg)] 

#observation: 19 GV genes w/ p<0.05

#3b) T test for young versus old MII mice 

pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(logdata)) { # For each gene : 
  x = logyoungMII[i,] # young GV expression of gene number i
  y = logagedMII[i,] # aged GV expression of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
sum(pvalue<0.05)   
MIIdeg <- (p.adjust(pvalue, "BH"))< 0.05   # adjusted for multiple corrections and filter P value
sum(MIIdeg) 
MIIgenes<- row.names(dat.a1_symbol)[which(MIIdeg)] 

#observation: 1120 MII genes with p<0.05

#3c) Venn diagram to compare the overlap of differentially expressed GV and MII genes as found through the T test

Ttestvenn<- venn(list(A=GVgenes,B=MIIgenes))

#3d) Using Limma w/ Contrast to identify differentially expressed genes:

#Differentially expressed GV genes ( young versus aged)
design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2,2)))
colnames(design) <- c("GV_Aged","GV_Young")
contrast.matrix <- makeContrasts(GV_Aged-GV_Young, levels =design)
fit<- lmFit(logGV, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
LimmaGVgenes<-topTable(fit2, adjust.method="BH", p.value=0.05, number=Inf)

##Differentially expressed MII genes ( young versus aged)
design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2,2)))
colnames(design) <- c("MII_Aged","MII_Young")
contrast.matrix <- makeContrasts(MII_Aged-MII_Young, levels =design)
fit<- lmFit(logMII, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
LimmaMIIgenes<-topTable(fit2, adjust.method="BH", p.value=0.05, number=Inf) 

#CONCLUSION: there are many more genes that are differentially expressed in the MII (egg) stage compared to GV(Oocyte) stage. Thus, the effect of maternal age appears to have a much more pronounced effect on the transcript profile in eggs than in oocytes.

#3f) Venn Diagram to show the overlap of differentially expressed genes found from the T test versus Limma for GV and MII stage

#GV stage Venn Diagram (diff. expressed genes from T test versus Limma )
LimmaGV <- cbind(genes=rownames(LimmaGVgenes),LimmaGVgenes) # add a column to the LimmaGV dataframe containing the 13 sig genes
LimmaGV<-as.character(LimmaGV[,1])# subset the data frame to show only that column, and convert to a character vector ( same  class as GVgenes)
GVvenn<- venn(list(A=GVgenes,B=LimmaGV))

#MII stage venn diagram  (diff. expressed genes from T test versus Limma )
LimmaMII <- cbind(genes=rownames(LimmaMIIgenes),LimmaMIIgenes) # add a column to the LimmaMII dataframe containing the 698 sig genes
LimmaMII<-as.character(LimmaMII[,1])# subset the data frame to show only that column, and convert to a character vector ( same  class as MIIgenes)
MIIvenn<- venn(list(A=MIIgenes,B=LimmaMII))

#############################################################################################

#PATHVAR for GV and MII pathway analysis using Reactome, Kegg, and Hallmark 
#part 1- GV stage using reactome ( via density and distribution counts)
#part 2- GV stage using kegg (via distribution counts)
#part 3- MII stage using kegg (via distribution counts)
#part 4- GV stage using hallmark (via distribution counts)
#part 5- MII stage using hallmark(via distribution counts)


#code from Sam to convert human reactome/kegg to mouse 

pways.kegg.mouse=NULL
pways.kegg.mouse$PATHID <- names(keggList("pathway", "mmu"))
pway.name.kegg.mouse <- keggList("pathway","mmu")
pways.kegg.mouse$PATHNAME=lapply(strsplit(pway.name.kegg.mouse,split=" -"), function(x) x[1])
pways.kegg.mouse$GENES=list()
for( i in 1:length(pways.kegg.mouse$PATHID) ){
  geneSymbol <- keggGet(pways.kegg.mouse$PATHID[[i]])[[1]]$GENE[c(FALSE,TRUE)]
  if(is.null(geneSymbol)) {
    pways.kegg.mouse$GENES[[i]] <- as.character()
  } else {
    geneSymbol <- sapply(strsplit(geneSymbol,split=";"), function(x) x[1])
    pways.kegg.mouse$GENES[[i]] <- unique(geneSymbol)
  }
}
names(pways.kegg.mouse$GENES) <- pways.kegg.mouse$PATHNAME
pways.kegg.mouse$SIZE <- as.numeric(lapply(pways.kegg.mouse$GENES, length))
pways.kegg.mouse$PATHNAME <- unlist(pways.kegg.mouse$PATHNAME)

# remove pathways that do not have any genes in them
pways.kegg.mouse$PATHID <- pways.kegg.mouse$PATHID[!pways.kegg.mouse$SIZE==0]
pways.kegg.mouse$PATHNAME <- pways.kegg.mouse$PATHNAME[!pways.kegg.mouse$SIZE==0]
pways.kegg.mouse$GENES <- pways.kegg.mouse$GENES[!pways.kegg.mouse$SIZE==0]
pways.kegg.mouse$SIZE <- pways.kegg.mouse$SIZE[!pways.kegg.mouse$SIZE==0]

# pathways below do not have any genes in them for mouse
#[1] "Metabolic pathways"                        "Carbon metabolism"
#[3] "2-Oxocarboxylic acid metabolism"           "Fatty acid metabolism"
#[5] "Biosynthesis of amino acids"               "EGFR tyrosine kinase inhibitor resistance"
#[7] "Endocrine resistance"                      "Antifolate resistance"
#[9] "Platinum drug resistance"

# now make mouse reactome pathways
library(AnnotationDbi)
reactome_db <- read.table("https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",comment="",sep="\t",stringsAsFactors=F,quote="") # May 23, 2018
reactome_db_mouse <- reactome_db[reactome_db[,6]=="Mus musculus",]
pways.reactome.mouse=NULL
pways.reactome.mouse$PATHID <-  unique(reactome_db_mouse[,2])
pways.reactome.mouse$PATHNAME = reactome_db_mouse[match(pways.reactome.mouse$PATHID,reactome_db_mouse[,2]),4]
pways.reactome.mouse_entrez <- lapply(split(reactome_db_mouse,reactome_db_mouse[,4]), function(x) x[,1])
pways.reactome.mouse$GENES=list()
for(i in 1:length(pways.reactome.mouse_entrez)) {
  myEntrezIDs <- pways.reactome.mouse_entrez[[i]]
  entrezToSym <- try(select(mouse4302.db, keys=myEntrezIDs,keytype="ENTREZID",columns=c("ENTREZID","SYMBOL")))
  if(class(entrezToSym) == "try-error") { # this happens if no probes match to a gene in the pathway
    entrezToSym <- as.character()
    pways.reactome.mouse$GENES[[i]] <- entrezToSym
  } else {
    mySyms <- unique(entrezToSym$SYMBOL)
    pways.reactome.mouse$GENES[[i]] <- mySyms
  }
}
names(pways.reactome.mouse$GENES) <- names(pways.reactome.mouse_entrez)
pways.reactome.mouse$GENES <- pways.reactome.mouse$GENES[pways.reactome.mouse$PATHNAME]
pways.reactome.mouse$SIZE <- as.numeric(lapply(pways.reactome.mouse$GENES, length))
# remove reactome pathways that do not have any genes in them. Only 1 reactome pathway was empty.
pways.reactome.mouse$PATHID <- pways.reactome.mouse$PATHID[!pways.reactome.mouse$SIZE==0]
pways.reactome.mouse$PATHNAME <- pways.reactome.mouse$PATHNAME[!pways.reactome.mouse$SIZE==0]
pways.reactome.mouse$GENES <- pways.reactome.mouse$GENES[!pways.reactome.mouse$SIZE==0]
pways.reactome.mouse$SIZE <- pways.reactome.mouse$SIZE[!pways.reactome.mouse$SIZE==0]
################################################################################################################

#Part 1: GV Stage Using Reactome Pathways: 

# Identify our two groups: GV aged versus GV young) 
samp.id <- colnames(logGV)
cell.id <- character(length(samp.id))
cell.id[grep("aged", samp.id)] <- "aged"
cell.id[grep("young", samp.id)] <- "young"

#filter lowly expressed genes (keep the transcripts that have expression >1 in at least 75% of the replicates)
qc.aged <- apply(logGV[,cell.id == "aged"], 1, function(x,ct){ sum(x >= ct) }, ct=1)
qc.young <- apply(logGV[,cell.id == "young"], 1, function(x,ct){ sum(x >= ct) }, ct=1)
logGV.age_young <- logGV[qc.aged >= .75*sum(cell.id == "aged") &
                         + qc.young >= .75*sum(cell.id == "young"),]

# choose variability statistic using diagnostics function.
     # choose the variability statistics that has the smallest correlation with the mean

diagnosticsVarPlotsTwoSample(logGV.age_young[1:5000,], groups=as.factor(c(rep(1,4),rep(2,4))))  

# find significant pathways using density and reactome:
grp<-c(rep(1, sum(cell.id == "aged")), rep(2, sum(cell.id == "young")))
set.seed(1)
resTwoSam<-pathVarTwoSamplesCont(logGV.age_young,pways.reactome.mouse,groups=as.factor(grp),varStat="sd") 
resTwoSam@tablePway[1:10]
sum(resTwoSam@tablePway[["APval"]]<0.05)
sigTwoSam<-sigPway(resTwoSam,0.05) # genes belonging to each sig pathway 

#show the density of expression variability for genes in the top 3 sig pathways for Group 1 and group 2
plotPway(resTwoSam,"PKMTs methylate histone lysines",sigTwoSam)
plotPway(resTwoSam,"HDMs demethylate histones",sigTwoSam)  
plotPway(resTwoSam,"SUMOylation of chromatin organization proteins", sigTwoSam) 

#what genes are the most different between the two groups for each pathway?
genes<-getGenes(resTwoSam,"PKMTs methylate histone lysines",c(0.2,0.6)) # what values should I put inside c()
setdiff(genes@genes1,genes@genes2) 

genes<-getGenes(resTwoSam,"HDMs demethylate histones",c(0.2,0.6))
setdiff(genes@genes1,genes@genes2) 

genes<-getGenes(resTwoSam,"SUMOylation of chromatin organization proteins",c(0.2,0.6))
setdiff(genes@genes1,genes@genes2) 

# find significant pathways using distribution counts and reactome (using the twoSampleDisc function to find pathways with significant difference in variability between the two groups using distribution counts)
resTwoSamDisc<-pathVarTwoSamplesDisc(logGV.age_young,pways.reactome.mouse,groups=as.factor(grp),test="exact",varStat="sd")
resTwoSamDisc@tablePway[1:20]
sum(resTwoSamDisc@tablePway[["APval"]]<0.01)
sigTwoSamDisc<-sigPway(resTwoSamDisc,0.01) 
plotPway(resTwoSamDisc,"PKMTs methylate histone lysines",sigTwoSamDisc)
plotPway(resTwoSamDisc,"HDMs demethylate histones",sigTwoSamDisc)
plotPway(resTwoSamDisc,"SUMOylation of chromatin organization proteins",sigTwoSamDisc)
plotAllTwoSampleDistributionCounts(logGV.age_young, resTwoSamDisc, perc=c(1/3,2/3), pvalue=0.05, NULL)

##########################################################################################################################

#PART 2: GV Stage Using Kegg Pathways:

#using distribution counts and Kegg
twoSamDisc_kegg<-pathVarTwoSamplesDisc(logGV.age_young,pways.kegg.mouse,groups=as.factor(grp),test="exact",varStat="sd")
twoSamDisc_kegg@tablePway[1:20]
sum(twoSamDisc_kegg@tablePway[["APval"]]<0.01)
sigTwoSamDisc_kegg<-sigPway(twoSamDisc_kegg,0.01) 
plotPway(twoSamDisc_kegg,"Spliceosome",sigTwoSamDisc_kegg)
plotPway(twoSamDisc_kegg,"DNA replication",sigTwoSamDisc_kegg)
plotPway(twoSamDisc_kegg,"Base excision repair",sigTwoSamDisc_kegg)
plotAllTwoSampleDistributionCounts(logGV.age_young, twoSamDisc_kegg, perc=c(1/3,2/3), pvalue=0.05, NULL)

###########################################################################################################################

#Part 3: MII stage using Kegg

# Identify our two groups: GV aged versus GV young) 
samp.id <- colnames(logMII)
cell.id <- character(length(samp.id))
cell.id[grep("aged", samp.id)] <- "aged"
cell.id[grep("young", samp.id)] <- "young"

#filter lowly expressed genes (keep the transcripts that have expression >1 in at least 75% of the replicates)
qc.aged <- apply(logMII[,cell.id == "aged"], 1, function(x,ct){ sum(x >= ct) }, ct=1)
qc.young <- apply(logMII[,cell.id == "young"], 1, function(x,ct){ sum(x >= ct) }, ct=1)
logMII.age_young <- logMII[qc.aged >= .75*sum(cell.id == "aged") &
                           + qc.young >= .75*sum(cell.id == "young"),]

#using distribution counts and Kegg 
twoSamDisc_keggMII<-pathVarTwoSamplesDisc(logMII.age_young,pways.kegg.mouse,groups=as.factor(grp),test="exact",varStat="sd")
View(twoSamDisc_keggMII@tablePway[1:20])
sum(twoSamDisc_keggMII@tablePway[["APval"]]<0.01)
sigTwoSamDisc_keggMII<-sigPway(twoSamDisc_keggMII,0.01) 
plotPway(twoSamDisc_keggMII,"Ascorbate and aldarate metabolism",sigTwoSamDisc_keggMII)
plotPway(twoSamDisc_keggMII,"Parkinson's disease",sigTwoSamDisc_keggMII)
plotPway(twoSamDisc_keggMII,"Oxidative phosphorylation",sigTwoSamDisc_keggMII)
plotAllTwoSampleDistributionCounts(logMII.age_young, twoSamDisc_keggMII, perc=c(1/3,2/3), pvalue=0.05, NULL)

#######################################################################################

#Part 4: GV Stage Using Hallmark:

hallmark_list <- gmt2list("h.all.v6.1.symbols.gmt")
pways.hallmark <- NULL
pways.hallmark$PATHID <- names(hallmark_list)
pways.hallmark$PATHNAME <- names(hallmark_list)
pways.hallmark$GENES <- hallmark_list
pways.hallmark$SIZE <- as.numeric(lapply(pways.hallmark$GENES, length))
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
} 
pways.hallmark$GENES <- lapply(pways.hallmark$GENES,capwords, strict=T)

#using distribution counts and hallmark
twoSamDisc_hallmark<-pathVarTwoSamplesDisc(logGV.age_young,pways.hallmark,groups=as.factor(grp),test="exact",varStat="sd")
View(twoSamDisc_hallmark@tablePway[1:20])
sum(twoSamDisc_hallmark@tablePway[["APval"]]<0.01)
sigTwoSamDisc_hallmark<-sigPway(twoSamDisc_hallmark,0.01) 
plotPway(twoSamDisc_hallmark,"HALLMARK_TNFA_SIGNALING_VIA_NFKB",sigTwoSamDisc_hallmark)
plotPway(twoSamDisc_hallmark,"HALLMARK_TGF_BETA_SIGNALING",sigTwoSamDisc_hallmark)
plotPway(twoSamDisc_hallmark,"HALLMARK_MITOTIC_SPINDLE",sigTwoSamDisc_hallmark)
plotAllTwoSampleDistributionCounts(logGV.age_young, twoSamDisc_hallmark, perc=c(1/3,2/3), pvalue=0.05, NULL)

#######################################################################################

#Part 5: MII Stage Using Hallmark:

twoSamDisc_hallmarkMII<-pathVarTwoSamplesDisc(logMII.age_young,pways.hallmark,groups=as.factor(grp),test="exact",varStat="sd")
View(twoSamDisc_hallmarkMII@tablePway[1:20])
sum(twoSamDisc_hallmarkMII@tablePway[["APval"]]<0.01)#returns none!  
sum(twoSamDisc_hallmarkMII@tablePway[["APval"]]<0.05)# only returns 1
sigTwoSamDisc_hallmarkMII<-sigPway(twoSamDisc_hallmarkMII,0.05) 
plotPway(twoSamDisc_hallmarkMII,"HALLMARK_UNFOLDED_PROTEIN_RESPONSE",sigTwoSamDisc_hallmarkMII)
plotAllTwoSampleDistributionCounts(logMII.age_young, twoSamDisc_hallmarkMII, perc=c(1/3,2/3), pvalue=0.05, NULL)




##2020.02.22
## Okamoto wheat mRNA-seq with udonko
## GOStat with the detected DEGs from Kallisto-Slueth

## Load libraries
library(AnnotationForge)
library(GSEABase)
library(GOstats)

## Load GO file
go <- read.table("IWGSC.GOslim.txt", header=T)
go <- go[,c(2,3,1)]
dim(go)
#[1] 666123      3
head(go)
#        GO_id Evidence            Gene_id
#1 GO:0003674      IEA TraesCS3A02G154800
#2 GO:0003824      IEA TraesCS3A02G154800
#3 GO:0016787      IEA TraesCS3A02G154800
#4 GO:0003674      IEA TraesCS2B02G516800
#5 GO:0003824      IEA TraesCS2B02G516800
#6 GO:0016787      IEA TraesCS2B02G516800


## Load Controls 
goFrame <- GOFrame(go, organism="wheat")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

########### Run GOstat
deg <- read.table("deg.input.ls.txt", header=T)
degs <- as.vector(deg[,1])
all <- as.vector(go$Gene_id)  ## a vector of all referenced list
length(unique(all))
#[1] 76029

p <- GSEAGOHyperGParams(
  name = "Paramaters",
  geneSetCollection = gsc,
  geneIds = degs,
  universeGeneIds = all,
  ontology = "MF",   ## "CC", "MF", "BP"
  pvalueCutoff = 0.05,
  conditional = FALSE,
  testDirection = "over"
)
result <- hyperGTest(p)
head(summary(result))
#GOMFID       Pvalue OddsRatio   ExpCount Count  Size
#1 GO:0003824 1.864996e-07  1.779756 154.803525   201 29656
#2 GO:0016787 2.213576e-06  1.905799  45.659105    77  8747
#3 GO:0140110 1.035536e-04  2.609848   8.863514    22  1698
#4 GO:0003700 1.035536e-04  2.609848   8.863514    22  1698
#5 GO:0030246 9.687753e-04  2.559318   6.493647    16  1244
#6 GO:0003677 3.950396e-03  1.573793  31.596498    47  6053
#Term
#1                        catalytic activity
#2                        hydrolase activity
#3          transcription regulator activity
#4 DNA-binding transcription factor activity
#5                      carbohydrate binding
#6                               DNA binding

write.table(summary(result), file = "1-dn1000.MF.txt",quote=FALSE,sep="\t",row.names=FALSE)


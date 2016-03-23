sed -i.bak 's/XGR/XGR/g' *.r
rm -rf *.r.bak

############################################################################
# Create a new R package called 'XGR'
############################################################################
cd ~/Sites/XGR
R

# base::dir.create("./XGR") # only for the very first time

# remove ./XGR/man/*.Rd
base::do.call(file.remove, list(list.files(path="./XGR/man", pattern="\\.Rd$", full.names=T)))
# remove ./XGR/R/*.r
base::do.call(file.remove, list(list.files(path="./XGR/R", pattern="\\.r$", full.names=T)))
# remove ./XGR/*.md
base::do.call(file.remove, list(list.files(path="./XGR", pattern="\\.md$", full.names=T)))

# remove and recreate the whole directory ./XGR/inst
base::unlink("./XGR/inst", recursive=T)
base::dir.create("./XGR/inst")
# remove and recreate the whole directory ./XGR/vignettes
base::unlink("./XGR/vignettes", recursive=T)
base::dir.create("./XGR/vignettes")

# create a skeleton for a new source package
utils::package.skeleton(name="XGR", path=".", code_files=c("xRDataLoader.r","xRdWrap.r","xFunArgs.r","xRd2HTML.r","xDAGanno.r","xDAGsim.r","xConverter.r","xEnricher.r","xEnricherGenes.r","xEnricherSNPs.r","xEnricherYours.r","xEnrichViewer.r","xSocialiser.r","xSocialiserGenes.r","xSocialiserSNPs.r","xCircos.r","xSubneterGenes.r","xSubneterSNPs.r","xVisNet.r","xPrioritiser.r","xPrioritiserGenes.r","xPrioritiserSNPs.r","xPrioritiserPathways.r","xPrioritiserManhattan.r"), force=T)

# do roxygenizing to document a package
devtools::document(pkg="XGR", roclets=c("collate", "namespace", "rd"))
# remove ./XGR/man/*-package.Rd
do.call(file.remove, list(list.files(path="./XGR/man", pattern="*-package.Rd", full.names=T)))

## wrap the lines with a long text for ./XGR/man/*.Rd files
source("./XGR/R/xRdWrap.r")
# use it when submitting it to CRAN
xRdWrap(path="./XGR/man", remove.dontrun=FALSE)

## create ./XGR/vignettes/*.html files
source("./XGR/R/xRd2HTML.r")
xRd2HTML(path.from="./XGR/man", path.to="./XGR/vignettes")

## extract the function title
#grep -E 'title' XGR/man/xRd2HTML.Rd | perl -lape '~s/\\title{//g;~s/}//g'

q('no')

############################################################################
# Build/Check/Install/Remove the R package called 'XGR'
############################################################################
cd ~/Sites/XGR

# Prepare the package to be built
## remove unnecessary files
find ./ -name '.DS_Store' -exec rm -rf {} \;
rm ./XGR/Read-and-delete-me
## copy necessary files
cp -r now.DESCRIPTION ./XGR/DESCRIPTION
cp -r now.README.md ./XGR/README.md
cp -r now.INSTALLATION.md ./XGR/inst/INSTALLATION.md
cp -r now.CITATIONS.md ./XGR/inst/CITATIONS.md
cp -r now.XGR.logo.png ./XGR/inst/XGR.logo.png
cp -r now.XGR.icon.png ./XGR/inst/XGR.icon.png
cp -r now.XGR_user.Rnw ./XGR/vignettes/XGR_user.Rnw
cp -r now.HTML.md ./XGR/vignettes/HTML.md

## copy other staffs to ./XGR/inst
cp -r XGR.sh ./XGR/inst/XGR.sh
#cp -r XGR.ppt ./XGR/inst/XGR.ppt
#cp -r XGR_talk.ppt ./XGR/inst/XGR_talk.ppt
#cp -r ~/Sites/XGR/XGR-site/xgr_app ./XGR/inst

# Build the package source
R CMD build --resave-data XGR
# Check the built package source
#R CMD check --timings --as-cran XGR_0.99.0.tar.gz
R CMD check --timings XGR_0.99.0.tar.gz
cat XGR.Rcheck/XGR-Ex.timings
rm -rf XGR.Rcheck

################
# manually build user and reference manuals
## user manual
R CMD Sweave --clean --pdf ./XGR/vignettes/XGR_user.Rnw
rm -rf *.tex
## reference manual
R CMD Rd2pdf --no-preview --force -o XGR_reference.pdf ./XGR
rm -rf .Rd2pdf*
## move both pdf files to ./XGR/vignettes
mv XGR_user.pdf XGR_reference.pdf ./XGR/vignettes
################

# Install from the package source
R CMD install XGR_0.99.0.tar.gz
# Remove the package source
#R CMD remove XGR

############################################################################
# Host the package source in github at https://github.com/hfang-bristol
############################################################################
# First, create/initialise the package in github
mkdir ~/Sites/SVN/github
cd ~/Sites/SVN/github

######
git config --global user.name "Hai Fang"
git config --global user.email hfang.shanghai@gmail.com
git init
# Git via Git for git updating
git clone git@github.com:git/git.git
cd git
make prefix=/usr/local all
sudo make prefix=/usr/local install
######

# in https://github.com/account, cat ~/.ssh/id_dsa.pub
# in https://github.com/hfang-bristol, create XGR
# 1. An R package to increase interpretability of genes and SNPs identified from GWAS and eQTLs
# 2. check the box 'Initialize this repository with a README' and choose license 'GNU General Public License v2.0'

# Then, clone into the local machine 
cd ~/Sites/SVN/github
git clone git@github.com:hfang-bristol/XGR.git

# Then, fill in the package source and update the package source in the github
cd ~/Sites/SVN/github/XGR
git pull origin master

############
# delete previous *.r *.Rd *.html files
find ./ -type f -name '*.r' -o -type f -name '*.Rd' -o -type f -name '*.html' -o -type f -name '*.ppt' | xargs rm -rf
############

############
# copy the whole package from ~/Sites/XGR/XGR to ~/Sites/SVN/github/XGR
cp -rf ~/Sites/XGR/XGR/* ./
#git add --all
############

git status | grep 'modified' | awk '{print $2}' | xargs git add
git status | grep 'deleted' | awk '{print $2}' | xargs git rm
find ./ -type f -name '*.r' -o -type f -name '*.Rd' -o -type f -name '*.html'  -o -type f -name '*.ppt' | xargs git add
git add --all # for recursively adding the entire directory 'xgr_app'
git commit -m 'Update to the version 0.99.0'
git push origin master


###################
# Tracking Branches (https://git-scm.com/book/en/v2/Git-Branching-Remote-Branches)
###################
## https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
# create a branch (local) and switch to this branch
git checkout -b hfang
# push this branch to github
git push origin hfang

# be sure being in the branch when committing something in this branch
git branch
# or: list out the local branches with more information
git branch -vv 

# delete a branch (remote)
git push origin --delete hfang
# delete a branch (local)
git branch -D hfang
##################

# being at the master branch: merge with hfang branch
git checkout master # git checkout --track origin/master
# git pull origin hfang
git fetch origin hfang
git merge origin hfang
git push origin master


############################################################################
# Obtain the package built-in RData
############################################################################

# in galahad
cd ~
git clone git@github.com:hfang-bristol/RDataCentre.git
cd ~/RDataCentre

a <- read.table('aaaa.txt',stringsAsFactors=F)
ind <- match(as.character(a[,1]), allSymbols)

mat_info <- allGenes[ind,]
priority_targets <- cbind(mat_info, pNode$priority, targets=vec_targets)
write.table(priority_targets, file="AS.genes_priority.targets.txt", sep="\t", row.names=FALSE)



############################################################################
############################################################################
cd /Users/hfang/Databases/grant_app/jobs/Oxford/Oxford_Wellcome_Trust_Centre_for_Human_Genetics/ULTRA-DD/TPN3rd_20160321/Data

R
library(XGR)
library(igraph)
library(dnet)
library(GenomicRanges)
library(ggbio)

RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/0.99.0"
EG <- xRDataLoader(RData=paste('org.Hs.eg', sep=''), RData.location=RData.location)
allGenes <- EG$gene_info[,1]
allSymbols <- as.vector(EG$gene_info$Symbol)
allDes <- as.vector(EG$gene_info$description)
allSynonyms <- as.vector(EG$gene_info$Synonyms)
allSynonyms <- gsub("\\|", ",", allSynonyms, perl=T)
allInfo <- cbind(allGenes, allSymbols, allDes, allSynonyms)

nomi <- read.delim(file="Targets_714.txt", stringsAsFactors=F)
ind <- match(nomi[,1], allGenes)
nomi_match <- allInfo[ind,]
df_nomi <- cbind(nomi_match, nomi)

## AS
AS <- read.delim(file="AS.genes_priority.txt", stringsAsFactors=F)
ind <- match(df_nomi$allSymbols, AS$name)
df_AS <- AS[ind,]

## Spondyloarthritis
SP <- read.delim(file="Spondyloarthritis.genes_priority.txt", stringsAsFactors=F)
ind <- match(df_nomi$allSymbols, SP$name)
df_SP <- SP[ind,]

## SLE
SLE <- read.delim(file="SLE.genes_priority.txt", stringsAsFactors=F)
ind <- match(df_nomi$allSymbols, SLE$name)
df_SLE <- SLE[ind,]

## combine
df <- cbind(AS=df_AS$rank, SP=df_SP$rank, SLE=df_SLE$rank, df_nomi)
df$AS[df$AS>500] <- NA 
df$SP[df$SP>500] <- NA
df$SLE[df$SLE>500] <- NA 
write.table(df, file="Nominated_Targets.txt", sep="\t", row.names=FALSE)



############################################################################
############################################################################




############################################################################
############################################################################
cp ~/Desktop/AS.txt ~/Sites/XGR/PI-site/app/examples/AS.txt
cp ~/Desktop/SLE.txt ~/Sites/XGR/PI-site/app/examples/SLE.txt
cp ~/Desktop/Spondyloarthritis.txt ~/Sites/XGR/PI-site/app/examples/Spondyloarthritis.txt
############################################################################
############################################################################

##########################
## for SLE
# http://www.ncbi.nlm.nih.gov/pubmed/26502338
##########################
# in galahad
mysql -uhfang -p5714fh -e "use ultraDDR; SELECT snp_id_current,pvalue FROM GWAS WHERE EF_name like '%Systemic Lupus Erythematosus%' and pvalue is not null;" > GWAS_SLE.txt
# in mac
cd ~/
scp galahad.well.ox.ac.uk:./GWAS_SLE.txt ./

R
# Load the library
library(XGR)
library(igraph)
library(dnet)
library(GenomicRanges)
library(ggbio)

RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/0.99.0"
## load ImmunoBase
ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
## get lead SNPs reported in AS GWAS and their significance info (p-values)
gr <- ImmunoBase$SLE$variant
seeds.snps <- as.matrix(mcols(gr)[,c(1,3)])
write.table(seeds.snps, file="ImmunoBase_SLE.txt", sep="\t", row.names=FALSE, quote=F)
q('no')

cat GWAS_SLE.txt ImmunoBase_SLE.txt | grep -E '^rs' | sort | uniq > All_SLE.txt

##########################
## for Spondyloarthritis
# http://www.ncbi.nlm.nih.gov/pubmed/25651891
# http://www.ncbi.nlm.nih.gov/pubmed/26974007
##########################
# in galahad
mysql -uhfang -p5714fh -e "use ultraDDR; SELECT snp_id_current,pvalue FROM GWAS WHERE EF_name like '%psoriatic arthritis%' and pvalue is not null;;" > GWAS_PSO.txt
# in mac
cd ~/
scp galahad.well.ox.ac.uk:./GWAS_PSO.txt ./

R
# Load the library
library(XGR)
library(igraph)
library(dnet)
library(GenomicRanges)
library(ggbio)

RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/0.99.0"
## load ImmunoBase
ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
## get lead SNPs reported in Psoriasis GWAS and their significance info (p-values)
gr <- ImmunoBase$PSO$variant
seeds.snps <- as.matrix(mcols(gr)[,c(1,3)])
write.table(seeds.snps, file="ImmunoBase_PSO.txt", sep="\t", row.names=FALSE, quote=F)

cat GWAS_PSO.txt ImmunoBase_PSO.txt | grep -E '^rs' | sort | uniq > All_PSO.txt








# in galahad
mysql -uhfang -p5714fh -e "use ultraDDR; SELECT snp_id_current,pvalue FROM GWAS WHERE EF_name like '%ankylosing spondylitis%' and pvalue is not null;" > GWAS_AS.txt

mysql -uhfang -p5714fh -e "use ultraDDR; SELECT snp_id_current,pvalue FROM GWAS WHERE EF_name like '%Systemic Lupus Erythematosus%' and pvalue is not null;" > GWAS_SLE.txt

mysql -uhfang -e "use ultraDDR; SELECT snps,pvalue FROM gwas_ebi WHERE mapped_trait like '%Systemic Lupus Erythematosus%';" > GWAS_SLE.txt

# in mac
cd ~/
scp galahad.well.ox.ac.uk:./GWAS_AS.txt ./
scp galahad.well.ox.ac.uk:./GWAS_SLE.txt ./

cat GWAS_AS.txt ImmunoBase_AS.txt | grep -E '^rs' | sort | uniq > All_AS.txt
cat GWAS_SLE.txt ImmunoBase_SLE.txt | grep -E '^rs' | sort | uniq > All_SLE.txt

cp ~/All_AS.txt ~/Sites/XGR/PI-site/app/examples/AS.txt
cp ~/All_SLE.txt ~/Sites/XGR/PI-site/app/examples/SLE.txt

# Load the library
library(XGR)
library(igraph)
library(dnet)
library(GenomicRanges)
library(ggbio)

RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/0.99.0"
ImmunoBase_LD <- xRDataLoader(RData.customised='ImmunoBase_LD', RData.location=RData.location)

# a) provide the seed SNPs with the weight info
## load ImmunoBase
ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
## get lead SNPs reported in AS GWAS and their significance info (p-values)
gr <- ImmunoBase$AS$variant
seeds.snps <- as.matrix(mcols(gr)[,c(1,3)])
write.table(seeds.snps, file="ImmunoBase_AS.txt", sep="\t", row.names=FALSE, quote=F)




##########################################################################################
##########################################################################################
##########################################################################################

EG <- xRDataLoader(RData=paste('org.Hs.eg', sep=''), RData.location=RData.location)
allSymbols <- as.vector(EG$gene_info$Symbol)
allGenes <- EG$gene_info[,c(3,2,1)]

############# Ankylosing Spondylitis (AS) #############
# perform priority analysis
pNode <- xPrioritiserSNPs(data="/Users/hfang/Sites/XGR/PI-site/app/examples/AS.txt", include.LD="EUR", include.eQTL=c("JKscience_TS2B","JKscience_TS3A"), network="STRING_high", significance.threshold=5e-5, distance.max=200000, restart=0.75, RData.location=RData.location)
write.table(pNode$priority, file="AS.genes_priority.txt", sep="\t", row.names=FALSE)

mp <- xPrioritiserManhattan(pNode, highlight.top=25, RData.location=RData.location, cex=0.4, highlight.col="deepskyblue", highlight.label.size=1.8, highlight.label.offset=0.01)
pdf('AS.genes_priority.pdf', height=5, width=10, compress=TRUE)
print(mp)
dev.off()

# derive pathway-level priority
eTerm <- xPrioritiserPathways(pNode=pNode, priority.top=100, ontology="MsigdbC2CPall", RData.location=RData.location)
res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
res_f <- data.frame(term=rownames(res), res)
write.table(res_f, file="AS.pathways_priority.txt", sep="\t", row.names=FALSE)

# append targets
targets <- read.table(file="~/Desktop/targets.txt", header=F, stringsAsFactors=F)
ind <- match(rownames(pNode$priority), targets[,1])
vec_targets <- rep('No', length(ind))
vec_targets[!is.na(ind)] <- 'Yes'
ind <- match(rownames(pNode$priority), allSymbols)
mat_info <- allGenes[ind,]
priority_targets <- cbind(mat_info, pNode$priority, targets=vec_targets)
write.table(priority_targets, file="AS.genes_priority.targets.txt", sep="\t", row.names=FALSE)



############# Spondyloarthritis (including AS and Psoriatic Arthritis) #############
# perform priority analysis
pNode <- xPrioritiserSNPs(data="/Users/hfang/Sites/XGR/PI-site/app/examples/Spondyloarthritis.txt", include.LD="EUR", include.eQTL=c("JKscience_TS2B","JKscience_TS3A"), network="STRING_high", significance.threshold=5e-5, distance.max=200000, restart=0.75, RData.location=RData.location)
write.table(pNode$priority, file="Spondyloarthritis.genes_priority.txt", sep="\t", row.names=FALSE)

mp <- xPrioritiserManhattan(pNode, highlight.top=25, RData.location=RData.location, cex=0.4, highlight.col="deepskyblue", highlight.label.size=1.8, highlight.label.offset=0.01)
pdf('Spondyloarthritis.genes_priority.pdf', height=5, width=10, compress=TRUE)
print(mp)
dev.off()

# derive pathway-level priority
eTerm <- xPrioritiserPathways(pNode=pNode, priority.top=100, ontology="MsigdbC2CPall", RData.location=RData.location)
res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
res_f <- data.frame(term=rownames(res), res)
write.table(res_f, file="Spondyloarthritis.pathways_priority.txt", sep="\t", row.names=FALSE)

# append targets
targets <- read.table(file="~/Desktop/targets.txt", header=F, stringsAsFactors=F)
ind <- match(rownames(pNode$priority), targets[,1])
vec_targets <- rep('No', length(ind))
vec_targets[!is.na(ind)] <- 'Yes'
ind <- match(rownames(pNode$priority), allSymbols)
mat_info <- allGenes[ind,]
priority_targets <- cbind(mat_info, pNode$priority, targets=vec_targets)
write.table(priority_targets, file="Spondyloarthritis.genes_priority.targets.txt", sep="\t", row.names=FALSE)



############# Systemic Lupus Erythematosus (SLE) #############
# perform priority analysis
pNode <- xPrioritiserSNPs(data="/Users/hfang/Sites/XGR/PI-site/app/examples/SLE.txt", include.LD="EUR", include.eQTL=c("JKscience_TS2B","JKscience_TS3A"), network="STRING_high", significance.threshold=5e-5, distance.max=200000, restart=0.75, RData.location=RData.location)
write.table(pNode$priority, file="SLE.genes_priority.txt", sep="\t", row.names=FALSE)

mp <- xPrioritiserManhattan(pNode, highlight.top=25, RData.location=RData.location, cex=0.4, highlight.col="deepskyblue", highlight.label.size=1.8, highlight.label.offset=0.01)
pdf('SLE.genes_priority.pdf', height=5, width=10, compress=TRUE)
print(mp)
dev.off()

# derive pathway-level priority
eTerm <- xPrioritiserPathways(pNode=pNode, priority.top=100, ontology="MsigdbC2CPall", RData.location=RData.location)
res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
res_f <- data.frame(term=rownames(res), res)
write.table(res_f, file="SLE.pathways_priority.txt", sep="\t", row.names=FALSE)

# append targets
targets <- read.table(file="~/Desktop/targets.txt", header=F, stringsAsFactors=F)
ind <- match(rownames(pNode$priority), targets[,1])
vec_targets <- rep('No', length(ind))
vec_targets[!is.na(ind)] <- 'Yes'
ind <- match(rownames(pNode$priority), allSymbols)
mat_info <- allGenes[ind,]
priority_targets <- cbind(mat_info, pNode$priority, targets=vec_targets)
write.table(priority_targets, file="SLE.genes_priority.targets.txt", sep="\t", row.names=FALSE)

##########################################################################################
##########################################################################################
##########################################################################################













a <- read.table('All_AS.txt', header=F)

# cis-eQTL
cis <- xRDataLoader(RData.customised='JKscience_TS2B', RData.location=RData.location)
minFDR <- apply(cis[,c(9:12)], 1, min, na.rm=T)
cis_df <- data.frame(SNP=cis[,1], Gene=cis[,4], FDR=minFDR, stringsAsFactors=F)
# trans-eQTL
trans <- xRDataLoader(RData.customised='JKscience_TS3A', RData.location=RData.location)
minFDR <- apply(trans[,c(9:12)], 1, min, na.rm=T)
trans_df <- data.frame(SNP=trans[,1], Gene=trans[,4], FDR=minFDR, stringsAsFactors=F)
# both eQTL
e_df <- rbind(cis_df, trans_df)
uid <- paste(e_df[,1], e_df[,2], sep='_')
df <- cbind(uid, e_df)
res_list <- split(x=df$FDR, f=df$uid)
res <- lapply(res_list, function(x){
	min(x)
})
vec <- unlist(res)
res_df <- do.call(rbind, strsplit(names(vec),'_'))
both_df <- data.frame(SNP=res_df[,1], Gene=res_df[,2], FDR=vec, row.names=NULL, stringsAsFactors=F)



Lead_Sig <- read.table('All_AS.txt', header=F, stringsAsFactors=F)
colnames(Lead_Sig) <- c('SNP','Sig')
leads <- Lead_Sig[,1]
sigs <- Lead_Sig[,2]

include.LD <- c('EUR','AMR','AFR','EAS','EUR')
LD.r2 <- 0.8
if(length(include.LD) > 0){
	GWAS_LD <- xRDataLoader(RData.customised='GWAS_LD', RData.location=RData.location, verbose=verbose)
	res_list <- lapply(include.LD, function(x){
		data_ld <- ''
		eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
		ind <- match(rownames(data_ld), leads)
		ind_lead <- which(!is.na(ind))
		ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
		
		sLL <- data_ld[ind_lead, ind_ld]
		summ <- summary(sLL)
		res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
	})
	## get data frame (Lead LD R2)
	LLR <- do.call(rbind, res_list)
	
	## get data frame (LD Sig)
	ld_list <- split(x=LLR[,-2], f=LLR[,2])
	res_list <- lapply(ld_list, function(x){
		ind <- match(x$Lead, leads)
		## power transformation of p-valuesn by R2, then keep the min
		min(sigs[ind] ^ x$R2)
	})
	vec <- unlist(res_list)
	LD_Sig <- data.frame(SNP=names(vec), Sig=vec, row.names=NULL, stringsAsFactors=F)
	
	## merge Lead and LD
	df <- rbind(Lead_Sig, as.matrix(LD_Sig))
	res_list <- split(x=df$Sig, f=df$SNP)
	res <- lapply(res_list, function(x){
		min(x)
	})
	vec <- unlist(res)
	SNP_Sig <- data.frame(SNP=names(vec), FDR=vec, row.names=NULL, stringsAsFactors=F)
	
}else{
	SNP_Sig <- Lead_Sig
}

pval <- as.numeric(SNP_Sig[,2])
names(pval) <- SNP_Sig[,1]
###########################
# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
significance.threshold=0.01
scores <- log10(2) * dFDRscore(pval, fdr.threshold=significance.threshold, scatter=F)
scores[scores<0] <- 0
seeds.snps <- scores



ind <- match(a[,1], res_df[,1])


res <- pRes$matrix

g <- pRes$g
priority <- pRes$priority
nodes_query <- names(pRes$seeds)

tmp <- priority[nodes_query]
tmp <- tmp[tmp > quantile(tmp,0.8)]

subg <- dNetInduce(g, nodes_query=names(tmp), knn=0)
xVisNet(g=subg, pattern=tmp, colormap="yellow-red", vertex.shape="sphere", zlim=c(min(tmp),max(tmp)))






# SNP-based enrichment analysis using GWAS Catalog traits (mapped to EF)
# a) provide the input SNPs of interest (eg 'EFO:0002690' for 'systemic lupus erythematosus')
## load GWAS SNPs annotated by EF (an object of class "dgCMatrix" storing a spare matrix)
anno <- xRDataLoader(RData='GWAS2EF', RData.location=RData.location)
ind <- which(colnames(anno)=='EFO:0002690')
data <- rownames(anno)[anno[,ind]==1]
data

eTerm <- xEnricherSNPs(data=data, ontology="EF", include.LD="EUR", path.mode=c("all_paths"), RData.location=RData.location)
xEnrichViewer(eTerm)

eT <- xEnricherSNPs(data=data, ontology="EF", include.LD="EUR", path.mode=c("all_paths"), RData.location=RData.location)
xEnrichViewer(eT)
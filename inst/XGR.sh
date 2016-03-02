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
utils::package.skeleton(name="XGR", path=".", code_files=c("xRDataLoader.r","xRdWrap.r","xFunArgs.r","xRd2HTML.r","xDAGanno.r","xDAGsim.r","xConverter.r","xEnricher.r","xEnricherGenes.r","xEnricherSNPs.r","xEnricherYours.r","xEnrichViewer.r","xSocialiser.r","xSocialiserGenes.r","xSocialiserSNPs.r","xCircos.r","xSubneter.r","xVisNet.r"), force=T)

# do roxygenizing to document a package
devtools::document(pkg="XGR", clean=FALSE, roclets=c("collate", "namespace", "rd"), reload=TRUE)
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
cp -r now.XGR.icon.png ./XGR/inst/XGR.logo.png
cp -r now.XGR_user.Rnw ./XGR/vignettes/XGR_user.Rnw
cp -r now.HTML.md ./XGR/vignettes/HTML.md

## copy other staffs to ./XGR/inst
cp -r XGR.sh ./XGR/inst/XGR.sh
cp -r XGR.ppt ./XGR/inst/XGR.ppt
cp -r XGR_talk.ppt ./XGR/inst/XGR_talk.ppt
cp -r ~/Sites/XGR-site/xor_app ./XGR/inst

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
# 1. An R package for deciphering genes and SNPs identified from GWAS and eQTLs
# 2. check the box 'Initialize this repository with a README' and choose license 'GNU General Public License v2.0'

# Then, clone into the local machine 
cd ~/Sites/SVN/github
git clone -b hfang git@github.com:hfang-bristol/XGR.git

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
git add --all # for recursively adding the entire directory 'xor_app'
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



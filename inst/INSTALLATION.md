# INSTALLATION 

## 1. R requirement

R (http://www.r-project.org) is a language and environment for statistical computing and graphics. We assume R (version 3.1.0 or higher) has been installed in your local machine. The latest version can be installed following instructions below under different platforms (Windows, Mac, and Linux).

* Quick link for `Windows`: [Download R for Windows](http://cran.r-project.org/bin/windows/base/R-3.2.4-win.exe).
* Quick link for `Mac`: [Download R for Mac OS X 10.6 (Snow Leopard or higher)](http://cran.r-project.org/bin/macosx).

* Below are `shell command lines in Terminal` (for `Linux`):

Assume you have a `ROOT (sudo)` privilege:
>
```ruby
sudo su
# here enter your password
wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.2.4.tar.gz
tar xvfz R-3.2.4.tar.gz
cd R-3.2.4
./configure
make
make check
make install
R # start R
```

Assume you do not have a `ROOT` privilege and want R installation under your home directory (`$HOME`):
>
```ruby
wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.2.4.tar.gz
tar xvfz R-3.2.4.tar.gz
cd R-3.2.4
./configure --prefix=$HOME/R-3.2.4
make
make check
make install
$HOME/R-3.2.4/bin/R # start R
```

## 2. Installation of the package

Notes: below are `R command lines (NOT shell command lines in Terminal)`.

First, install the package `devtools` (to help install packages directly from github) and other suggested packages:
>
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite(c("devtools","dnet","RCircos","ggbio"))
```

Second, install the package `XGR` from [GitHub](https://github.com/hfang-bristol/XGR):
>
```{r}
library(devtools)
install_github(c("hfang-bristol/dnet","hfang-bristol/XGR"), ref="master", dependencies=T)
```
# Running STRAF on a local computer {#local-computer}

## Intro

STRAF is a web application
It is based on the Shiny web application framework, R programming language.

It runs on a server and can be accessed through the browser.

But it can also be ran on a local computer, provided that you have relevant
dependencies installed on your computer.

This can be useful if you have any reason not to upload data on an external server.


## Prerequisites

### Installing R

If your system is running under Windows, you can follow the following link:

[https://cran.r-project.org/bin/windows/base/](http://cran.r-project.org/bin/windows/base/)

You can dowmnload the installation software from here. you can then follow the
instruction to install R on your system, keeping default options.

[https://cran.r-project.org/bin/macosx/](http://cran.r-project.org/bin/macosx/)

If you are using a Linux system, vous devriez pouvoir trouver R via votre gestionnaire de paquets, cela pouvant dépendre d’une distribution de Linux à une autre
As it is common, prior to installing R, let us update the system package index and upgrade all our installed packages using the following two commands:

sudo apt update

sudo apt -y upgrade

After that, all that you have to do is run the following in the command line to install base R.

sudo apt -y install r-base


Now, you should be able to run R-based software on your computer.

## Install RStudio

RStudio is a free an open-source user-friendly environment built on top of R

You can download RStudio Desktop from https://www.rstudio.com/products/rstudio/download

RStudio: RStudio is a free and open source integrated development environment (IDE) for R. While you can write and use Shiny apps with any R environment (including R GUI and ESS), RStudio has some nice features specifically for authoring, debugging, and deploying Shiny apps. We recommend giving it a try, but it’s not required to be successful with Shiny or with this book. 
    R packages: This book uses a bunch of R packages. You can install them all at once by running:

    install.packages(c(
      "gapminder", "ggforce", "gh", "globals", "openintro", "profvis", 
      "RSQLite", "shiny", "shinycssloaders", "shinyFeedback", 
      "shinythemes", "testthat", "thematic", "tidyverse", "vroom", 
      "waiter", "xml2", "zeallot" 
    ))


## Install the shiny and straf packages

## Launch the application



_Et voilà!_

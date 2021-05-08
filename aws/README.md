## Install STRAF on AWS

### Install shiny-server

```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install -y r-base r-base-dev libudunits2-dev libcurl4-openssl-dev libgdal-dev
sudo R -e "install.packages('shiny', INSTALL_opts = '--no-lock')""
sudo apt-get install -y gdebi-core
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.15.953-amd64.deb
sudo gdebi -y shiny-server-1.5.15.953-amd64.deb

sudo rm /srv/shiny-server/index.html
sudo rm -R /srv/shiny-server/sample-apps
```

### Install straf

```
echo $'library(straf)\ndir <- system.file("application", package = "straf")\nsetwd(dir)\nshiny::shinyAppDir(".")' > /srv/shiny-server/app.R

sudo cp /srv/shiny-server/straf/aws/shiny-server.conf /etc/shiny-server/shiny-server.conf
sudo systemctl reload shiny-server

R -e "install.packages(c('remotes'))"
R -e "remotes::install_github('agouy/straf')"

sudo -su shiny
R -e "install.packages(c('remotes', 'dbplyr', 'RSQLite'))"
R -e "remotes::install_github('agouy/straf')"
exit

sudo systemctl reload shiny-server
```

### Set up reverse proxy

```
sudo apt install -y nginx
sudo cp straf/aws/shiny.conf /etc/nginx/sites-available/shiny.conf

sudo nginx -t

cd /etc/nginx/sites-enabled
sudo ln -s /etc/nginx/sites-available/shiny.conf .
```

### certificates


## install packages

```
sudo add-apt-repository ppa:marutter/c2d4u
sudo apt-get update
apt install --no-install-recommends r-cran-rcppeigen
```

## Update STRAF on the server

One simple needs to update the R package from GitHub.

```
sudo R -e "remotes::install_github('agouy/straf')"
```

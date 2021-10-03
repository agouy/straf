## Install STRAF on AWS

### Install shiny-server

```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install -y r-base r-base-dev libudunits2-dev libcurl4-openssl-dev libgdal-dev
sudo R -e "install.packages('shiny', INSTALL_opts = '--no-lock')"
sudo apt-get install -y gdebi-core

wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.15.953-amd64.deb

sudo gdebi -n shiny-server-1.5.15.953-amd64.deb

sudo rm /srv/shiny-server/index.html
sudo rm -R /srv/shiny-server/sample-apps
```

### Install straf

Don't forget to add a 3838 rule.

```
sudo mkdir /srv/shiny-server/app
echo $'library(straf)\ndir <- system.file("application", package = "straf")\nsetwd(dir)\nshiny::shinyAppDir(".")' > /srv/shiny-server/app/app.R

git clone https://github.com/agouy/straf
sudo cp ./straf/app.R /srv/shiny-server/app/app.R
sudo cp ./straf/aws/shiny-server.conf /etc/shiny-server/shiny-server.conf
sudo cp ./straf/landing/index.html /srv/shiny-server/index.html

sudo systemctl reload shiny-server

sudo -su shiny
R -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
R -e "install.packages(c('remotes', 'dbplyr', 'RSQLite'), lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cran.rstudio.com/')"


R -e "remotes::install_github('agouy/straf')"
exit

sudo systemctl reload shiny-server
```

### Set up reverse proxy

```
sudo apt install -y nginx
sudo cp ./straf/aws/shiny.conf /etc/nginx/sites-available/shiny.conf

sudo nginx -t

cd /etc/nginx/sites-enabled
sudo ln -s /etc/nginx/sites-available/shiny.conf .
```

## Install certificates



## Update STRAF on the server

One simply needs to update the R package from GitHub.

```
sudo -su shiny
sudo R -e "remotes::install_github('agouy/straf')"
sudo systemctl reload shiny-server
```

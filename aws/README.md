## Install STRAF on AWS

```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install -y r-base r-base-dev libudunits2-dev libcurl4-openssl-dev libgdal-dev
sudo R -e "install.packages('shiny', INSTALL_opts = '--no-lock')""
sudo apt-get install -y gdebi-core
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.15.953-amd64.deb
sudo gdebi -y shiny-server-1.5.15.953-amd64.deb

cd /srv/shiny-server
sudo git clone https://github.com/agouy/straf
sudo cp /srv/shiny-server/straf/* /srv/shiny-server/ -r

cd /srv/shiny-server
sudo ln -s ~/straf .
sudo rm index.html
sudo rm -R sample-apps

sudo cp /srv/shiny-server/straf/aws/shiny-server.conf /etc/shiny-server/shiny-server.conf
sudo systemctl reload shiny-server

sudo -su shiny
R -e "install.packages(c('adegenet', 'ade4', 'pegas', 'hierfstat', 'DT', 'car', 'shinythemes', 'colourpicker', 'plotly', 'ggrepel'))"
R -e "install.packages(c('dbplyr', 'RSQLite'))"
exit

sudo systemctl reload shiny-server
sudo apt install -y nginx
sudo cp straf/aws/shiny.conf /etc/nginx/sites-available/shiny.conf

sudo nginx -t

cd /etc/nginx/sites-enabled
sudo ln -s /etc/nginx/sites-available/shiny.conf .
```

## Update

```
cd /srv/shiny-server/straf
sudo git pull
sudo cp /srv/shiny-server/straf/* /srv/shiny-server/ -r
```

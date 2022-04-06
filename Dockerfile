FROM rocker/shiny-verse:4.1.2

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libgdal-dev \
    libglpk-dev \
    libmysqlclient-dev

# RUN apt-get update && \
#    apt-get upgrade -y && \
#    apt-get clean

RUN Rscript -e "install.packages(c('remotes'))"
RUN Rscript -e "remotes::install_github('agouy/straf', ref = '2.0.8')"

RUN Rscript -e "install.packages(c('markdown'))"

COPY /app.R ./app.R

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('./', host = '0.0.0.0', port = 3838)"]
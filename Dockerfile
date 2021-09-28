FROM rocker/shiny:latest

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

COPY /app.R ./app.R

RUN Rscript -e "install.packages(c('remotes'))"
RUN Rscript -e "remotes::install_github('agouy/straf')"

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('./', host = '0.0.0.0', port = 3838)"]
FROM rocker/r-ver:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN apt-get install -y curl

RUN apt-get install -y libcurl4-openssl-dev

RUN apt-get install -y libssl-dev

RUN apt-get install -y libxml2-dev 

RUN apt-get install -y libgeos-dev

# Install R packages
RUN R -e 'install.packages("data.table", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("dplyr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("DT", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("Matrix", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-01-20")'
RUN R -e 'install.packages("Seurat", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2022-05-20")'
RUN R -e 'install.packages("shiny", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinydashboard", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinyjs", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinythemes", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("stringi", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("stringr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("vroom", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'

RUN R -e 'install.packages("BiocManager")'

RUN R -e 'BiocManager::install("GEOquery")'

# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Expose the port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

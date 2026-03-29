FROM rocker/r-ver:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

ENV RSPM=https://packagemanager.rstudio.com/cran/__linux__/focal/2022-05-20

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN apt-get install -y curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgeos-dev \
    libglpk40 \
    libglpk-dev

# CRAN packages (consistent snapshot)
RUN R -e 'install.packages(c( \
"cowplot","data.table","dplyr","DT","ff","ggfun","ggplot2","ggplotify","ggpubr", \
"ggrepel","ggridges","gridExtra","magrittr","Matrix","matrixStats","msigdbr", \
"patchwork","plotly","reshape2","Seurat","shiny","shinycssloaders", \
"shinydashboard","shinyjs","shinythemes","stringi","stringr","tidytree","gplots","remotes" \
), repos=Sys.getenv("RSPM"))'

# Bioconductor (locked to R 4.1)
RUN R -e 'install.packages("BiocManager", repos=Sys.getenv("RSPM"))'
RUN R -e 'BiocManager::install(version="3.14", ask=FALSE)'

RUN R -e 'BiocManager::install(c("ggtree","scuttle","BiocSingular","fgsea"), ask=FALSE)'

# GitHub package without pulling incompatible deps
RUN R -e 'remotes::install_github("XiaoLuo-boy/ggheatmap", dependencies=FALSE)'

# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Expose the port
EXPOSE 3838

# Run app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

FROM rocker/r-ver:4.5.2

MAINTAINER Marta Bica <mbica.compbio@gmail.com>

# System dependencies
RUN apt-get update && apt-get install -y \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgeos-dev \
    libglpk40 \
    libglpk-dev \
    libicu-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Force R to build from source (avoids binary mismatch issues)
RUN R -e "options(pkgType='source')"

# CRAN packages (one command, proper quoting)
RUN R -e "install.packages(c('cowplot','data.table','dplyr','DT','ff','ggfun','ggplot2','ggplotify','ggpubr','ggrepel','ggridges','gridExtra','magrittr','Matrix','matrixStats','msigdbr','patchwork','plotly','reshape2','Seurat','shiny','shinycssloaders','shinydashboard','shinyjs','shinythemes','stringr','tidytree','gplots'))"

# Bioconductor packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('scuttle','BiocSingular','fgsea','markeR'))"

# GitHub packages
RUN R -e "install.packages('devtools')"
RUN R -e "remotes::install_github('YuLab-SMU/ggtree')"
RUN R -e "devtools::install_github('XiaoLuo-boy/ggheatmap')"

# App setup
WORKDIR /home/app
COPY . /home/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

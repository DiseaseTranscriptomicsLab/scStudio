FROM rocker/r-ver:4.5.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN apt-get install -y curl

RUN apt-get install -y libcurl4-openssl-dev

RUN apt-get install -y libssl-dev

RUN apt-get install -y libxml2-dev 

RUN apt-get install -y libgeos-dev

RUN apt-get install -y libglpk40 libglpk-dev

 RUN apt-get install -y libicu-dev

RUN R -e "install.packages('curl', type='source')"

# Install R packages
RUN R -e 'install.packages("cowplot", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-29")'
RUN R -e 'install.packages("data.table", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("dplyr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("DT", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ff", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggfun", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggplot2", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggplotify", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggpubr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggrepel", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("ggridges", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("gridExtra", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("magrittr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("Matrix", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("matrixStats", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("msigdbr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("patchwork", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("plotly", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("reshape2", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("Seurat", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("shiny", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("shinycssloaders", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("shinydashboard", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("shinyjs", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("shinythemes", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("stringr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'
RUN R -e 'install.packages("tidytree", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'

RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("scuttle")'
RUN R -e 'BiocManager::install("BiocSingular")'
RUN R -e 'BiocManager::install("fgsea")'
RUN R -e 'BiocManager::install("markeR")'

RUN R -e 'install.packages("devtools")'
RUN R -e 'remotes::install_github("YuLab-SMU/ggtree")'
RUN R -e 'devtools::install_github("XiaoLuo-boy/ggheatmap")' 
RUN R -e 'install.packages("gplots", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2026-01-24")'

# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Expose the port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

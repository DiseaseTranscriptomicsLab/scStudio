FROM r-base:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN apt-get install -y curl

RUN apt-get install -y libcurl4-openssl-dev

RUN apt-get install -y libssl-dev

RUN apt-get install -y libxml2-dev 

# Install R packages
RUN R -e 'install.packages("shiny", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinycssloaders", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("ggplot2", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("markdown", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-03-20")'
RUN R -e 'install.packages("plotly", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("DT", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("shinydashboard", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinyjs", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("scico", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("matrixStats", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("Seurat", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2022-05-20")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'install.packages("stringr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("Matrix", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-01-20")'
RUN R -e 'install.packages("stringi", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("rPanglaoDB", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("readxl", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("testthat", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("data.table", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("ggrepel", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("shinythemes", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'

RUN R -e 'install.packages("msigdbr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("dplyr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("vroom", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("reshape2", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("gridExtra", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("umap", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("devtools", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("ff", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("colourpicker", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'
RUN R -e 'install.packages("remotes", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'

RUN R -e 'BiocManager::install("GEOquery")'

RUN R -e 'BiocManager::install("scater")'

RUN R -e 'BiocManager::install("scDblFinder")'

RUN R -e 'BiocManager::install("scran")'

RUN R -e 'BiocManager::install("scuttle")'

RUN R -e 'BiocManager::install("limma")'

RUN R -e 'BiocManager::install("BiocSingular")'

RUN R -e 'BiocManager::install("PCAtools")'

RUN R -e 'BiocManager::install("fgsea")'

RUN R -e 'remotes::install_github("daqana/dqshiny")'

RUN R -e 'remotes::install_github("ChristophH/gplots")'

RUN apt-get install -y libfontconfig1-dev

RUN R -e 'install.packages("clustree", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-04-20")'

RUN apt-get install -y cmake

RUN R -e 'install.packages("nloptr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-01-20")'

RUN R -e 'install.packages("ggpubr", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-01-20")'

RUN R -e 'BiocManager::install("MAST")'


RUN R -e "install.packages(c('withr'), repos='https://cloud.r-project.org/')"
RUN R -e "withr::with_makevars(c(PKG_CPPFLAGS='-DHTTP_MAX_HEADER_SIZE=0x7fffffff'),  {install.packages(c('shiny'), repos='https://cloud.r-project.org/')}, assignment = '+=')"
          
        
# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Expose the port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]
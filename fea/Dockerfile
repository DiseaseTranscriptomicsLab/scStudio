FROM rocker/tidyverse:4.5.2

LABEL maintainer="Marta Bica <mbica.compbio@gmail.com>"

# ---- Extra system dependencies ----
RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libgeos-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---- Force source installs ----
RUN R -e "options(pkgType='source')"

# ---- Install CRAN packages ----
RUN R -e "install.packages(c('cowplot','data.table','DT','ff','ggfun','ggplotify','ggpubr','ggrepel','ggridges','gridExtra','magrittr','Matrix','matrixStats','msigdbr','patchwork','plotly','reshape2','Seurat','shiny','shinycssloaders','shinydashboard','shinyjs','shinythemes','stringr','tidytree','gplots','BiocManager','devtools'))"

# ---- Install Bioconductor packages ----
RUN R -e "BiocManager::install(c('scuttle','BiocSingular','fgsea','markeR'), ask=FALSE, update=FALSE)"

# ---- Install GitHub packages ----
RUN R -e "remotes::install_github('YuLab-SMU/ggtree')"
RUN R -e "devtools::install_github('XiaoLuo-boy/ggheatmap')"

# ---- App setup ----
WORKDIR /home/app
COPY . /home/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

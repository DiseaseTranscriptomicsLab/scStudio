FROM r-base:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update -qq && apt-get -y upgrade --fix-missing && \
    apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e 'install.packages("shiny", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
RUN R -e 'install.packages("shinyjs", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'
          
# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Expose the port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

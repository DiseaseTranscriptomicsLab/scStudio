FROM r-base:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

RUN apt-get update -qq || (echo "apt-get update failed" && tail -n 50 /var/log/apt/*log && exit 1)
RUN apt-get install -y curl || (echo "Failed to install curl" && exit 1)
RUN apt-get install -y libcurl4-openssl-dev || (echo "Failed to install libcurl4-openssl-dev" && exit 1)
RUN apt-get install -y libssl-dev || (echo "Failed to install libssl-dev" && exit 1)
RUN apt-get install -y libxml2-dev || (echo "Failed to install libxml2-dev" && exit 1)

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

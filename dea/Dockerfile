FROM rocker/r-ver:4.1.2
MAINTAINER Marta Bica <mbica.compbio@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN apt-get install -y curl

RUN apt-get install -y libcurl4-openssl-dev

RUN apt-get install -y libssl-dev

RUN apt-get install -y libxml2-dev 

RUN apt-get install -y libgeos-dev

RUN apt-get install -y libglpk40 libglpk-dev

# Install renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"
      
# Set the working directory
WORKDIR /home/app

# Copy the app directory into the image
COPY . /home/app

# Set GitHub PAT (token with no scopes)
ENV GITHUB_PAT=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

# Restore the environment exactly
RUN R -e "renv::restore()"

# Expose the port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3838)"]

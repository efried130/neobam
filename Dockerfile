# STAGE 0 - Ubuntu packages and R repository
FROM ubuntu as stage0
RUN echo "America/New_York" | tee /etc/timezone \
	&& apt update \
	&& DEBIAN_FRONTEND=noninteractive apt install -y \
		build-essential \
		gcc \
		gfortran \
        locales \
		libcurl4-gnutls-dev \
		libfontconfig1-dev \
		libfribidi-dev \
		libgit2-dev \
		libharfbuzz-dev \
		libnetcdf-dev \
		libnetcdff-dev \
		libssl-dev \
		libtiff5-dev \
		libxml2-dev \
		tzdata \
		libv8-dev \
    && locale-gen en_US.UTF-8 \
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
	&& . /etc/lsb-release \
	&& echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${DISTRIB_CODENAME}-cran40/" >> /etc/apt/sources.list

# STAGE 1 - R and R packages
FROM stage0 as stage1
RUN apt update && apt -y install \
		r-base \
		r-base-dev \
	&& rm -rf /var/lib/apt/lists/* \
	&& /usr/bin/Rscript -e "install.packages('devtools')" \
	# && /usr/bin/Rscript -e "devtools::install_version('rstan',version='2.21.0', dependencies=TRUE, repos='http://cran.rstudio.com/')"	
	&& /usr/bin/Rscript -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.21.0-7.tar.gz')" \
	&& /usr/bin/Rscript -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.21.7.tar.gz', dependencies=TRUE)" \
	

	# && /usr/bin/Rscript -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.21.7.tar.gz', args = paste0('--library=', .libPaths()[2]))"
	&& /usr/bin/Rscript -e "install.packages('RNetCDF', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('parallel', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('doParallel', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('rjson', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('dplyr', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('tidyr', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('stringr', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('R.utils', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& /usr/bin/Rscript -e "install.packages('reticulate', dependencies=TRUE, repos='http://cran.rstudio.com/')"
	

# STAGE 2 - Python and python packages for S3 functionality
FROM stage1 as stage2
RUN apt update && apt -y install python3 python3-dev python3-pip python3-venv python3-boto3

# STAGE 3 set up I/O directories, copy geobamdata installer and R script
FROM stage2 as stage3
RUN mkdir -p /app/data/input \
	&& mkdir /app/data/output 
COPY ./run_neobam.R /app/
COPY ./neobam /app/neobam
COPY ./sos_read /app/sos_read/

# STAGE 3 - Execute algorithm
FROM stage3 as stage4
LABEL version="1.0" \
	description="Containerized neoBAM algorithm." \
	"confluence.contact"="ntebaldi@umass.edu" \
	"algorithm.contact"="cjgleason@umass.edu,cbrinkerhoff@umass.edu"
ENTRYPOINT [ "/usr/bin/Rscript",  "/app/run_neobam.R" ]
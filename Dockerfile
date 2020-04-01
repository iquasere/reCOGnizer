FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels anaconda \
&& git clone -b development https://github.com/iquasere/reCOGnizer.git \
&& conda install -c anaconda pandas \
&& conda install -c bioconda blast \
&& conda install -c anaconda lxml \
&& conda install -c anaconda openpyxl \
&& git clone https://github.com/marbl/Krona.git \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P /reCOGnizer \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/reCOGnizer/recognizer.py" ]
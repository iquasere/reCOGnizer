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
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P /reCOGnizer \
&& mkdir /Databases \
&& cd /Databases \
&& wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz \
&& tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp' \
&& rm cdd.tar.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz \
&& gunzip /Databases/cddid.tbl.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/reCOGnizer/recognizer.py" ]
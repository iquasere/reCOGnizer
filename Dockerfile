FROM ubuntu

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& git clone https://github.com/iquasere/reCOGnizer.git \
&& pip install -r requirements.txt \
&& cd /Databases \
&& wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz \
&& tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp' \
&& rm cdd.tar.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz \
&& gunzip /Databases/cddid.tbl.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog \
&& cd ../../.. \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P /reCOGnizer \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/reCOGnizer/recognizer.py" ]
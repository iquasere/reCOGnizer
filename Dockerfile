FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/reCOGnizer.git \
&& conda env create -f envs/environment.yml \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/reCOGnizer/recognizer.py" ]
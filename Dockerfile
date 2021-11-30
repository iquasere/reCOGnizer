FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/reCOGnizer.git \
&& conda install -c conda-forge -y mamba \
&& mamba env update --file reCOGnizer/envs/environment.yml --name base \
&& bash reCOGnizer/envs/ci_build.sh \
&& conda clean --all -y \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "share/recognizer.py" ]
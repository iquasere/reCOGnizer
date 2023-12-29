FROM continuumio/miniconda3:23.3.1-0

RUN git clone https://github.com/iquasere/reCOGnizer.git \
&& conda install -c conda-forge mamba libarchive=3.6.2=h039dbb9_1 \
&& mamba env update --file reCOGnizer/envs/environment.yml --name base \
&& bash reCOGnizer/envs/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "share/recognizer.py" ]
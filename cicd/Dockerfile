FROM continuumio/miniconda3:23.3.1-0

RUN git clone https://github.com/iquasere/reCOGnizer.git \
&& conda install -c conda-forge mamba \
&& mamba env update --file reCOGnizer/cicd/environment.yml --name base \
&& bash reCOGnizer/cicd/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "share/recognizer.py" ]
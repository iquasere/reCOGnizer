sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev --no-install-recommends
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels anaconda
conda install -y -c anaconda pandas lxml openpyxl
conda install -y -c bioconda blast krona
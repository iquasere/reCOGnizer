sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev --no-install-recommends
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels anaconda
conda install -y -c anaconda pandas
conda install -y -c bioconda blast
conda install -y -c anaconda lxml
conda install -y -c anaconda openpyxl
cd reCOGnizer
git clone https://github.com/marbl/Krona.git
wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl

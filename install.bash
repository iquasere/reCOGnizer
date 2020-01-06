sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev --no-install-recommends
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels anaconda
conda install -c anaconda pandas
conda install -c bioconda blast
conda install -c anaconda lxml
wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P /reCOGnizer
mkdir Databases
cd Databases
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp'
rm cdd.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
gunzip /Databases/cddid.tbl.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog
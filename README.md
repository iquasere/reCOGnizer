# reCOGnizer

A tool for domain based annotation with the COG database.

## Features

reCOGnizer takes a FASTA file as input and produces two main outputs into the output directory:
* protein2cog.xlsx, an Excel file assigning COGs to the proteins
* cog_quantification.tsv and respective Krona representation, which describes the functinoal landscape of the proteins in the input file

## Installation

To install reCOGnizer, simply clone this repository and run install.bash!
```
git clone https://github.com/iquasere/reCOGnizer.git
sudo reCOGnizer/install.bash
```

## Usage

reCOGnizer needs an input file, but that is all it needs!
```
usage: recognizer.py [-h] -f [FILE [FILE ...]] [-t THREADS] [-o OUTPUT]
                     [-db DATABASE] [-seqs MAX_TARGET_SEQS]

reCOGnizer - a tool for domain based annotation with the COG database

optional arguments:
  -h, --help            show this help message and exit
  -f [FILE [FILE ...]], --file [FILE [FILE ...]]
                        Fasta file with protein sequences for annotation
  -t THREADS, --threads THREADS
                        Number of threads for reCOGnizer to use. Default is
                        number of CPUs available minus 2.
  -o OUTPUT, --output OUTPUT
                        Output directory
  -db DATABASE, --database DATABASE
                        Basename of COG database for annotation
  -seqs MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                        Number of maximum identifications for each protein.
                        Default is 1.
```

## Docker

reCOGnizer already has its own image! To use it, just pull the image and run it!.
```
docker pull iquasere/recognizer:latest
docker run iquasere/recognizer [arguments]
```
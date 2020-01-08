# reCOGnizer

A tool for domain based annotation with the COG database.

## Features

reCOGnizer is a user-friendly implementation of protein functional identification using COG database. It builds a splitted version of the COG database with which PSI-BLAST can run in multithread, significantly speeding the time intensive step of protein annotation. After COG assignment to proteins, reCOGnizer makes use of perl2cog to convert CDD IDs to their respective COGs, before organizing those COGs into a relational table of protein to COG, with the inclusion of the three levels of functional classification from COG.

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
The simplest way to run reCOGnizer just needs the specification of the fasta file and an output directory
```
recognizer.py -f input_file.fasta -o output_folder
```

## Docker

reCOGnizer already has its own image! To use it, just pull the image and run it!.
```
docker pull iquasere/recognizer:latest
docker run -it -v absolute/path/to/fasta_folder:/input_folder /absolute/path/to/output_folder:/output_folder --rm iquasere/recognizer -f /input_folder/input_file.fasta -o /output_folder [other arguments]
```
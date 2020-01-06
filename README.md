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

## Docker

reCOGnizer already has its own image! To use it, just pull the image and run it!.
```
docker pull iquasere/recognizer:latest
docker run iquasere/recognizer [arguments]
```
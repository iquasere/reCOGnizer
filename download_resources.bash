#!/bin/bash

if [ "$1" != "" ]; then
    cd $1
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
    tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp'
    rm cdd.tar.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    gunzip /Databases/cddid.tbl.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog
else
    echo "Positional parameter 1 is empty"
fi
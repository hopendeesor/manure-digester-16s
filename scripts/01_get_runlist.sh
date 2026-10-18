#!/usr/bin/env bash
base="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,sample_accession,sample_alias,read_count,fastq_ftp,fastq_md5"
curl -s "$base?accession=PRJEB6969&result=read_run&fields=$fields&format=tsv" > data/processed/runlist.tsv

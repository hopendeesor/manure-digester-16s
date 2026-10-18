#!/usr/bin/env bash
set -e
mkdir -p data/raw/ref
f=silva_nr99_v138.2_toGenus_trainset.fa.gz
u="https://zenodo.org/records/14169026/files/$f?download=1"
curl -L -s -f -o data/raw/ref/$f "$u"
md5sum data/raw/ref/$f
ls -l data/raw/ref/$f

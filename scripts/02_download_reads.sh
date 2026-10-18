#!/usr/bin/env bash
set -e
list=data/processed/runlist.tsv
mkdir -p data/raw
while IFS=$'\t' read run samp alias n url md5; do
  f="data/raw/$run.fastq.gz"
  if ! echo "$md5  $f" | md5sum -c --quiet - 2>/dev/null; then
    curl -s -f -o "$f" "https://$url"
    echo "$md5  $f" | md5sum -c --quiet -
  fi
  echo "ok $alias $run"
done < <(tail -n +2 $list | tr -d '\r')

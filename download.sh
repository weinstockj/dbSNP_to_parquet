#!/bin/bash

set -e

echo "Starting dbSNP download..."

wget -c -t 0 --retry-connrefused --waitretry=30 \
  https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz

wget -c -t 0 --retry-connrefused --waitretry=30 \
  https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

echo "Download completed successfully!"

#!/bin/bash
snakemake -p --profile slurm -j 23 1>stdout.log 2>stderr.log

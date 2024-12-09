
CONDA="/home/jweins17/data-abattle4/jweins17/mambaforge/etc/profile.d/conda.sh"
MAMBA="/home/jweins17/data-abattle4/jweins17/mambaforge/etc/profile.d/mamba.sh"

CHROMS = [f"chr{x}" for x in range(1, 23)] + ['chrX']

print(CHROMS)

rule all:
    input:
        "dbSNP_156.bcf",
        expand("dbSNP_156.{CHROM}.bcf", CHROM = CHROMS),
        expand("dbSNP_156.{CHROM}.lookup.parquet", CHROM = CHROMS)

rule convert_contig:
    input:
        "GCF_000001405.40.gz"
    output:
        "dbSNP_156.bcf"
    threads: 3 
    resources:
        mem = "6G",
        partition = "shared",
        time = "07:30:00"
    shell:
        """
        source {CONDA}
        source {MAMBA}

        mamba activate bcftools
            
        bcftools view -i "INFO/VC='SNV'" {input} | 
            bcftools annotate --rename-chrs contig_map.tsv -Ob --threads {threads} > {output}

        bcftools index {output}
        """

rule subset_BCF:
    input:
        "dbSNP_156.bcf"
    output:
        "dbSNP_156.{CHROM}.bcf"
    threads: 2 
    resources:
        mem = "3G",
        partition = "shared",
        time = "03:30:00"
    shell:
        """
        source {CONDA}
        source {MAMBA}

        mamba activate bcftools
        bcftools view -Ob --threads {threads} {input} {wildcards.CHROM} > {output}
        bcftools index {output}
        """
        
rule convert_parquet:
    input:
        "dbSNP_156.{CHROM}.bcf"
    output:
        "dbSNP_156.{CHROM}.lookup.parquet"
    resources:
        mem = "8G",
        partition = "shared",
        time = "05:30:00"
    shell:
        """
        source {CONDA}
        source {MAMBA}

        mamba activate cyvcf2

        python convert_to_parquet.py {input} {output}
        """
    

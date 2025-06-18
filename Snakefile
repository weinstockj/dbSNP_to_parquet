
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

CHROMS = [f"chr{x}" for x in range(1, 23)] + ['chrX']

print(CHROMS)

rule all:
    input:
        f"{config['output_dir']}/dbSNP_156.bcf",
        expand(f"{config['output_dir']}/dbSNP_156.{{CHROM}}.bcf", CHROM = CHROMS),
        expand(f"{config['output_dir']}/dbSNP_156.{{CHROM}}.lookup.parquet", CHROM = CHROMS)

rule convert_contig:
    input:
        "GCF_000001405.40.gz"
    output:
        f"{config['output_dir']}/dbSNP_156.bcf"
    threads: config['resources']['convert_contig']['threads']
    resources:
        mem = config['resources']['convert_contig']['mem'],
        partition = config['cluster']['partition'],
        time = config['resources']['convert_contig']['time']
    shell:
        """
        bcftools view -i "INFO/VC='SNV'" {input} | 
            bcftools annotate --rename-chrs contig_map.tsv -Ob --threads {threads} > {output}

        bcftools index {output}
        """

rule subset_BCF:
    input:
        f"{config['output_dir']}/dbSNP_156.bcf"
    output:
        f"{config['output_dir']}/dbSNP_156.{{CHROM}}.bcf"
    threads: config['resources']['subset_BCF']['threads']
    resources:
        mem = config['resources']['subset_BCF']['mem'],
        partition = config['cluster']['partition'],
        time = config['resources']['subset_BCF']['time']
    shell:
        """
        bcftools view -Ob --threads {threads} {input} {wildcards.CHROM} > {output}
        bcftools index {output}
        """
        
rule convert_parquet:
    input:
        f"{config['output_dir']}/dbSNP_156.{{CHROM}}.bcf"
    output:
        f"{config['output_dir']}/dbSNP_156.{{CHROM}}.lookup.parquet"
    resources:
        mem = config['resources']['convert_parquet']['mem'],
        partition = config['cluster']['partition'],
        time = config['resources']['convert_parquet']['time']
    shell:
        """
        source .venv/bin/activate

        python convert_to_parquet.py {input} {output}
        """
    

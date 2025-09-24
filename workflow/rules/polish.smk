#  Bring assembly to local disk to avoid
# I/O bottleneck
rule bring_assembly:
    input:
        draft = f'{base_dir}/assembly/flye/{{sample}}/assembly.fasta',
    output:
        draft = temp("data/assembly/{sample}.fasta"),
        index = temp("data/assembly/{sample}.fasta.fai")
    conda:
        "../envs/contigging.yaml"
    shell:
        """
        cp {input.draft} {output.draft} && \
        samtools faidx {output.draft}
        """
# Use the input ONT reads to polish
# the assembly
rule polish_medaka:
    input:
        basecalls = "data/raw_reads/{sample}.fastq.gz",
        draft = "data/assembly/{sample}.fasta",
        index = "data/assembly/{sample}.fasta.fai"
    output:
        consensus = f"{base_dir}/polishing/medaka/{{sample}}/consensus.fasta",
        mappings = f"{base_dir}/polishing/medaka/{{sample}}/calls_to_draft.bam"
    threads:
        resources["medaka"]["threads"]
    log:
        'logs/{sample}/polishing/medaka.log'
    conda:
        config["polish"]["medaka"]["env"]
    shell:
        """
        medaka_consensus \
            -i {input.basecalls} \
            -d {input.draft}\
            -o tmp_polishing/{wildcards.sample} \
            -t {threads} \
            -m {config[polish][medaka][model]} 2> {log} && \
        cp -r tmp_polishing/{wildcards.sample} {base_dir}/polishing/medaka/ && \
        rm -r tmp_polishing/{wildcards.sample}
        """



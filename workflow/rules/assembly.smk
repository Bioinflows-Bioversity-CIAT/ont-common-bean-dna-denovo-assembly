""" Run assembly using flye Assembler.
 store the assembly in remote disk """
rule flye_denovo_assembly:
    input:
        raw_reads = get_raw_reads
    output:
        draft = f"{base_dir}/assembly/flye/{{sample}}/assembly.fasta",
        graph = multiext(f"{base_dir}/assembly/flye/{{sample}}/assembly_graph",".gfa",".gv"),
        info = f"{base_dir}/assembly/flye/{{sample}}/assembly_info.txt",
    log:
        'logs/{sample}/assembly/flye.log'
    conda:
        "../envs/contigging.yaml"
    threads:
        resources['flye']['threads']
    shell:
        """
        flye --nano-raw {input.raw_reads} \
        --out-dir {base_dir}/assembly/flye/{wildcards.sample} \
        --threads {threads} 2> {log}
        """

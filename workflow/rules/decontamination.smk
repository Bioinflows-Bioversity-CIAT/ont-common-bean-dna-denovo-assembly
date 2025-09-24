rule blast:
    """Perform a BLAST search on the assembled contigs to filter out
     potential contaminants by identifying taxonomic matches."""
    input:
      curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
    output:
        out_blast = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.asm.vs.nt.max10.1e25.blastn.out",
    log:
        'logs/{sample}/decontamination/blast.log'
    conda:
        "../envs/decontamination.yaml"
    threads:
        resources["decontamination"]["blastn_threads"]
    resources: 
        mem_mb = resources["decontamination"]["blastn_mem_mb"]
    params: 
        nt_db = config['nt_db']
    shell:
        """
        (blastn -query {input.curated_draft} \
        -db {params.nt_db} \
        -outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {threads} \
        -out {output.out_blast}) 2>> {log}
        """

rule coverage_estimate:
    input:
        curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
        mappings = f"{base_dir}/polishing/medaka/{{sample}}/calls_to_draft.bam",
    output:
        cov = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.coverage",
    conda:
        "../envs/decontamination.yaml"
    log:
        'logs/{sample}/decontamination/map2cov.log'
    shadow:
        "shallow"
    shell:
        """
        blobtools map2cov -i {input.curated_draft} -b {input.mappings} -o ./coverage_tmp 2> {log} && \
        mv ./coverage_tmp.calls_to_draft.bam.cov {output.cov}
        """

rule blobtools_create_db:
    input:
        curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
        out_blast = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.asm.vs.nt.max10.1e25.blastn.out",
        cov = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.coverage",
    output:
        f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobDB.json",
    conda:
        "../envs/decontamination.yaml"
    params:
        taxdump_db = config['taxdump_db']
    shadow:
        "shallow"
    log:
        'logs/{sample}/decontamination/bloptools_create.log'
    shell:
        """
        blobtools create -i {input.curated_draft} -t {input.out_blast} \
        -c {input.cov} \
        --nodes {params[taxdump_db]}/nodes.dmp \
        --names {params[taxdump_db]}/names.dmp \
        -x bestsumorder -o tmp_db 2> {log} && \
        mv tmp_db.blobDB.json {output}
        """

rule blobtools_plot_table:
    input:
        blob_db = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobDB.json",
    output:
        read_cov_plot = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobplot.read_cov.cov0.png",
        blob_plot = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobplot.cov0.png",
        table = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobplot.stats.txt",
    conda:
        "../envs/decontamination.yaml"
    log:
        'logs/{sample}/decontamination/blobplot_table.log'
    shadow:
        "shallow"
    shell:
        """
        blobtools plot -i {input.blob_db} --sort count --hist count -x bestsumorder -o plot 2> {log} && \
        blobtools view -i {input.blob_db} --hits --rank all -x bestsumorder -o table 2>> {log} && \
        mv table.{wildcards.sample}*.table.txt {output.table} && \
        mv *.read_cov.cov0.png {output.read_cov_plot} && \
        mv *.blobplot.cov0.png {output.blob_plot}
        """

rule get_contamination:
    input:
        table = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.blobplot.stats.txt",
        curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
    output:
        contamination_contigs = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.contamination.txt",
        decontaminated_draft = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.decontaminated.fa",
    conda:
        "../envs/decontamination.yaml"
    log:
        'logs/{sample}/decontamination/seqkit_grep.log'
    shell:
        """
        cat {input.table}| tail -n+12| \
        awk '{{if($10 != "{config[target_taxa_id]}" && $10 != "no-hit"){{print $1}}}}' > {output.contamination_contigs} && 
        if [ "$(wc -l < {output.contamination_contigs})" -gt 1 ]; then
        seqkit grep -f {output.contamination_contigs} -v -o {output.decontaminated_draft} {input.curated_draft} 2> {log}
        else
        cp {input.curated_draft} {output.decontaminated_draft}
        fi
        """


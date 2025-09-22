rule fix_space_repeat_masker:
    """ 
    RepeatMasker produces an gff
    with TE annotations but integers have an extra space
    after the \t this rule fix it.
    """
    input:
        f"{base_dir}/annotation/repeats/RepeatMasker/{{sample}}/{{sample}}.fa.out.gff"
    output:
        f"{base_dir}/annotation/genes/Maker/{{sample}}/RepeatMasker/{{sample}}.gff"
    shell:
        """
        cat {input}| sed 's/\t[[:space:]]\+/\t/g' > {output}
        """

rule hardmask_assembly:
    """ 
    Maker predicted a more rasonable number of genes
    when hardmask the assembly
    """
    input:
        assembly = rules.get_contamination.output.decontaminated_draft,
        repeats = rules.fix_space_repeat_masker.output
    output:
        f"{base_dir}/annotation/genes/Maker/{{sample}}/RepeatMasker/{{sample}}_hardmasked.fa"
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        bedtools maskfasta -fi {input.assembly} \
        -bed {input.repeats} -fo {output}
        """

rule set_maker_config_files:
    """ 
    This rule initialize the maker2 configuration
    files and make some modifications on the parameters such as:
    Indicate the path of a priori transcripts & proteins sequences.
    Similarly, specify the name of augustus species model named:
    phaseouls_vulgaris_UI-111, keep_preds=1 and specify blast parameters:
    depth_blastx = 10, depth_blastn = 10 and depth_tblast = 10
    """
    input:
        assembly_masked = rules.hardmask_assembly.output
    output:
        multiext(f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/maker_',
         'bopts.ctl',
          'evm.ctl',
            'exe.ctl',
             'opts.ctl')
    params:
        opts_maker_script = "./workflow/scripts/maker_opts_converter.awk",
        bopts_maker_script = "./workflow/scripts/maker_bopts_converter.awk"
    shadow:
        'shallow'
    conda:
        config["annotation"]["maker_env"]
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        maker -CTL && \
        awk -v est_path="{config[annotation][reference_sequences][transcriptome]}" \
            -v protein_path="{config[annotation][reference_sequences][proteins]}" \
            -v genome_path="{input.assembly_masked}" \
            -v species_name="{config[annotation][augustus_specie_name]}" -f {params.opts_maker_script} \
            ./maker_opts.ctl > ./tmp.ctl && \
        rm ./maker_opts.ctl && \
        mv ./tmp.ctl ./maker_opts.ctl && \
        awk -f {params.bopts_maker_script} ./maker_bopts.ctl > ./tmp.ctl && \
        mv ./tmp.ctl ./maker_bopts.ctl && \
        mv ./maker_* {base_dir}/annotation/genes/Maker/{wildcards.sample}/maker2_run/ 
        """

rule run_maker:
    input:
        rules.set_maker_config_files.output
    output:
        index = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.maker.output/{{sample}}_hardmasked_master_datastore_index.log'
    conda:
        config["annotation"]["maker_env"]
    threads:
        resources["annotation"]["genes"]["maker2_threads"]
    log:
        f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/maker2_run.log'
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        {config[mpiexec]} -n {threads} maker -c 1 -q {base_dir}/annotation/genes/Maker/{wildcards.sample}/maker2_run/maker_opts.ctl \
        {base_dir}/annotation/genes/Maker/{wildcards.sample}/maker2_run/maker_bopts.ctl \
        {base_dir}/annotation/genes/Maker/{wildcards.sample}/maker2_run/maker_exe.ctl 2> {log} && \
        mv ./{wildcards.sample}_hardmasked.maker.output {base_dir}/annotation/genes/Maker/{wildcards.sample}/maker2_run/
        """

rule collect_gff_annotations:
    input:
        index_path = rules.run_maker.output.index
    output:
        merged_gff = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.gff'
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        gff3_merge -s -d {input.index_path} > \
            {output.merged_gff}
        """

rule get_annotation_fasta:
    input:
        index_path = rules.run_maker.output.index
    output:
        proteins = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.maker.proteins.fasta',
        transcripts = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.maker.transcripts.fasta'
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        fasta_merge -d {input.index_path} && \
            mv ./{wildcards.sample}_hardmasked.all.maker.proteins.fasta {output.proteins} && \
            mv ./{wildcards.sample}_hardmasked.all.maker.transcripts.fasta {output.transcripts}
        """

rule rename_annotation_names:
    input:
        merged_gff = rules.collect_gff_annotations.output.merged_gff
    output:
        map_ids = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.map'
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        maker_map_ids --prefix {wildcards.sample}- --justify 6 {input.merged_gff} > {output.map_ids}
        """

rule rename_proteins_transcripts:
    input:
        proteins = rules.get_annotation_fasta.output.proteins,
        transcripts = rules.get_annotation_fasta.output.transcripts,
        merged_gff = rules.collect_gff_annotations.output.merged_gff,
        map_ids = rules.rename_annotation_names.output.map_ids
    output:
        proteins = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.maker.renamed.proteins.fasta',
        transcripts = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.maker.renamed.transcripts.fasta',
        merged_gff = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.renamed.gff',
    conda:
        config["annotation"]["maker_env"]
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        cp {input.proteins} ./proteins.fa && \
        cp {input.transcripts} ./transcripts.fa && \
        cp {input.merged_gff} ./annotations.gff && \
        map_gff_ids {input.map_ids} ./annotations.gff && \
        map_fasta_ids {input.map_ids} ./proteins.fa && \
        map_fasta_ids {input.map_ids} ./transcripts.fa && \
        mv ./annotations.gff {output.merged_gff} && \
        mv ./transcripts.fa {output.transcripts} && \
        mv ./proteins.fa {output.proteins}
        """

rule filter_annotations_default:
    input:
        merged_gff = rules.rename_proteins_transcripts.output.merged_gff
    output:
        merged_gff = f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.renamed.default.gff'
    conda:
        config["annotation"]["maker_env"]
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        quality_filter.pl -d {input.merged_gff} > {output.merged_gff}
        """

rule blast_uniprot:
    input:
        proteins = rules.rename_proteins_transcripts.output.proteins
    output:
        f'{base_dir}/annotation/genes/Maker/{{sample}}/functional/uniprot/{{sample}}_hardmasked.all.maker2.renamed.maker2uniprot.blastp',
    log:
        f'{base_dir}/annotation/genes/Maker/{{sample}}/functional/uniprot/blastp.log',
    conda:
        config["annotation"]["maker_env"]
    threads:
        resources["annotation"]["genes"]["blastp_uniprot_threads"]
    shell:
        """
        blastp -db {config[annotation][reference_sequences][uniprot_db]} \
        -query {input.proteins} \
        -out {output} \
        -evalue 0.000001 \
        -outfmt 6 \
        -max_hsps 1 \
        -num_alignments 1 \
        -seg yes \
        -lcase_masking \
        -soft_masking true \
        -num_threads {threads} 2> {log}
        """

rule annotate_functional_uniprot:
    input:
        merged_gff = rules.filter_annotations_default.output.merged_gff,
        uniprot_blastp = rules.blast_uniprot.output,
    output:
        f'{base_dir}/annotation/genes/Maker/{{sample}}/maker2_run/{{sample}}_hardmasked.all.renamed.default.uniprot.gff'
    conda:
        config["annotation"]["maker_env"]
    shadow:
        "shallow"
    shell:
        """
        export PATH=$PATH:{config[annotation][maker_bin]} && \
        maker_functional_gff {config[annotation][reference_sequences][uniprot_db]} {input.uniprot_blastp} {input.merged_gff} > {output}
        """


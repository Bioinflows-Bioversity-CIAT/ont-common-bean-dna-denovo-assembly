rule make_blast_db:
    """ Create a blast db with the decontaminated
    draft for search the genetic map markers """
    input:
        decontaminated_draft = rules.get_contamination.output.decontaminated_draft,
    output:
        db = multiext( f"{base_dir}/scaffolding/{{sample}}/blastdb/{{sample}}",
            ".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto")
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        makeblastdb -in {input.decontaminated_draft} \
        -dbtype nucl \
        -input_type fasta -blastdb_version 5 -parse_seqids \
        -out {base_dir}/scaffolding/{wildcards.sample}/blastdb/{wildcards.sample}
        """

rule map_markers_to_assembly:
    input:
        markers_sequence = config['scaffolding']['linkage_map_markers'],
        db = multiext( f"{base_dir}/scaffolding/{{sample}}/blastdb/{{sample}}",
            ".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto")
    output:
        f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_markers.blastn"
    params:
        evalue = 1e-10,
        idy = 95,
        word_size =11,
        max_target_seqs = 1
    conda:
        "../envs/decontamination.yaml"
    threads:
        resources["scaffolding"]["blastn_linkage_map_threads"]
    log:
        'logs/{sample}/scaffolding/{sample}_mapmarkers.log'
    shell:
        """
        blastn -db {base_dir}/scaffolding/{wildcards.sample}/blastdb/{wildcards.sample} \
        -query {input.markers_sequence} \
        -evalue {params.evalue} \
        -perc_identity {params.idy} \
        -word_size {params.word_size} \
        -max_target_seqs {params.max_target_seqs} \
        -num_threads {threads} \
        -max_hsps 1 \
        -outfmt 6 -out {output} 2> {log}
        """

rule get_genetic_map_bed:
    input:
        genetic_map = config['scaffolding']['linkage_map'],
        blast_out = f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_markers.blastn"
    output:
        f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_map.bed"
    params:
        script = "./workflow/scripts/convert_blast2csv.R",
        cutoff = 100,
    conda:
        "../envs/r-env.yaml"
    shell:
        """
        Rscript {params.script} -blastn {input.blast_out} \
        -cutoff {params.cutoff} -map {input.genetic_map} \
        -out {output}
        """

rule make_ref_blast_db:
    """ Create a blast db with the reference
    peptide sequences to infer the syntenylogs """
    input:
        ref_prots = get_prot_ref_assembly
    output:
        db = multiext( f"{base_dir}/scaffolding/ref_blastdb/{{ref_assembly}}",
            ".pdb", ".phr",".pin", ".pot", ".psq", ".ptf", ".pto")
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        makeblastdb -in {input.ref_prots} \
        -dbtype prot \
        -input_type fasta -blastdb_version 5 -parse_seqids \
        -out {base_dir}/scaffolding/ref_blastdb/{wildcards.ref_assembly}
        """

rule blastp_job_ref_assembly:
  input:
    query = f'{base_dir}/quality/genes/{{sample}}/{{sample}}_hardmasked.all.renamed.default_proteins.fa',
    blastpdb = multiext(f"{base_dir}/scaffolding/ref_blastdb/{{ref_assembly}}",
        ".pdb", ".phr", ".pin", ".pot", ".psq", ".ptf", ".pto")
  output:
    f"{base_dir}/scaffolding/{{sample}}/synteny/{{ref_assembly}}/{{sample}}_2_{{ref_assembly}}.blastp.txt"
  log:
    'logs/{sample}/synteny/blastp/{ref_assembly}.log'
  threads:
    resources["scaffolding"]["blastp_ref_proteins_threads"]
  conda:
        "../envs/decontamination.yaml"
  shell:
    """
    blastp -query {input.query} \
    -db {base_dir}/scaffolding/ref_blastdb/{wildcards.ref_assembly} \
    -evalue 1e-5 \
    -num_threads {threads} \
    -outfmt 6 \
    -num_alignments 5 \
    -out {output} 2> {log}
    """

rule extract_gene_ids:
  input:
      gff = rules.annotate_functional_uniprot.output
  output:
      f"{base_dir}/scaffolding/{{sample}}/synteny/{{sample}}_gene_ids.txt"
  shell:
    """
    cat {input}| awk '{{if($3=="mRNA"){{print $9 }}}}'|cut -f1 -d';'| sed 's/ID=//g' > {output}
    """

rule get_gene_directory_query:
  input:
      gff = rules.annotate_functional_uniprot.output
  output:
      f"{base_dir}/scaffolding/{{sample}}/synteny/{{sample}}_gene_directory.txt"
  shell:
    """
    cat {input}| awk '{{if($3=="mRNA"){{print $1,$4,$5,$9 }}}}'|cut -f1 -d';'| sed 's/ID=//g'| sed 's/\s/\t/g' > {output}
    """

rule get_gene_directory_ref:
  input:
      get_gene_annot_ref_assembly
  output:
      f"{base_dir}/scaffolding/{{sample}}/synteny/{{ref_assembly}}/{{ref_assembly}}_gene_directory.txt"
  shell:
    """
    zcat {input}| grep -v "^#" | awk '{{if($3=="gene"){{print $1,$4,$5,$9 }}}}'| sed 's/;/\\t/g'| sed 's/\s/\t/g'| cut -f1,2,3,5| sed 's/Name=//g' > {output}
    """

rule format_blastp2dagchainer:
    input:
        blastp = rules.blastp_job_ref_assembly.output,
        r_dict = rules.get_gene_directory_ref.output,
        qry_dict = rules.get_gene_directory_query.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/synteny/{{sample}}_{{ref_assembly}}.matchlist.txt"
    params:
        script = "./workflow/scripts/blastp2dagchainer.R",
        suffix = lambda wildcards: get_protein_suffix(wildcards)
    conda:
        "../envs/r-env.yaml"
    shell:
        """
        Rscript {params.script} -blastp {input.blastp} \
        -ref_dic {input.r_dict} \
        -qry_dic {input.qry_dict} \
        -p {params.suffix} \
        -out {output}
        """

rule dagchainer_filter:
    input:
        rules.format_blastp2dagchainer.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/synteny/{{sample}}_{{ref_assembly}}.filter.txt"
    params:
        neighbor_genes = config["scaffolding"]["dagchainer"]["neighbor_genes_filter"],
        script =  config["scaffolding"]["dagchainer"]["filter_script"]
    conda:
        "../envs/dagchainer.yaml"
    shell:
        """
        {params.script} {params.neighbor_genes} < {input} > {output}
        """

rule dagchainer_job:
    input:
        rules.dagchainer_filter.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/synteny/{{sample}}_{{ref_assembly}}.filter.txt.aligncoords"
    log:
        'logs/{sample}/synteny/dagchainer/{sample}_{ref_assembly}.filter.aligncoords.log'
    conda:
        "../envs/dagchainer.yaml"
    shell:
        """
        run_DAG_chainer.pl -i {input} -Z 12 -D 10 -g 1 -A 5 1> {log} 
        """

rule get_syntenylogs:
    input:
        synteny_blocks = rules.dagchainer_job.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_sinteny_blocks.anchors"
    shell:
        """
        cat {input.synteny_blocks}| grep -v '^#'| awk "{{print \$6,\$2,\$9}}"| sed 's/\s/\t/g' > {output}
        """

rule get_query_gene_set_bed:
    input:
        gff = rules.annotate_functional_uniprot.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_query_gene_set.bed"
    shell:
        """
        cat {input.gff}| awk '{{if($3=="mRNA"){{print $1,$4,$5,$7,$9 }}}}' |\
        cut -f1 -d';'| sed 's/ID=//g' |\
        sed 's/\s/\t/g'|\
        awk '{{print $1,$2,$3,$5, 1000,$4}}' | sed 's/\s/\t/g' > {output}
        """

rule get_ref_gene_set_bed:
    input:
        gff = get_gene_annot_ref_assembly
    output:
        f"{base_dir}/scaffolding/{{sample}}/allmaps/{{ref_assembly}}_gene_set.bed"
    shell:
        """
        zcat {input}| grep -v "^#" | awk '{{if($3=="gene"){{print $1,$4,$5,$7,$9 }}}}' |\
        sed 's/;/\\t/g'| sed 's/\s/\t/g'|sed 's/Name=//g'| cut -f1,2,3,4,6 |\
        awk '{{print $1,$2,$3,$5, 1000,$4}}'| sed 's/\s/\t/g' > {output}
        """

rule get_synteny_map_bed:
    input:
        qbed = rules.get_query_gene_set_bed.output,
        sbed = rules.get_ref_gene_set_bed.output,
        syntenylogs = rules.get_syntenylogs.output
    output:
        f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_synteny.bed"
    conda:
        config["scaffolding"]["allmaps_env"]
    shell:
        """
        python -m jcvi.assembly.syntenypath bed {input.syntenylogs} --qbed {input.qbed} --sbed {input.sbed} --scale 100000 -o {output}
        """
rule merge_beds:
    input:
        synteny_bed = rules.get_synteny_map_bed.output,
        genetic_map_bed = rules.get_genetic_map_bed.output
    output:
        bed = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_merged.bed",
        weights = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{ref_assembly}}_weights.txt",
    conda:
        config["scaffolding"]["allmaps_env"]
    shell:
        """
        python -m jcvi.assembly.allmaps mergebed {input.genetic_map_bed} {input.synteny_bed} -o {output.bed} && \
        cp weights.txt {output.weights} && rm weights.txt 
        """

rule all_maps_scaffolding:
    input:
        bed = rules.merge_beds.output.bed,
        decontaminated_draft = rules.get_contamination.output.decontaminated_draft,
        weights = rules.merge_beds.output.weights
    output:
        summary = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_merged.summary.txt",
        chain = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_merged.chain",
        scaffolded_asm = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_merged.fasta",
    params:
        ngen = 1000
    conda:
        config["scaffolding"]["allmaps_env"]
    threads:
        resources["scaffolding"]["allmaps_threads"]
    log:
        'logs/allmaps_{{sample}}_{{ref_assembly}}.log'
    shell:
        """
        python -m jcvi.assembly.allmaps path {input.bed} {input.decontaminated_draft} -w {input.weights} --cpus {threads} --ngen {params.ngen} \
        2> {log} && \
        mkdir {config[base_dir]}/scaffolding/{wildcards.sample}/allmaps/{wildcards.ref_assembly}_plots && \
        mv *.pdf {config[base_dir]}/scaffolding/{wildcards.sample}/allmaps/{wildcards.ref_assembly}_plots/
        """

rule migrate_maker_annotations:
    input:
        gff = rules.annotate_functional_uniprot.output,
        chain = rules.all_maps_scaffolding.output.chain,
    output:
        gff = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_decontaminated_hardmasked.all.renamed.default.uniprot_maker_lifted.gff",
        unmapped = f"{base_dir}/scaffolding/{{sample}}/allmaps/{{sample}}_{{ref_assembly}}_decontaminated_hardmasked.all.renamed.default.uniprot_maker_lifted.unmapped",
    conda:
        config["annotation"]["repeats_env"]
    shell:
        """
        cat {input.gff} | awk '{{if($2 == "maker"){{print $0}}}}' \
        > {config[base_dir]}/scaffolding/{wildcards.sample}/allmaps/{wildcards.sample}_decontaminated_hardmasked.all.renamed.default.uniprot_maker.gff && \
        {config[liftOver]} -gff {config[base_dir]}/scaffolding/{wildcards.sample}/allmaps/{wildcards.sample}_decontaminated_hardmasked.all.renamed.default.uniprot_maker.gff \
        {input.chain} {output.gff} {output.unmapped}
        """

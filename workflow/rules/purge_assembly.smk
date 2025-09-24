# with raw reads allign to draft
# the resultig bam is the input for polish
rule overlaps_for_purge: 
    input:
        polish_draft = f"{base_dir}/polishing/medaka/{{sample}}/consensus.fasta",
        reads = get_read_list_per_sample
    output:
        overlaps = f"{base_dir}/purge/purge_haplotigs/{{sample}}/overlaps/{{sample}}.overlaps.bam"
    conda:
        "../envs/purge_haplotigs.yaml"
    threads:
        resources["purge_haplotigs"]["minimap2_threads"]
    shell:
        """
        minimap2 -t {threads} -a -x map-ont \
         {input.polish_draft} {input.reads}| samtools sort -o {output.overlaps} --write-index - 
        """
# Generate with readhist a coverage
# histogram to define the thresholds
# of junk, haploid and diploid peaks
rule get_read_histogram:
    input:
        overlaps =  f"{base_dir}/purge/purge_haplotigs/{{sample}}/overlaps/{{sample}}.overlaps.bam",
        polish_draft = f"{base_dir}/polishing/medaka/{{sample}}/consensus.fasta",
    output:
        gen_coverage_stats =  f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}_polished.gencov",
        histogram = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}_polished_histogram.png",
    threads:
        resources["purge_haplotigs"]["readhist_threads"]
    conda:
        "../envs/purge_haplotigs.yaml"
    shell:
        """
        purge_haplotigs readhist -t {threads} \
            -b {input.overlaps} \
            -g {input.polish_draft} && \
            mv {wildcards.sample}.overlaps.bam.200.gencov {output.gen_coverage_stats} && \
            mv {wildcards.sample}.overlaps.bam.histogram.200.png {output.histogram}
        """
# With contigcov classify each
# contig in categories junk or suspect
rule coverage_stats_for_purge:
    input:
        gen_coverage_stats =  f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}_polished.gencov",
    output:
        contig_coverage =  f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}_polished_contig_coverage.csv",
    params:
        cutoffs = lambda wildcards: get_read_depth_cutoffs(wildcards),
    conda:
        "../envs/purge_haplotigs.yaml"
    log:
        'logs/{sample}/purge/coverage_stats.log'
    shell:
        """
        purge_haplotigs contigcov -i {input.gen_coverage_stats} \
            -l {params[cutoffs][l]} \
            -m {params[cutoffs][m]} \
            -h {params[cutoffs][h]} && \
        mv coverage_stats.csv {output.contig_coverage}
        """
# Identify the primary and haplotigs
# generate a fasta with only primary contigs
rule purge_assembly:
    input:
        polish_draft = f"{base_dir}/polishing/medaka/{{sample}}/consensus.fasta",
        overlaps = f"{base_dir}/purge/purge_haplotigs/{{sample}}/overlaps/{{sample}}.overlaps.bam",
        contig_coverage =  f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}_polished_contig_coverage.csv",
    output:
        curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
        artefacts = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.artefacts.fasta",
        log = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.contig_associations.log",
        haplotigs = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.haplotigs.fasta",
        reassignments = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.reassignments.tsv",
    threads:
        resources["purge_haplotigs"]["purge_threads"]
    log:
        'logs/{sample}/purge/purge_assembly.log'
    conda:
        "../envs/purge_haplotigs.yaml"
    shadow:
        "shallow"
    shell:
        """
        purge_haplotigs purge -g {input.polish_draft} \
            -t {threads} \
            -c {input.contig_coverage} \
            -o {wildcards.sample} \
            -d -b {input.overlaps} 2> {log} && \
        mv {wildcards.sample}* {base_dir}/purge/purge_haplotigs/{wildcards.sample}/ && \
        mv dotplots* {base_dir}/purge/purge_haplotigs/{wildcards.sample}/ 
        """

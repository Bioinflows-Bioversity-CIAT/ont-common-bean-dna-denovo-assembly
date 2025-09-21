rule build_repeat_db:
    input:
        assembly = rules.get_contamination.output.decontaminated_draft
    output:
        multiext(f"{base_dir}/annotation/repeats/RepeatModeler/{{sample}}/{{sample}}",
                 ".fa","-families.fa","-families.stk",
                 ".nhr",".nin",".njs",".nnd",".nni",
                 ".nog",".nsq","-rmod.log",".translation")
    resources:
        tmp_dir = get_big_temp
    conda:
        config["annotation"]["repeats_env"]
    threads:
        resources["annotation"]["repeats"]["BuildDatabase_threads"]
    shell:
        """
        export BLAST_USAGE_REPORT=false && \
        mkdir -p {resources.tmp_dir} && \
        cp {input.assembly} {resources.tmp_dir}{wildcards.sample}.fa && \
        cd {resources.tmp_dir} && \
        BuildDatabase -name "{wildcards.sample}" ./{wildcards.sample}.fa && \
        RepeatModeler -database {wildcards.sample} -threads {threads} && \
        mkdir -p {base_dir}/annotation/repeats/RepeatModeler/{wildcards.sample}/ && \
        mv {resources.tmp_dir}* {base_dir}/annotation/repeats/RepeatModeler/{wildcards.sample}/
        """

rule mask_repeats:
    input:
        assembly = rules.get_contamination.output.decontaminated_draft,
        repeat_db = f"{base_dir}/annotation/repeats/RepeatModeler/{{sample}}/{{sample}}-families.fa",
    output:
        multiext(f"{base_dir}/annotation/repeats/RepeatMasker/{{sample}}/{{sample}}.fa",
                 ".cat.gz",".masked",".out",".out.gff",".tbl")
    resources:
        tmp_dir = get_big_temp
    conda:
        config["annotation"]["repeats_env"]
    threads:
        resources["annotation"]["repeats"]["RepeatMasker_threads"]
    shell:
        """
        mkdir -p {resources.tmp_dir} && \
        export BLAST_USAGE_REPORT=false && \
        cp {input.assembly} {resources.tmp_dir}{wildcards.sample}.fa && \
        RepeatMasker -pa {threads} -lib {input.repeat_db} -dir {resources.tmp_dir} \
        -cutoff 255 -gff -excln -frag 40000 -xsmall {resources.tmp_dir}{wildcards.sample}.fa && \
        mkdir -p {base_dir}/annotation/repeats/RepeatMasker/{wildcards.sample}/ && \
        mv {resources.tmp_dir}* {base_dir}/annotation/repeats/RepeatMasker/{wildcards.sample}/
        """

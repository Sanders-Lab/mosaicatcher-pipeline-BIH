
rule change_id_and_sam:
    input:
        bam_orig = config["input_bam_location"]+ "{sample}/raw/{cell}.bam"
        #bam_orig = "{SM}/{ID}_sorted.bam"
    output:
        bam_out = config["input_bam_location"]+ "{sample}/all/{cell}.bam"
    conda:
        "mosaicatcher_env"
    shell:
        """
        # First, the 'ID' tag
        samtools view -H {input.bam_orig} | sed "s/ID:.*\t/ID:{wildcards.cell}\t/" | samtools reheader - {input.bam_orig} > {output.bam_out}pre1
        # next, 'SM' tag. We also have to make a new header
        samtools view -H {output.bam_out}pre1 | sed "s/SM:.*$/SM:{wildcards.sample}/" | samtools reheader - {output.bam_out}pre1 > {output.bam_out}pre2
        # then, RG:Z tag.
        samtools view -h {output.bam_out}pre2 | sed "s/RG:Z:.*/RG:Z:{wildcards.cell}/" > {output.bam_out}.sam
        # sam to bam
        samtools view -Sb {output.bam_out}.sam > {output.bam_out}
        # Remove intermediate files
        rm {output.bam_out}pre1 {output.bam_out}pre2 {output.bam_out}.sam
        """

rule add_idx:
    input:
        bam = config["input_bam_location"]+ "{sample}/all/{cell}.bam"
    output:
        bai = config["input_bam_location"]+ "{sample}/all/{cell}.bam.bai"
    conda:
        "mosaicatcher_env"
    shell:
        """
        samtools index {input.bam}
        """

rule symlink_all_to_select:
    input:
        bam = config["input_bam_location"]+ "{sample}/all/{cell}.bam",
        bai = config["input_bam_location"]+ "{sample}/all/{cell}.bam.bai"
    output:
        bam = config["input_bam_location"]+ "{sample}/selected/{cell}.bam",
        bai = config["input_bam_location"]+ "{sample}/selected/{cell}.bam.bai"
    conda:
        "mosaicatcher_env"
    shell:
        """
        ln -s ../all/{wildcards.cell}.bam {output.bam}
        ln -s ../all/{wildcards.cell}.bam.bai {output.bai}
        """

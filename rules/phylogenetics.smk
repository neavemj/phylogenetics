"""
These will build a phylogenetic tree

The inputs are:
    - assembled AIV genomes in fasta format
The outputs are:
    - a rendered phylogenetic tree in png format
"""

import os, sys

HA_types = ["H5", "H7"]


rule tmp_all:
    input:
        expand("04_phylogenetics/{sample}.tree_finished.txt", sample=config["samples"])



# made the IRMA assembly rule a checkpoint
# because we don't know what subtype will be assembled
def determine_tree(wildcards):
    # need a rule to return either a tree file (for subtypes which we have backbone sets)
    # or return a file stating that no tree will be drawn
    # the below ensures that irma_scan becomes a dependancy and is executed
    checkpoint_output = checkpoints.irma_scan.get(**wildcards).output[0]

    # now we can glob the files containing HA genes that IRMA produced
    file_names = expand("02_irma_assembly/{sample}/irma_output/A_HA_{subtype}.fasta", sample=wildcards.sample,
                  subtype=glob_wildcards("02_irma_assembly/{sample}/irma_output/A_HA_{subtype}.fasta").subtype)
    
    # get subtype as a string from the files produced
    # this relies entirely on IRMA outputting file names in the format: A_HA_H5.fasta
    # note: this determines subtype based on the first file name in the list
    # if other subtypes are assembled, would probably cause issues
    subtype = os.path.basename(file_names[0]).split("_")[2].rstrip(".fasta")
    
    # if we have a backbone for the particular subtype detected, then
    # return a tree png file - this makes the tree rules run
    if any(HA_type in subtype for HA_type in HA_types):        
        return "04_phylogenetics/{}_sequences.tree.png".format(subtype)       
    # or else, produce an empty file and the tree rules won't run
    else:
        return "04_phylogenetics/no_tree_for_this_subtype"


def aggregate_input(wildcards):
    # also need a function to just return the full IRMA HA file names
    # this is used as input to the tree building rules below
    return expand("02_irma_assembly/{sample}/irma_output/A_HA_{subtype}.fasta", sample=config["samples"],
                  subtype=glob_wildcards("02_irma_assembly/{sample}/irma_output/A_HA_{subtype}.fasta").subtype)
                  

# need an aggregate rule, which determines which tree to build
# based on the return of the determine_tree function above
rule aggregate:
    input:
        determine_tree
    output:
        "04_phylogenetics/{sample}.tree_finished.txt"
    shell:
        "ls {input} > {output}"


rule create_no_tree_file:
    output: "04_phylogenetics/no_tree_for_this_subtype"
    shell: "touch {output}"


rule tidy_fasta_names:
    message:
        """
        ** phylogenetics **
        Replacing brackets with underscores in fasta names
        IQ-Tree doesn't like unusual characters
        """
    input:
        backbone = config["program_dir"] + "introduction/raw_data/NCBI_{subtype}_backbone.fasta"
    output:
        backbone_clean = "04_phylogenetics/NCBI_{subtype}_backbone.clean.fasta"
    shell:
        """
        # use sed to replace brackets with underscores
        # however, if the bracket is the final character in the header,
        # this looks weird as an underscore
        # so, just remove the bracket if its at the line-end
        # also replace any spaces with underscores
        sed "s/)$//g" {input.backbone} | \
            sed "s/)/_/g" | \
            sed "s/(/_/g" | \
            sed "s/ /_/g" | \
            sed "s/:/_/g" | \
            sed "s/#//g" \
            > {output.backbone_clean}
        """ 


rule combine_sample_backbone:
    message:
        """
        ** phylogenetics **
        Combining assembled sequences with backbone reference sequences
        for phylogenetic analysis
        """
    input:
        dummy_file = expand("02_irma_assembly/{sample}/IRMA_COMPLETE", sample=config["samples"]),
        # the aggregate input function will return files produced by IRMA
        # could be any subtype
        HA = aggregate_input,
        backbone_clean = "04_phylogenetics/NCBI_{subtype}_backbone.clean.fasta"
    output:
        HA_sequences = "04_phylogenetics/{subtype}_sequences.fasta"
    shell:
        """
        cat \
            {input.HA} \
            {input.backbone_clean} \
            > {output.HA_sequences}
        """


rule grab_metadata:
    message:
        """
        ** phylogenetics **
        Compiling simple metadata for IQ-Tree
        """
    input:
        HA = aggregate_input,        
        backbone_clean = "04_phylogenetics/NCBI_{subtype}_backbone.clean.fasta"      
    output:
        tree_metadata = "04_phylogenetics/{subtype}_sequences.metadata.txt"
    shell:
        """
        # want to get the headers from both files and in the second column
        # indicate if the sequence is 'backbone' or the new sample
        # then can color accordingly in the tree
        # use awk to grab only rows starting with ">", then remove the ">",
        # then print the remaining text, plus add a second column of where
        # the sample came from
        
        {{ awk '/^>/ {{ gsub(">",""); print $1 "\tNEW" }}' {input.HA} ; \
          awk '/^>/ {{ gsub(">",""); print $1 "\tBACKBONE" }}' {input.backbone_clean} ; }} \
          > {output.tree_metadata}

        """
        
 
rule mafft_align:
    message:
        """
        ** phylogenetics **
        Aligning combined sequences using MAFFT
        """
    input:
        HA_sequences = "04_phylogenetics/{subtype}_sequences.fasta"
    output:
        HA_alignment = "04_phylogenetics/{subtype}_sequences.align.fasta"
    threads: 8
    shell:
        """
        mafft \
            --thread 8 \
            --adjustdirection \
            --auto \
            {input.HA_sequences} \
            > {output.HA_alignment}
        """

rule iqtree:
    message:
        """
        ** phylogenetics **
        Creating phylogenetic tree from aligned sequences
        """
    input:
        HA_alignment = "04_phylogenetics/{subtype}_sequences.align.fasta"
    output:
        HA_tree = "04_phylogenetics/{subtype}_sequences.treefile"
    params:
        model = "GTR+G"
    threads: 8
    shell:
        """
        # adding the redo flag in case snakemake wants to rerun the tree
        # otherwise, iqtree throws an error
        iqtree2 \
            -nt 8 \
            -m {params.model} \
            -s {input.HA_alignment} \
            -redo \
            -pre "04_phylogenetics/{wildcards.subtype}_sequences"
        """
        
rule draw_ggtree:
    message:
        """
        ** phylogenetics **
        Drawing the tree using ggtree in R
        """
    input:
        HA_tree = "04_phylogenetics/{subtype}_sequences.treefile",
        tree_metadata = "04_phylogenetics/{subtype}_sequences.metadata.txt"
    output:
        HA_tree_png = "04_phylogenetics/{subtype}_sequences.tree.png"
    shell:
        """
        Rscript {config[program_dir]}/phylogenetics/scripts/draw_ggtree.R \
            --tree {input.HA_tree} \
            --metadata {input.tree_metadata} \
            --output {output.HA_tree_png}
        """
        



        
        
    

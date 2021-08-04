"""
These will build a phylogenetic tree

The inputs are:
    - assembled AIV genomes in fasta format
The outputs are:
    - a rendered phylogenetic tree in png format
"""

import os, sys

configfile: "config.yaml"

rule tidy_fasta_names:
    message:
        """
        ** phylogenetics **
        Replacing brackets with underscores in fasta names
        IQ-Tree doesn't like unusual characters
        """
    input:
        backbone = "../introduction/raw_data/H7_NCBI_backbone.fasta",
    output:
        backbone_clean = "04_phylogenetics/H7_NCBI_backbone.clean.fasta"
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
            sed "s/ /_/g" \
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
        HA = "../introduction/raw_data/21-02023-03_H7.fasta",
        backbone_clean = "04_phylogenetics/H7_NCBI_backbone.clean.fasta"
    output:
        HA_sequences = "04_phylogenetics/HA-sequences.fasta"
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
        HA = "../introduction/raw_data/21-02023-03_H7.fasta",
        backbone_clean = "04_phylogenetics/H7_NCBI_backbone.clean.fasta"
    output:
        tree_metadata = "04_phylogenetics/HA-sequences.metadata.txt"
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
        HA_sequences = "04_phylogenetics/HA-sequences.fasta"
    output:
        HA_alignment = "04_phylogenetics/HA-sequences.align.fasta"
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
        HA_alignment = "04_phylogenetics/HA-sequences.align.fasta"
    output:
        HA_tree = "04_phylogenetics/HA-sequences.treefile"
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
            -pre "04_phylogenetics/HA-sequences"
        """
        
rule draw_ggtree:
    message:
        """
        ** phylogenetics **
        Drawing the tree using ggtree in R
        """
    input:
        HA_tree = "04_phylogenetics/HA-sequences.treefile",
        tree_metadata = "04_phylogenetics/HA-sequences.metadata.txt"
    output:
        HA_tree_png = "04_phylogenetics/HA-sequences.tree.pdf"
    shell:
        """
        Rscript {config[program_dir]}/phylogenetics/scripts/draw_ggtree.R
        """
        
        
        
    

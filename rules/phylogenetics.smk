"""
These will build a phylogenetic tree

The inputs are:
    - assembled AIV genomes in fasta format
The outputs are:
    - a rendered phylogenetic tree in png format
"""

import os, sys

configfile: "config.yaml"

rule combine_sample_backbone:
    message:
        """
        ** phylogenetics **
        Combining assembled sequences with backbone reference sequences
        for phylogenetic analysis
        """
    input:
        HA = "../introduction/raw_data/21-02023-03_H7.fasta",
        backbone = "../introduction/raw_data/H7_NCBI_backbone.fasta",
    output:
        HA_sequences = "04_phylogenetics/HA-sequences.fasta"
    shell:
        """
        cat \
            {input.HA} \
            {input.backbone} \
            > {output.HA_sequences}
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
        model = "GTR"
    threads: 8
    shell:
        """
        iqtree2 \
            -nt 8 \
            -m {params.model} \
            -s {input.HA_alignment} \
            -pre "04_phylogenetics/HA-sequences"
        """

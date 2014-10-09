a filtering procedure based on E-value calculated on domain architecture using MEME and MAST software v4.4. A cut-off is automatically calculated for sequences with a too low E-value with regards to the median E-value of sequences composing the cluster. The filtering procedure also removes alternatively spliced products as well as mis-annotated (e.g. truncated genes) and too highly divergent sequences. Filtered sequences are not excluded from the cluster but they are not taken into account in the next steps to avoid long-branch attraction effect and misleading ortholog identification.

Multiple alignment using MAFFT v6.240 of the full length protein-coding genes.

Masking of the multiple alignments using AL2CO to optimize the alignment for phylogenetic construction, by removing poorly informative amino acid positions. Then, the number of conserved sites after masking and the ratio of the masked columns on total columns are assessed to check the quality of the masking step.

Phylogenetic construction using PhyML v3 with 100 bootstraps.

Gene rooting and orthologous scoring using RIO implemented in a new version (v2 alpha) of the forester package that produces gene trees at the phyloXML standard. Ortholog predictions are based on the concepts (e.g. orthology, sub-tree neighbor, super-orthologs) described by Zmasek and Eddy. Bootstrapped rooted trees are used to assign an orthologs score for each sequence.

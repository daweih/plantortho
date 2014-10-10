## To Do List.
- [ ] Achieve the pipeline from GreenDB.
  - [ ] Run the whole pipeline.
    - [ ] phyml bootstrap 10 times  to test RIO. nanlandetian
  - [ ] Test the pipeline with gramene ortholog data.
    - [x] Extract proteins for test in fasta format: [parser](https://github.com/daweih/plantortho/blob/master/parser/leaf_nodes_species2pep_fa.emf.parser.pl), [data](https://github.com/daweih/plantortho/tree/master/pipeline) by daweih
    - [ ] MAFFT: multiple alignment
    - [ ] [AL2CO](https://github.com/daweih/plantortho/blob/master/bin/al2co.zip): masking of the multiple alignments to optimize the alignment for phylogenetic construction
    - [ ] PhyML: phylogenetic construction with 100 bootstraps
    - [ ] RIO: gene rooting and orthologous scoring
  - [ ] Find other database for testing.
- [ ] Run inparanoid pipeline on server.  nanlandetian
- [ ] OrthoMCL on Qomo. dotswing

## How to write in to do list
 - [Markdown on github](https://help.github.com/articles/writing-on-github/)
 - [Github flavored markdown](https://help.github.com/articles/github-flavored-markdown/)
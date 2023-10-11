#### Chimeric gene-TE file description

Column description:

- **TE**: name of the TE (Rech et al. 2022)
- **Strain TE**: AKA-017, JUT-011, MUN-016, SLA-001, TOM-007 - strain with the TE insertion
- **Tissue**: head, gut or ovary - tissue that was analized the ChIP-seq and RNA-seq 
- **Strain breakpoint**: strain without the TE insertion
- **Epigenetic effect**: 0, 1: boolean indicating if the TE has epigenetic effect or not
- **Spread**: contiguous kb up- down-stream the TE insertion showing enrichment or depletion for a histone mark
- **Average**: % of enrichment or depletion 
- **Threshold**: consistent, threshold=1, indicates if the comparisons of the epigenetic status between alleles with and without the insertion was consistent or not
- **Histone**: histone affected (H3K9me3, H3K27ac or both [bivalent])
- **Gene distance**: distance to the closer TE (retrieved from Rech et al. 2022)
- **Number of genomes present (47)**: number from 1-47 indicating in how many genomes the insertion was found
- **class**: TE class
- **order**: TE order
- **superfamily**: TE superfamily
- **family**: TE family
- **length**: TE length

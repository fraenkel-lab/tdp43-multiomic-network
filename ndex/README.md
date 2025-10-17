| Node Attribute | Definition |
|----------------|-------------|
| 00_leiden_cluster | Leiden cluster assignment |
| 01_meta_cluster | Functional cluster annotation based on GO term enrichment of Leiden clusters |
| 02_prize | Scaled prize value (for nodes included in input sets only) |
| 03_type | Type of node (based on input source type if available, Steiner nodes are specified as proteins or metabolites) |
| 04_in_prize_set | Whether node was included in input sets |
| 05_degree | Degree of node (number of edges) |
| 06_betweenness | Betweenness centrality of node |
| 07_robustness | Robustness of node based on OI randomizations, as defined in manuscript methods |
| 08_inv_specificity | Inverse specificity of node based on OI randomizations, as defined in manuscript methods |
| 09_source | Source set for input nodes (Steiner node otherwise) |
| 10_location | Subcellular location of node, derived from Gene Ontology |
| 11_general_process | General biological process of node, derived from GO |
| 12_specific_process | Specific biological process of node, derived from GO |
| 13_general_function | General molecular function of node, derived from GO |
| 14_specific_function | Specific molecular function of node, derived from GO |
| 15_log2FC_protein | Log2 fold change in iNeuron proteomics (TDP-43 KD vs. control) |
| 16_log2FC_RNA | Log2 fold change in iNeuron RNAseq (TDP-43 KD vs. control) |
| 17_combined_DE_iNeurons | Whether node has significant differential expression on RNA (gene) and/or protein level, as defined in manuscript methods |
| 18_Fig1_plot_group | Group assignment for differentially expressed proteins (DEPs), as defined in manuscript Fig. 1 |
| 19_DS_iNeurons | Whether gene has significant differential splicing as determined in iNeuron MAJIQ/RMATS SE analyses (TDP-43 KD vs. control) |
| 20_highest_magn_IncLevelDifference_rmats | Highest magnitude IncLevelDifference event detected in RMATS |
| 21_highest_magn_deltaPSI_majiq | Highest magnitude deltaPSI event detected in MAJIQ |
| 22_TDP43_protein_interactor | Whether protein has experimental evidence of physical interaction with TDP-43 (in iNeurons or from databases) |
| 23_reported_TDP43_CLIP_target | Whether gene has evidence of TDP-43 binding in exonic, intronic, or UTR regions (POSTAR3 CLIP database, neuronal tissues/cell lines) |
| 24_num_datasets_reported_misspliced | Number of total datasets/analyses where gene has reported differential or mis-splicing (inclusive of cryptic events, as defined in manuscript methods) |
| 25_num_datasets_reported_cryptic | Number of total datasets/analyses where gene has reported cryptic splicing (as defined in manuscript methods) |
| 26_num_datasets_reported_APA | Number of total datasets/analyses where gene has reported cryptic or alternative polyadenylation |
| 27_misspliced_postmortem | Whether gene has missplicing or DS detected in at least one post-mortem dataset/analysis |
| 28_cryptic_postmortem | Whether gene has cryptic splicing detected in at least one post-mortem dataset/analysis |
| 29_APA_postmortem | Whether gene has cryptic/APA detected in at least one post-mortem dataset/analysis |
| 30_num_postmortem_proteomics_datasets_DE | Number of post-mortem proteomics datasets where gene is significantly differentially expressed |
| 31_num_postmortem_RNA_datasets_DE | Number of post-mortem RNAseq datasets where gene is significantly differentially expressed |
| 32_genetic_ALS_association_OpenTargets | Quantitative genetic ALS association score from OpenTargets |

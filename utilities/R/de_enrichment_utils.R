# Functions related to diff analysis and enrichment for ineuron and postmortem
# datasets

library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(tidyr)


# small function to ensure universe is set exactly as intended
add_universe_set = function(m_go_df, universe_add) {
  new_df = as.data.frame(list(
    gs_name = 'universe_set',
    gene_symbol = universe_add
  ))
  
  for (col in setdiff(colnames(m_go_df), colnames(new_df))) {
    new_df[[col]] = 'universe'
  }
  
  return(rbind(m_go_df, new_df))
}


# This function assumes df format for enricher
get_enriched_terms = function(term2gene, limma = NULL, deseq = NULL, network_df = NULL,
                              genename_col = 'gene_name', dataset_name = 'dataset',
                              max_size = 300, de_qval_cutoff = 0.05, go_qval_cutoff = 0.05,
                              min_count = 5,  use_min_count = F,  # this only applies for the network results
                              filter_redundant_genesets = F, overlap_threshold = 0.5,
                              universe = NULL, add_universe_set = F, return_single_df = F
) {
  network_style = FALSE
  if (!is.null(limma)) {
    df_use = limma %>% mutate(
      adjpval = adj.P.Val,
      log2fc = logFC
    )
  } else if (!is.null(deseq)) {
    df_use = deseq %>% mutate(
      adjpval = padj,
      log2fc = log2FoldChange
    )
  } else if (!is.null(network_df)) {
    df_use = network_df
    network_style = TRUE
  } else {
    stop('Must pass either limma or deseq table, or network_df.')
  }
  
  if (network_style) {
    
    if (is.null(universe)) {
      stop('Must pass universe when using network_df.')
    }
    
    # If setting strict universe/background without dropping missing genes
    if (add_universe_set) {
      term2gene = add_universe_set(term2gene, universe)
    }
    
    enrich_results = list()
    # Get enrichments for each cluster individually and then finally for whole network
    for (clust in sort(unique(network_df_ann$leiden_clusters))) {
      print(clust)
      fore_genes = row.names(network_df[network_df$leiden_clusters == clust,])
      
      em2 <- clusterProfiler::enricher(fore_genes, 
                                       TERM2GENE = term2gene, 
                                       universe = universe,
                                       minGSSize = 5,
                                       maxGSSize = max_size, qvalueCutoff = go_qval_cutoff)
      
      # also include some other columns for completeness
      df1 = em2@result[c('ID', 'qvalue', 'Count', 'geneID', 'GeneRatio', 'BgRatio', 'pvalue')]
      df1[['dataset']] = dataset_name
      df1[['subset']] = paste0('Cluster_', clust)
      df1 = df1 %>% filter(qvalue < go_qval_cutoff ) %>% 
        left_join(term2gene[,c('gs_name', 'gs_subcat', 'gs_exact_source')],
                  by = c('ID' = 'gs_name')) %>% distinct()
      
      if (use_min_count) {
        df1 = df1 %>% filter(Count > min_count)
      }
      
      df1 = df1[c('ID', 
                  'dataset', 'subset', 'qvalue', 'Count', 'geneID', 
                  'GeneRatio', 'BgRatio', 'pvalue',
                  'gs_subcat', 'gs_exact_source'
      )]
      
      if (filter_redundant_genesets) {
        df1 = filterRedundantGenesets(df1, overlap_threshold = overlap_threshold)
      }
      
      enrich_results[[paste0('Cluster_', clust)]] = df1
      
    }
    # Finally do whole network enrichment
    fore_genes = row.names(network_df)
    em2 <- clusterProfiler::enricher(fore_genes, TERM2GENE = term2gene, 
                                     universe = universe,
                                     minGSSize = 5,
                                     maxGSSize = max_size, qvalueCutoff = go_qval_cutoff)
    
    # also include some other columns for completeness
    df1 = em2@result[c('ID', 'qvalue', 'Count', 'geneID', 'GeneRatio', 'BgRatio', 'pvalue')]
    df1[['dataset']] = dataset_name
    df1[['subset']] = 'Full_Network'
    df1 = df1 %>% filter(qvalue < go_qval_cutoff ) %>% 
      left_join(term2gene[,c('gs_name', 'gs_subcat', 'gs_exact_source')],
                by = c('ID' = 'gs_name')) %>% distinct()
    
    if (use_min_count) {
      df1 = df1 %>% filter(Count > min_count)
    }
    
    df1 = df1[c('ID', 
                'dataset', 'subset', 'qvalue', 'Count', 'geneID', 
                'GeneRatio', 'BgRatio', 'pvalue',
                'gs_subcat', 'gs_exact_source'
    )]
    
    if (filter_redundant_genesets) {
      df1 = filterRedundantGenesets(df1, overlap_threshold = overlap_threshold)
    }
    
    enrich_results[['Full_Network']] = df1
    
  } else {
    # THIS IS WHEN NOT USING NETWORK
    if (is.null(universe)) {
      universe = df_use[[genename_col]]
    }
    
    # If setting strict universe/background without dropping missing genes
    if (add_universe_set) {
      term2gene = add_universe_set(term2gene, universe)
    }
    
    de_neg = df_use %>% filter(adjpval < de_qval_cutoff & log2fc < 0)
    # print(dim(limma_deps_neg))
    de_pos = df_use %>% filter(adjpval < de_qval_cutoff & log2fc > 0)
    # print(dim(limma_deps_pos))
    de_all = df_use %>% filter(adjpval < de_qval_cutoff)
    
    fore_sets = list(de_neg[[genename_col]],
                     de_pos[[genename_col]],
                     de_all[[genename_col]]
    )
    # print(str(fore_sets))
    fore_sets_names = c('negative', 'positive', 'all')
    enrich_results = list()
    
    for (j in seq_along(fore_sets)) {
      fore_genes=fore_sets[[j]]
      em2 <- clusterProfiler::enricher(fore_genes, TERM2GENE = term2gene, 
                                       universe = universe,
                                       minGSSize = 5, pvalueCutoff = 1,
                                       maxGSSize = max_size, qvalueCutoff = 1)
      
      # also include some other columns for completeness
      df1 = em2@result[c('ID', 'qvalue', 'Count', 'geneID', 'GeneRatio', 'BgRatio', 'pvalue')]
      df1[['dataset']] = dataset_name
      df1[['subset']] = fore_sets_names[j]
      df1 = df1 %>% filter(qvalue < go_qval_cutoff ) %>% 
        left_join(term2gene[,c('gs_name', 'gs_subcat', 'gs_exact_source')],
                  by = c('ID' = 'gs_name')) %>% distinct()
      
      df1 = df1[c('ID', 
                  'dataset', 'subset', 'qvalue', 'Count', 'geneID', 
                  'GeneRatio', 'BgRatio', 'pvalue',
                  'gs_subcat', 'gs_exact_source'
      )]
      
      if (filter_redundant_genesets) {
        df1 = filterRedundantGenesets(df1, overlap_threshold = overlap_threshold)
      }
      
      enrich_results[[fore_sets_names[j]]] = df1
    }
    
  }
  if (return_single_df) {
    return(Reduce(rbind, enrich_results))
  }
  return(enrich_results)
}
###################################################################################

# similar to function above, but runs GSEA instead of GO
# right now only does GSEA on genes sorted by pval*log2FC
get_enriched_gsea <- function(term2gene,
                              limma = NULL,
                              deseq = NULL,
                              genename_col = 'gene_name',
                              dataset_name = 'dataset',
                              max_size = 300,
                              de_qval_cutoff = 0.05,
                              go_qval_cutoff = 0.05,
                              min_count = 5,
                              use_min_count = FALSE,
                              filter_redundant_genesets = FALSE,
                              overlap_threshold = 0.5,
                              universe = NULL,
                              return_single_df = FALSE) {
  
  network_style <- FALSE
  
  if (!is.null(limma)) {
    df_use <- limma %>% mutate(
      pval    = P.Value,
      adjpval = adj.P.Val,
      log2fc  = logFC
    )
  } else if (!is.null(deseq)) {
    df_use <- deseq %>% mutate(
      pval    = pvalue,
      adjpval = padj,
      log2fc  = log2FoldChange
    )
  } else if (exists("network_df")) {          # keeps previous behaviour
    df_use <- network_df
    network_style <- TRUE
  } else {
    stop("Must pass either limma, deseq, or a network_df object.")
  }
  
  ## -------------------------------------------------------------------- ##
  ##  Zero-p-value fixer: replace pval == 0 with a tiny value slightly
  ##  smaller (more significant) than the next-smallest non-zero p-value  ##
  ## -------------------------------------------------------------------- ##
  if (any(df_use$pval == 0, na.rm = TRUE)) {
    
    smallest_nz <- min(df_use$pval[df_use$pval > 0], na.rm = TRUE)
    
    # Fallback in the (unlikely) event *all* pvals were zero
    if (is.infinite(smallest_nz)) smallest_nz <- 1e-300
    
    new_p <- smallest_nz * 0.99         # “a little smaller”
    
    df_use <- df_use %>%
      mutate(pval = ifelse(pval == 0, new_p, pval))
  }
  
  ## -------------------------------------------------------------------- ##
  ##  downstream steps remain unchanged                                   ##
  ## -------------------------------------------------------------------- ##
  df_use <- df_use %>%
    mutate(signed_pval = -log10(pval) * log2fc)
  
  sorted_gene_list <- df_use$signed_pval
  names(sorted_gene_list) <- df_use[[genename_col]]
  sorted_gene_list <- na.omit(sorted_gene_list)
  sorted_gene_list <- sort(sorted_gene_list, decreasing = TRUE)
  
  em2 <- clusterProfiler::GSEA(sorted_gene_list,
                               TERM2GENE     = term2gene,
                               eps           = 1e-50,
                               maxGSSize     = max_size,
                               pvalueCutoff  = 1,
                               seed          = TRUE,
                               nPermSimple   = 10000)
  
  df1 <- em2@result %>%
    select(ID, qvalue, setSize, core_enrichment,
           NES, enrichmentScore, pvalue) %>%
    mutate(dataset = dataset_name) %>%
    filter(qvalue < go_qval_cutoff) %>%
    left_join(term2gene[, c("gs_name", "gs_subcat", "gs_exact_source")],
              by = c("ID" = "gs_name")) %>%
    distinct() %>%
    mutate(geneID = core_enrichment) %>%
    select(ID, dataset, qvalue, geneID,
           NES, setSize, enrichmentScore, pvalue,
           gs_subcat, gs_exact_source)
  
  if (filter_redundant_genesets) {
    df1 <- filterRedundantGenesets(df1,
                                   overlap_threshold = overlap_threshold)
  }
  
  return(df1)
}




###############################################################################

compile_DE_results = function(results_list, names, genecols, modalities = NULL, 
                              drop_missing_genes = T, drop_duplicated_genes = T) {
  if (is.null(modalities)) {
    # Determine modalities through the column names of results -- auto version
    modalities = sapply(results_list, function(x) {
      if ('adj.P.Val' %in% colnames(x)) {
        # limma specific column
        'prot'
      } else if ('padj' %in% colnames(x)) {
        # deseq2 specific column
        'RNA'
      } else {
        # deal with other modalities that might also use deseq2 results format
        'unspecified'
      }
    })
  }
  
  for (i in seq_along(results_list)) {
    print(names[i])
    results_df = results_list[[i]]
    if (modalities[i] == 'prot') {
      results_df = results_df[,c('logFC', 't', 'P.Value', 'adj.P.Val', genecols[i])]
    } else if (modalities[i] == 'RNA') {
      results_df = results_df[,c('log2FoldChange', 'stat', 'pvalue', 'padj', genecols[i])]
    } else {
      stop('No other modalities currently supported.')
    }
    
    # Rename the columns 
    colnames(results_df) = c(paste0(c('log2fc', 't', 'pval', 'adjpval'), '__', modalities[i], '__', names[i]), 'external_gene_name') 
    
    # Drop missing genes (by name)
    # TODO: Consider other forms of missingness
    if (drop_missing_genes) {
      results_df = results_df[results_df$external_gene_name != '',]
      results_df = results_df[results_df$external_gene_name != 0,]
      results_df = results_df[!is.na(results_df$external_gene_name),]
    }
    
    # Keep only the most significant gene change (when multiple)
    if (drop_duplicated_genes) {
      print('Dropping duplicated genes; keeping smallest adj. p-values')
      adjp_col = paste0( 'adjpval', '__', modalities[i], '__', names[i])
      print(paste0('Original nrows: ', nrow(results_df)))
      results_df = results_df %>%
        group_by(external_gene_name) %>% 
        slice_min(!!sym(adjp_col)) %>%
        ungroup()
      
      print(paste0('New nrows: ', nrow(results_df)))
    }
    
    # Merge if not first 
    if (i == 1) {
      results_df_full = results_df
      
    } else {
      results_df_full = merge(results_df_full, results_df, by = 'external_gene_name', all=T)
    }
  }
  
  return(results_df_full)
}






############################################################################3


######################################################################################
# Also try version of this from new package called simplifyEnrichment 
# https://github.com/jokergoo/simplifyEnrichment

library(simplifyEnrichment)
# mat = GO_similarity(go2, ont = 'MF', measure = "Wang")
# mat = GO_similarity(go2, ont = 'BP', measure = "Wang")


filterRedundantGenesets <- function(df, overlap_threshold = 0.5, gene_delim = "/") {
  # Ensure the required columns are present
  if (!all(c("geneID", "qvalue") %in% colnames(df))) {
    stop("The dataframe must contain 'geneID' and 'qvalue' columns.")
  }
  
  # Sort the dataframe by qvalue (lowest first)
  df_sorted <- df[order(df$qvalue), ]
  
  # Initialize an index vector to keep track of selected genesets
  selected_indices <- c()
  
  # Helper function to get unique genes from the geneID string
  getGenes <- function(genes_str) {
    unique(unlist(strsplit(genes_str, split = gene_delim)))
  }
  
  # Pre-calculate gene lists for all rows to avoid repeated splitting
  geneLists <- lapply(df_sorted$geneID, getGenes)
  
  # Iterate over the sorted dataframe rows
  for (i in seq_len(nrow(df_sorted))) {
    currentGenes <- geneLists[[i]]
    keep <- TRUE
    
    # Compare with each already selected geneset
    for (j in selected_indices) {
      otherGenes <- geneLists[[j]]
      # Calculate Jaccard similarity: intersection / union
      inter <- length(intersect(currentGenes, otherGenes))
      uni <- length(union(currentGenes, otherGenes))
      jaccard <- if (uni > 0) inter / uni else 0
      
      # If the similarity exceeds the threshold, mark this geneset to be dropped
      if (jaccard >= overlap_threshold) {
        keep <- FALSE
        break
      }
    }
    
    # If it is not redundant, add its index to the selected list
    if (keep) {
      selected_indices <- c(selected_indices, i)
    }
  }
  
  # Return the filtered dataframe
  return(df_sorted[selected_indices, ])
}

# test = filterRedundantGenesets(leiden_enrichments_full[['12']],overlap_threshold = 0.4)

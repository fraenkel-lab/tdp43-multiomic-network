# Plotting related utilities for iNeuron multiomic figures

library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rlang)
library(forcats)
library(ComplexHeatmap)

# Volcano plotting for DE results

volcano_limma <- function(
    results_table,
    prot_trans_table,
    using_uniprot      = TRUE,
    names_from_column  = NULL,
    ylim               = NULL,
    center_y = TRUE,   # center around y axis (symm. x lims)
    ## aesthetics
    text_size          = 20,
    label_size         = 4,
    max_overlaps       = 20,
    legend_position    = 'none',
    ## thresholds
    adjp_threshold     = 0.10,
    abslogfc_threshold = 0,
    # Black color for selected points
    color_black = NULL,    # pass in vector of gene names -- only colors significant ones
    always_label       = NULL,          # character vector of gene/protein names
    shadow_colour      = "grey",       # outline colour for shadow text
    always_nudge_x  = 0.1,   # min segment length to push away from other labels
    always_nudge_y = 0.1,
    always_label_box = TRUE,     # make F to remove boxes from always label points
    always_label_color = 'black',
    ## outputs
    return_plot        = FALSE,
    return_table       = FALSE
) {
  
  
  res <- results_table %>%
    mutate(
      neglog_qval = -log10(adj.P.Val),
      sig         = adj.P.Val < adjp_threshold & abs(logFC) > abslogfc_threshold,
      sign_adjp_only = adj.P.Val < adjp_threshold,
      sign_label  = case_when(
        sig &  logFC >  abslogfc_threshold ~ "POS",
        sig &  logFC < -abslogfc_threshold ~ "NEG",
        TRUE                                ~ "NS"
      )
    )
  
  
  ## -------------------------------------------------------------------- ##
  ## 2. Attach readable names                                              ##
  ## -------------------------------------------------------------------- ##
  if (!is.null(names_from_column)) {
    res$name   <- res[[names_from_column]]
    using_uniprot <- FALSE
    
    # still need to set uniprot_id column
    # TODO: Throw error if no uniprot column present at all
    # right now assuming row names 
    res$uniprot_id = row.names(res)
  } else {
    res$name <- rownames(res)
  }
  
  if (using_uniprot) {
    ids <- sub("([_\\|]).*$", "", rownames(res))     # keep first token before “_” or “|”
    names(ids) = row.names(res)
    res <- res %>%
      mutate(
        uniprot_id = ids,
        name       = translate_uniprot(ids, prot_trans_table)
      )
  }
  
  if (!is.null(color_black)) {
    res = res %>% mutate(
      sign_label  = case_when(
        sign_label %in% c('POS', 'NEG') & name %in% color_black ~ "LABEL",
        TRUE ~ sign_label
      )
    )
  }
  
  ## -------------------------------------------------------------------- ##
  ## 3. Base volcano plot                                                  ##
  ## -------------------------------------------------------------------- ##
  g <- ggplot(res, aes(logFC, neglog_qval, colour = sign_label)) +
    geom_point() +
    geom_hline(yintercept = c(-log10(adjp_threshold), 0), linetype = c(2, 1), colour = "grey30") +
    geom_vline(xintercept = c(-abslogfc_threshold, 0, abslogfc_threshold),
               linetype   = c(2, 1, 2), colour = "grey30") +
    scale_colour_manual(values = c(NS = "grey", POS = "red", NEG = "blue", LABEL = 'black')) +
    labs(x = "Log2 Fold Change",
         y = "-log10(adj. p-value)",
         colour = NULL) +
    theme_bw(base_size = text_size) +
    theme(legend.position = legend_position)
  
  if (!is.null(ylim)) g <- g + ylim(ylim)
  
  if (center_y) {
    xlim_max = max(max(res$logFC), abs(min(res$logFC))) + 0.1
    g = g + xlim(c(-xlim_max, xlim_max))
  }
  
  ## -------------------------------------------------------------------- ##
  ## 4. Repelled labels for significant hits (minus always_label genes)   ##
  ## -------------------------------------------------------------------- ##
  label_data <- res %>%
    filter(sig) %>%
    { if (!is.null(always_label)) filter(., !name %in% always_label) else . }
  
  g <- g +
    ggrepel::geom_text_repel(
      data              = label_data,
      aes(label         = name),
      size              = label_size,
      box.padding       = 0.5,
      max.overlaps      = max_overlaps,
      min.segment.length= 0,
      show.legend       = FALSE
    )
  
  ## -------------------------------------------------------------------- ##
  ## 5. Always-present shadow labels (no overlap filtering)               ##
  ## -------------------------------------------------------------------- ##
  if (!is.null(always_label)) {
    keep <- filter(res, name %in% always_label)
    
    if (nrow(keep)) {
      if (always_label_box) {
        g <- g +
          geom_label_repel(
            data        = keep,
            aes(label   = name,
                segment.color = 'black'),
            size        = label_size,
            # bg.colour   = shadow_colour,
            fill = 'white',
            min.segment.length = 0,
            nudge_x  = always_nudge_x,
            nudge_y = always_nudge_y,
            # colour      = "black",
            show.legend = FALSE
          )
      } else {
        g <- g +
          geom_text_repel(
            data        = keep,
            aes(label   = name,
                segment.color = 'black'),
            size        = label_size,
            min.segment.length = 0,
            nudge_x  = always_nudge_x,
            nudge_y = always_nudge_y,
            colour      = "black",
            show.legend = FALSE
          )
      }
      
    }
  }
  
  if (return_table) return(res)
  if (return_plot)  return(g)
  
  print(g)
}

# Helper function for getting protein names
# prot_trans_table must be df returned by biomart with appropriate columns
translate_uniprot = function(uniprot_list, prot_trans_table,
                              key_col = 'uniprot_gn_id', out_col = 'external_gene_name') {
  tbl = prot_trans_table[!duplicated(prot_trans_table[[key_col]]), , drop = FALSE]
  m = match(uniprot_list, tbl[[key_col]])
  res = ifelse(is.na(m), uniprot_list, tbl[[out_col]][m])
  miss = setdiff(uniprot_list, tbl[[key_col]])
  if (length(miss)) warning(paste0(length(miss), ' uniprot IDs not found'))
  if (!is.null(names(uniprot_list))) res[is.na(m)] = names(uniprot_list)[is.na(m)]
  names(res) = uniprot_list
  res[res == ''] = names(res)[res == '']
  res
}


# Panel padding so every facet contains the same set of pathway IDs
# ids_union must have columns: associated_cluster2, ID
pad_panel = function(ids_union, df) {
  ids_union %>%
    left_join(df, by = c("associated_cluster2", "ID"))
}

# Build the faceted enrichment dot plot
base_dotplot = function(dat, fill_aes, highlight_splice){
  
  ggplot(dat, aes(xlabels, ID, size = Count,
                  fill = {{fill_aes}}, color = {{highlight_splice}})) +
    geom_point(shape = 21, stroke = 1, na.rm = TRUE) +
    facet_grid(rows  = vars(associated_cluster2),
               cols  = vars(dataset_type),
               scales = "free", space = "free", switch = "y") +
    scale_color_manual(values = c("TRUE"  = "black",
                                  "FALSE" = "transparent"), guide = "none") +
    theme_bw() +
    theme(axis.text.x       = element_text(angle = 45, hjust = 1),
          axis.title        = element_blank(),
          strip.placement   = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background  = element_rect(fill = "grey90", colour = NA))
}

# Function to generate horizontal heatmap with ComplexHeatmap

plot_multiomic_heatmap_genes_x = function(
    df,
    genes,
    gene_groups = NULL,
    gene_group_col = NULL,
    gene_order = NULL,
    col_anno_cols = NULL,
    custom_anno_colors = list(),
    modality_colours = c(RNA = "#377EB8", Proteomics = "#E41A1C"),
    star_levels = c("***"=0.001,"**"=0.01,"*"=0.05),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    gene_label_angle = 45,
    legends_below = F,
    ...
) {
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  
  # Identify log2FC columns
  log2_cols = grep("^log2fc__", names(df), value=TRUE)
  meta = tibble(col=log2_cols) %>%
    separate(col, into=c("metric","modality","dataset"), sep="__", extra="merge", remove=FALSE) %>%
    mutate(
      modality=ifelse(modality=="prot","Proteomics","RNA"),
      group=case_when(
        str_detect(dataset, regex("iNeuron", ignore_case=TRUE)) ~ "iNeuron",
        modality=="RNA" ~ "RNAseq",
        TRUE ~ "Proteomics"
      ),
      adj_col=paste0("adjpval__", ifelse(modality=="Proteomics","prot","RNA"), "__",dataset),
      original_order=row_number()
    ) %>%
    arrange(
      factor(group, levels=c("iNeuron","RNAseq","Proteomics")),
      factor(modality, levels=c("RNA","Proteomics")),
      original_order
    )
  log2_cols = meta$col
  
  df_sub = df %>%
    filter(.data$external_gene_name %in% genes) %>%
    distinct(external_gene_name, .keep_all=TRUE)
  
  mat = df_sub %>%
    dplyr::select(external_gene_name, all_of(log2_cols)) %>%
    column_to_rownames("external_gene_name") %>%
    as.matrix() %>%
    t()
  
  pmat = df_sub %>%
    dplyr::select(all_of(meta$adj_col)) %>%
    as.matrix() %>%
    t()
  
  # Add clean names 
  meta = meta %>%
    mutate(
      clean_dataset = str_remove(dataset, "_RNAseq|_proteomics")
    )
  rownames(mat) = meta$clean_dataset
  rownames(pmat) = meta$clean_dataset
  
  star_mat = matrix("", nrow=nrow(pmat), ncol=ncol(pmat), dimnames=dimnames(pmat))
  for (sym in names(star_levels)) {
    star_mat[pmat < star_levels[[sym]]] = sym
  }
  
  rng = max(abs(mat), na.rm=TRUE)
  if (!is.finite(rng)) stop("log2FC matrix empty or invalid.")
  col_fun = colorRamp2(c(-rng,0,rng), c("blue","white","red"))
  
  # Column annotation (genes)
  col_ann = NULL
  if (!is.null(col_anno_cols)) {
    ca_df = df_sub %>%
      dplyr::select(external_gene_name, all_of(col_anno_cols)) %>%
      column_to_rownames("external_gene_name")
    anno_colors = list()
    for (col in col_anno_cols) {
      if (col %in% names(custom_anno_colors)) {
        anno_colors[[col]] = custom_anno_colors[[col]]
      } else {
        vals = unique(ca_df[[col]])
        n_vals = length(vals)
        colors = if (n_vals <=8) RColorBrewer::brewer.pal(n_vals,"Set2") else grDevices::rainbow(n_vals)
        anno_colors[[col]] = setNames(colors, vals)
      }
    }
    col_ann = HeatmapAnnotation(df=ca_df, col=anno_colors, show_annotation_name=TRUE)
  }
  
  # Column split (genes)
  col_split = NULL
  if (!is.null(gene_group_col)) {
    col_split = factor(df_sub[[gene_group_col]])
  } else if (!is.null(gene_groups)) {
    membership = unlist(lapply(names(gene_groups),
                                \(grp) setNames(rep(grp,length(gene_groups[[grp]])), gene_groups[[grp]])))
    col_split = factor(membership[colnames(mat)], levels=names(gene_groups))
  } else if (!is.null(gene_order)) {
    mat = mat[,gene_order, drop=FALSE]
  }
  
  # Row split (datasets)
  row_split = factor(meta$group, levels=c("iNeuron","RNAseq","Proteomics"))
  row_ann = rowAnnotation(
    modality=meta$modality,
    col=list(modality=modality_colours),
    show_annotation_name=FALSE
  )
  
  ht = Heatmap(
    mat,
    name="log2FC",
    col=col_fun,
    cluster_rows=cluster_rows,
    cluster_columns=cluster_columns,
    show_column_names=show_column_names,
    column_names_rot=gene_label_angle,
    row_split=row_split,
    column_split=col_split,
    rect_gp = gpar(col = "grey30", lwd = 0.5),  # add border to all cells
    top_annotation=col_ann,
    left_annotation=row_ann,
    cell_fun=function(j,i,x,y,w,h,fill){
      label = star_mat[i,j]
      if(label!="") grid.text(label,x,y,gp=gpar(fontsize=9))
    },
    ...
  )
  if (legends_below) {
    draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = TRUE,       # merge all legends in one area
      # legend_direction = "horizontal"
    )
  } else {
    draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  }
  invisible(ht)
}
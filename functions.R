plot_rich_reads_samlenames_lm <- function(physeq, group = "Site", label = "Repeat", delimiter="_"){
  rish <- estimate_richness(physeq, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(physeq))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeat"] <-unlist(purrr::map(stringr::str_split(rownames(physeq@sam_data), delimiter, 2), function(x) x[[2]]))
  reads.summary["Site"] <- physeq@sam_data[[group]]
  library(ggrepel)
  require(ggforce)
  p1 <- ggplot(data=reads.summary) + 
    geom_point(aes(y=otus, x=log2(reads), color=Site),size=3) + 
    geom_text_repel(aes(y=otus, x=log2(reads), label=paste0(Repeat))) + 
    theme_bw() +
    geom_smooth(aes(y=otus, x=log2(reads), fill=Site, color=Site),method=lm, se=FALSE, ymin = 1) + 
    scale_x_continuous(sec.axis = sec_axis(sec.axis ~ 2**.)) 

  return(p1)
}


beta_custom_norm_NMDS_elli_w <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  ps@otu_table[ps@otu_table < 0] <- 0
  ordination.b <- ordinate(ps, "NMDS", "bray")
  mds <- as.data.frame(ordination.b$points)
  p  <-  plot_ordination(ps,
                         ordination.b,
                         type="sample",
                         color = Color,
                         title="NMDS - Bray-Curtis",
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    annotate("text",
    x=min(mds$MDS1) + abs(min(mds$MDS1))/4,
    y=max(mds$MDS2),
    label=paste0("Stress -- ", round(ordination.b$stress, 3))) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") +
    scale_colour_viridis_d(option = "magma", 
                       aesthetics = "color", 
                       begin = 0, 
                       end = 0.8)
  
  return(p)
}

plot_alpha_w_toc <- function(ps, group, metric) {
  
  require(phyloseq)
  require(ggplot2)
  
  ps_a <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  er <- estimate_richness(ps_a)
  df_er <- cbind(ps_a@sam_data, er)
  df_er <- df_er %>% dplyr::select(c(group, metric))
  stat.test <- aov(as.formula(paste0(metric, "~", group)), data = df_er) %>%
    rstatix::tukey_hsd()
  y <-  seq(max(er[[metric]]), length=length(stat.test$p.adj), by=max(er[[metric]]/20))

  plot_richness(ps_a, x=group, measures=metric) + 
    geom_boxplot() +
    geom_point(size=1.2, alpha=0.3) +
    ggpubr::stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", 
      y.position = y) +
    theme_light() + 
    scale_color_brewer(palette="Dark2") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank()) +
    labs(y=paste(metric, "index")) 
}

plot_alpha_w_toc_2 <- function(ps, group, metric) {
  
  require(phyloseq)
  require(ggplot2)
  
  ps_a <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  er <- estimate_richness(ps_a)
  df_er <- cbind(ps_a@sam_data, er)
  df_er <- df_er %>% dplyr::select(c(group, metric))
  stat.test <- aov(as.formula(paste0(metric, "~", group)), data = df_er) %>%
    rstatix::tukey_hsd()
  y <- seq(max(er[[metric]]), length=length(stat.test$p.adj.signif[stat.test$p.adj.signif != "ns"]), by=max(er[[metric]]/20))
  
  plot_richness(ps_a, x=group, measures=metric) + 
    geom_boxplot() +
    geom_point(size=1.2, alpha=0.3) +
    ggpubr::stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", 
      y.position = y,
      hide.ns=TRUE) +
    theme_light() + 
    scale_color_brewer(palette="Dark2") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank()) +
    labs(y=paste(metric, "index")) 
}

phyloseq_to_ampvis2 <- function(physeq) {
  #check object for class
  if(!any(class(physeq) %in% "phyloseq"))
    stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
  
  #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
  if(is.null(physeq@tax_table))
    stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
  
  #OTUs must be in rows, not columns
  if(phyloseq::taxa_are_rows(physeq))
    abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
  else
    abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
  
  #tax_table is assumed to have OTUs in rows too
  tax <- phyloseq::tax_table(physeq)@.Data
  
  #merge by rownames (OTUs)
  otutable <- merge(
    abund,
    tax,
    by = 0,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  colnames(otutable)[1] <- "OTU"
  
  #extract sample_data (metadata)
  if(!is.null(physeq@sam_data)) {
    metadata <- data.frame(
      phyloseq::sample_data(physeq),
      row.names = phyloseq::sample_names(physeq), 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    #check if any columns match exactly with rownames
    #if none matched assume row names are sample identifiers
    samplesCol <- unlist(lapply(metadata, function(x) {
      identical(x, rownames(metadata))}))
    
    if(any(samplesCol)) {
      #error if a column matched and it's not the first
      if(!samplesCol[[1]])
        stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
    } else {
      #assume rownames are sample identifiers, merge at the end with name "SampleID"
      if(any(colnames(metadata) %in% "SampleID"))
        stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
      metadata$SampleID <- rownames(metadata)
      
      #reorder columns so SampleID is the first
      metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
    }
  } else
    metadata <- NULL
  
  #extract phylogenetic tree, assumed to be of class "phylo"
  if(!is.null(physeq@phy_tree)) {
    tree <- phyloseq::phy_tree(physeq)
  } else
    tree <- NULL
  
  #extract OTU DNA sequences, assumed to be of class "XStringSet"
  if(!is.null(physeq@refseq)) {
    #convert XStringSet to DNAbin using a temporary file (easiest)
    fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
    Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
  } else
    fastaTempFile <- NULL
  
  #load as normally with amp_load
  ampvis2::amp_load(
    otutable = otutable,
    metadata = metadata,
    tree = tree,
    fasta = fastaTempFile
  )
}

detachAllPackages <- function() {

  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

  package.list <- setdiff(package.list,basic.packages)

  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE, force = TRUE)

}

ps_vst <- function(ps, factor){
  diagdds = phyloseq::phyloseq_to_deseq2(ps, as.formula(paste( "~", factor)))              
  diagdds = DESeq2::estimateSizeFactors(diagdds, type="poscounts")
  diagdds = DESeq2::estimateDispersions(diagdds, fitType = "local") 
  pst <- DESeq2::varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(SummarizedExperiment::assay(pst))) 
  # pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  return(ps.varstab)
} 

filter_chrmitnas <- function(physeq){
  ps_object <- physeq 
  ps_object <- subset_taxa(ps_object, Phylum != "NA")
  
  ps_object@tax_table[is.na(ps_object@tax_table)] <- ps@tax_table@.Data
  ps_object <- subset_taxa(ps_object,
                           !(Family  == "Mitochondria" |
                               Class   == "Chloroplast" |
                               Order   == "Chloroplast"))
  ps_object@tax_table <- dplyr::na_if(ps_object@tax_table, TRUE)
  ps.f <- ps_object
  
  out_old <- capture.output(physeq)[4] %>% 
    str_extract(regex("([1-9]).*(?= by)"))
  out_new <- capture.output(ps.f)[4] %>% 
    str_extract(regex("([1-9]).*(?= by)"))
  message(paste0("old = ", out_old))
  message(paste0("filtered = ", out_new))
  
  return(ps.f)
  
}

beta_custom_norm_PCA_elli_w <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  ps@otu_table[ps@otu_table < 0] <- 0
  ordination.b <- ordinate(ps, "PCoA", "euclidean")
  mds <- as.data.frame(ordination.b$points)
  p  <-  plot_ordination(ps,
                         ordination.b,
                         type="sample",
                         color = Color,
                         title="PCA - clr-transformed",
                         # title=NULL,
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") +
    scale_colour_viridis_d(option = "magma", 
                       aesthetics = "color", 
                       begin = 0, 
                       end = 0.8)
  
  return(p)
}


filter_chrmitnas <- function(physeq){
  ps_object <- physeq 
  ps_object <- subset_taxa(ps_object, Phylum != "NA")
  
  ps_object@tax_table <- ps_object@tax_table %>% 
    as.data.frame() %>% 
    mutate_all(~ replace_na(., "unknown")) %>% 
    as.matrix() %>% 
    tax_table()
  
  ps_object <- subset_taxa(ps_object,
                           !(Family  == "Mitochondria" | 
                               Class   == "Chloroplast" | 
                               Order   == "Chloroplast"))
  
  ps_object@tax_table <- ps_object@tax_table %>% 
    as.data.frame() %>% 
    mutate_all(~ na_if(., "unknown")) %>% 
    as.matrix() %>% 
    tax_table()
  
  ps.f <- ps_object
  
  out_old <- capture.output(physeq)[4] %>%
    str_extract(regex("([1-9]).*(?= by)"))
  out_new <- capture.output(ps.f)[4] %>%
    str_extract(regex("([1-9]).*(?= by)"))
  message(paste0("old = ", out_old))
  message(paste0("filtered = ", out_new))
  
  return(ps.f)
  
}

color_filt <- function(ps, df){
  
  library(tidyverse)
  library(reshape2)
  library(gridExtra)
  
  l = list()
  for (i in levels(df$module)){
    message(i)
    color_name <-  df %>% 
      filter(module == i) %>% 
      pull(asv) %>% 
      unique()
    ps.col <- prune_taxa(color_name, ps)
    amp.col <- phyloseq_to_ampvis2(ps.col)
    heat <- amp_heatmap(amp.col, tax_show = 60, 
                        group_by = "Day", 
                        tax_aggregate = "OTU",
                        tax_add = "Genus", 
                        normalise=FALSE, 
                        showRemainingTaxa = TRUE)
    
    ps.rel  <-  phyloseq::transform_sample_counts(ps, function(x) x / sum(x) * 100)
    ps.rel.col <- prune_taxa(color_name, ps.rel)
    amp.r <- phyloseq_to_ampvis2(ps.rel.col)
    heat.rel <- amp_heatmap(amp.r, tax_show = 60, 
                            group_by = "Day", 
                            tax_aggregate = "OTU",
                            tax_add = "Genus", 
                            normalise=FALSE, 
                            showRemainingTaxa = TRUE)
    
    tree <- ps.col@phy_tree
    taxa <- as.data.frame(ps.col@tax_table@.Data)
    p1 <- ggtree(tree) + 
      geom_tiplab(size=2, align=TRUE, linesize=.5) +
      theme_tree2()
    
    taxa[taxa == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium"
    taxa[taxa == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
    
    tx <- taxa %>% 
      rownames_to_column("id") %>% 
      mutate(id = factor(id, levels = rev(get_taxa_name(p1)))) %>% 
      dplyr::select(-c(Kingdom, Species, Order)) %>% 
      melt(id.var = 'id')
    
    p2 <- ggplot(tx, aes(variable, id)) + 
      geom_tile(aes(fill = value), alpha = 0.4) +
      geom_text(aes(label = value), size = 3) +
      theme_bw() + 
      theme(legend.position = "none",
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) 
    
    p <- ggpubr::ggarrange(p1, p2, widths = c(0.9, 1))
    l[[i]] <- list("ps" = ps.col, 
                   "amp" = amp.col,
                   "heat" = heat,
                   "heat_rel" = heat.rel,
                   "tree" = p,
                   "taxa" = knitr::kable(taxa) %>%
                     kableExtra::kable_paper() %>%
                     kableExtra::scroll_box(width ="100%", height = "200px"))
    
  }
  return(l)
}


get_list_from_ps <- function(physeq, amaz_fac) {
  my_keys <- levels(sample_data(physeq)[[amaz_fac]])
  mylist <- vector(mode="list", length=length(my_keys))
  names(mylist) <- my_keys
  for (i in levels(sample_data(physeq)[[amaz_fac]])) { 
    x <- prune_samples(sample_data(physeq)[[amaz_fac]] %in% i, physeq)
    x <- prune_taxa(taxa_sums(x) > 0, x)
    mylist[[i]] <- taxa_names(x)
  }
  return(mylist)
}


beta_custom_norm_PCoA_elli_w <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  ps@otu_table[ps@otu_table < 0] <- 0
  ordination.b <- ordinate(ps, "PCoA", "bray")
  mds <- as.data.frame(ordination.b$points)
  p  <-  plot_ordination(ps,
                         ordination.b,
                         type="sample",
                         color = Color,
                         title="PCoA - bray - rarefied",
                         # title=NULL,
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    theme(legend.position = "none") +
    scale_colour_viridis_d(option = "magma", 
                       aesthetics = "color", 
                       begin = 0, 
                       end = 0.8)
  
  return(p)
}


bargraph <- function(ps, rank, threshold=0.05){
  require(dplyr)
  require(ggplot2)
  require(phyloseq)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps2 <- tax_glom(ps, taxrank = rank)
  ps3 = transform_sample_counts(ps2, function(x) x / sum(x) )
  data <- psmelt(ps3) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) # convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- data %>% group_by(Plot) %>% mutate(median=median(data$Abundance))
  remainder <- medians[medians$median <= threshold,]$Plot
  
  # create palette long enough for our data
  base.palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", 
                    "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", 
                    "darksalmon", "dodgerblue3", "steelblue1", "darkgoldenrod1", "brown1", "cyan1", "darkgrey")
  required.colors <- nlevels(factor(data$Plot))
  repeats = required.colors %/% length(base.palette) + 1
  palette <- rep(base.palette, length.out = repeats * length(base.palette))
  
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = palette) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle('ITS2 - Phylum level')
  
}

norm_anc <- function(ps, anc){
  samp_frac <-  anc$samp_frac
  samp_frac[is.na(samp_frac)] <-  0 
  log_obs_abn <-  log(anc$feature_table + 1)
  log_corr_abn <-  t(t(log_obs_abn) - samp_frac)
  ps.out <- ps
  otu_table(ps.out) <- otu_table(t(log_corr_abn), taxa_are_rows = FALSE)
  return(ps.out)
}

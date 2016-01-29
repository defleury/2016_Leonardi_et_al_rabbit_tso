#!/usr/bin/Rscript
################################################################################
#Project: Irina, rabbit gut microbiota
#Full R analysis script
#
#(i)		Load & process raw data (OTU tables, taxonomy, tree, sample metadata)
#(ii)		Get alpha-div per sample
#(iii)	Get beta-div between all samples
#(iv)		Get taxonomies per sample
#(v)		Get significantly up- and down-regulated OTUs (EdgeR & phyloseq)
#(vi)		Compare different sequenced subregions
#
#
#2015-05-01
#sebastian.schmidt@imls.uzh.ch
################################################################################

################################################################################
################################################################################
#Load packages
library("bigmemory");
library("biganalytics");
library("foreach");
library("doMC")
library("doParallel")
library("iterators")
library("parallel")
library("BiocParallel")
library("raster")
library("pracma")
library("RColorBrewer")
library("reshape2")
library("ape")
library("ggplot2")
library("plyr")
library("car")
library("amap")
library("DESeq")
library("DESeq2")
library("edgeR")
library("vegan")
library("phyloseq")
################################################################################
################################################################################


################################################################################
################################################################################
#Load function definitions
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#Convert phyloseq object into edgeR object
#=> from http://joey711.github.io/phyloseq-extensions/edgeR.html
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
  }

#Convert phyloseq object to DESeq object
#=> from http://joey711.github.io/phyloseq-extensions/DESeq.html
phyloseq_to_DESeq = function(physeq, designFactor, fitType = "local", locfunc = median, ...) {
  # Enforce Orientation
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # Convert to matrix, round up to nearest integer
  x = ceiling(as(otu_table(physeq), "matrix")) + 1L
  # Add taxonomy data, if present.
  if (!is.null(tax_table(physeq, FALSE))) {
    taxADF = as(data.frame(as(tax_table(physeq), "matrix")), "AnnotatedDataFrame")
  } else {
    taxADF = NULL
  }
  # Add sample data if present
  if (!is.null(sample_data(physeq, FALSE))) {
    samplesADF = as(data.frame(sample_data(physeq)), "AnnotatedDataFrame")
  } else {
    samplesADF = NULL
  }
  # Initalize the count data sets.
  if (identical(length(designFactor), 1L)) {
    if (designFactor %in% sample_variables(physeq)) {
      designFactor <- get_variable(physeq, designFactor)
    } else {
      stop("You did not provide an appropriate `designFactor` argument. Please see documentation.")
    }
  }
  cds = newCountDataSet(x, conditions = designFactor, phenoData = samplesADF, 
                        featureData = taxADF)
  # First, estimate size factors, then estimate dispersion.
  cds = estimateSizeFactors(cds, locfunc = locfunc)
  # Now dispersions/variance estimation, passing along additional options
  cds = estimateDispersions(cds, fitType = fitType, ...)
  return(cds)
}

################################################################################
################################################################################


################################################################################
################################################################################
#Get parameters
PARAM <- list();
PARAM$cutoff.levels <- c("F", "G", "R", "S", "T");                            #Corresponding to 90%, 96%, 97%, 98% and 99% OTU clustering thresholds
PARAM$sample.size.cutoff <- 500;
PARAM$variance.threshold <- 1;
PARAM$alpha.cutoff <- 0.01;
PARAM$correlation.cutoff <- c(-0.5, 0.5);
PARAM$sig_hits.cutoff <- 15;
PARAM$tax.levels <- c("phylum", "class", "order", "family", "genus", "species");
PARAM$use.cores <- 80;
#Get working directories & file paths
PARAM$folder <- list();
PARAM$folder$base <- ""; #PUT FOLDER NAME HERE
PARAM$folder$data <- ""; #PUT FOLDER NAME HERE
PARAM$folder$results <- ""; #PUT FOLDER NAME HERE
PARAM$file <- list();
PARAM$file$sample.data <- ""; #PUT SAMPLE METADATA FILE NAME HERE
PARAM$clinical.parameters <- c("Weight_Change", "Disease_Activity_Index", "Histo_Villous_Stunting", "Histo_Epithelial_Injury", "Histo_Crypt_distortion", "Histo_IEL", "Histo_LPL_Plasma", "Histo_Eosinophils", "Histo_Morpho", "Histo_Infl", "Histo_Total", "Eosino", "IEL_50");
#Load sample metadata
sample.data <- read.table(file = PARAM$file$sample.data, header = T, sep = "\t");
rownames(sample.data) <- make.names(sample.data$Microsynth_Name);
sample.data$Time <- relevel(sample.data$Time, "post");
################################################################################
################################################################################


################################################################################
################################################################################
#Analyse "plain" taxonomy
################################################################################
#Loop through taxonomic levels
for (tax.lev in PARAM$tax.levels) {
  ################################
  #Load current taxonomy-per-sample tables
  file.tax_table <- ""; #PUT TAX TABLE FILE NAME HERE
  tmp.tax.table <- as.matrix(read.delim(file = file.tax_table, header = T, row.names = 1, sep = "\t", fill = T));
  #Normalize taxon tables
  size.sample <- colSums(tmp.tax.table);
  tax.table <- t(t(tmp.tax.table) / size.sample);
  ################################
  #Make new, global data.frame and export
  tmp.sample.data <- data.frame(sample.data, Size=size.sample[rownames(sample.data)]);
  tax.frame <- data.frame(tmp.sample.data, t(tax.table[,rownames(sample.data)]));
  write.table(tax.frame, file=paste(PARAM$folder$results, "taxonomy_overview.", tax.lev, ".tsv", sep=""), sep="\t", col.names=NA, quote=F);
  ################################
  #Filter for 8 largest taxa (incl unclassified)
  size.tax <- rowSums(tax.table, na.rm=T);
  order.tax <- unique(c(1, order(size.tax, decreasing = T)[1:7]));
  filtered.frame <- data.frame(tmp.sample.data, t(tax.table[order.tax,rownames(sample.data)]), Others=colSums(tax.table[! seq(1,nrow(tax.table)) %in% order.tax,rownames(sample.data)]));
  ################################
  #Plot stacked bar charts
  #=> plot per drug treatment, parasite-state gridded, time point gridded in 2nd dimension
  #=> plot per drug treatment, immune-state gridded, time point gridded in 2nd dimension
  #=> plot per parasite treatment, immune-state gridded, time point gridded in 2nd dimension
  ################################
  ################################
  #No drug treatment & parasite grid
  #Get current data
  curr.frame <- subset(filtered.frame, Drug == "water");
  curr.order <- order(curr.frame$Parasite, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Parasite", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.parasite_effect.treatment_water.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Parasite, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Parasite ~ Time);
  print(curr.plot);
  dev.off();
  ################################
  #DSS treatment & parasite grid
  #Get current data
  curr.frame <- subset(filtered.frame, Drug == "DSS");
  curr.order <- order(curr.frame$Parasite, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Parasite", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.parasite_effect.treatment_DSS.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Parasite, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Parasite ~ Time);
  print(curr.plot);
  dev.off();
  ################################
  ################################
  #No drug treatment & immune-state grid
  #Get current data
  curr.frame <- subset(filtered.frame, Drug == "water");
  curr.order <- order(curr.frame$Immune_State, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Immune_State", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.immune_effect.treatment_water.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Immune_State, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Immune_State ~ Time);
  print(curr.plot);
  dev.off();
  ################################
  #DSS treatment & parasite grid
  #Get current data
  curr.frame <- subset(filtered.frame, Drug == "DSS");
  curr.order <- order(curr.frame$Immune_State, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Immune_State", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.immune_effect.treatment_DSS.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Immune_State, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Immune_State ~ Time);
  print(curr.plot);
  dev.off();
  ################################
  ################################
  #No parasite treatment & immune-state grid
  #Get current data
  curr.frame <- subset(filtered.frame, Parasite == "vehicle");
  curr.order <- order(curr.frame$Immune_State, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Immune_State", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.immune_effect.treatment_vehicle.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Immune_State, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Immune_State ~ Time);
  print(curr.plot);
  dev.off();
  ################################
  #TSO treatment & parasite grid
  #Get current data
  curr.frame <- subset(filtered.frame, Parasite == "TSO");
  curr.order <- order(curr.frame$Immune_State, curr.frame$Time);
  curr.frame <- curr.frame[curr.order, ];
  curr.frame <- transform(curr.frame, Microsynth_Name=factor(Microsynth_Name, levels=curr.frame$Microsynth_Name[curr.order]));
  #Melt data in plottable "long" format
  long.frame <- melt(curr.frame, id.vars=c("Microsynth_Name", "Immune_State", "Time"), measure.vars=colnames(curr.frame)[12:length(colnames(curr.frame))], variable.name="Taxon", value.name="Rel_Abundance");
  #Generate plot
  pdf(file=paste(PARAM$folder$results, "tax_chart.immune_effect.treatment_TSO.", tax.lev, ".pdf", sep=""), width=20, height=20);
  curr.plot <- ggplot(long.frame, aes(x=interaction(Immune_State, Microsynth_Name), y=Rel_Abundance, fill=Taxon)) + geom_bar(stat="identity") + coord_flip() + facet_grid(Immune_State ~ Time);
  print(curr.plot);
  dev.off();
}
################################################################################

################################################################################
################################################################################
#Analyze data biologically
#=>	iterate through cutoff levels
#=> load data: OTU tables, OTU data, trees of representatives (coming soon...)
#=>	generate & save phyloseq objects (for later re-use)
#=> analyze alpha diversity
#=>	analyze beta diversity
#=> analyze significant associations of OTUs (EdgeR)
################################################################################
#Loop through cutoff levels
for (c in PARAM$cutoff.levels) {
  #Preallocate
  data <- list(); data$sample <- data$otu <- list(); data$otu$tax <- list();
  n <- list();
  
  ################################
  #Load current OTU table
  file.otu_table <- ""; #PUT OTU TABLE FILE NAME HERE
  otu.table <- as.matrix(read.delim(file = file.otu_table, header = T, row.names = 1));
  #Get current OTU and sample sizes
  n$otu <- nrow(otu.table);
  data$otu$size <- rowSums(otu.table);
  data$sample$size <- colSums(otu.table);
  data$sample$use <- rownames(sample.data)[data$sample$size > PARAM$sample.size.cutoff];
  #Get current OTU table w/ relative abundances
  otu.table.rel <- apply(otu.table, 1, function(i) {i / data$sample$size});
  #Load current OTU data
  file.otu_data <- ""; #PUT OTU DATA FILE NAME HERE
  otu.data <- as.matrix(read.table(file = file.otu_data, header = F, na.strings = "", stringsAsFactors=F, row.names=1, col.names = c("Name", PARAM$tax.levels, "full_taxonomy"), sep = "\t"));
  #Load current OTU rep tree
  file.otu_tree <- ""; #PUT OTU PHYLOGENY FILE NAME HERE
  otu.tree <- read_tree(file.otu_tree);
  ################################
  
  ################################
  #Create phyloseq object
  tmp.physeq <- phyloseq(otu_table(otu.table, taxa_are_rows = T), sample_data(sample.data), tax_table(otu.data));
  phy_tree(tmp.physeq) <- root(otu.tree, sample(taxa_names(tmp.physeq), 1), resolve.root=T);
  #Prune data down to include only samples that are large enough
  my.data <- prune_samples(data$sample$use, tmp.physeq);
  ################################
  
  ################################
  #Alpha diversity 
  #=>	ACE index (see http://www.mothur.org/wiki/Ace)
  #=>	Gini-Simpson index (see http://en.wikipedia.org/wiki/Diversity_index)
  #=>	Shannon index (see http://en.wikipedia.org/wiki/Diversity_index, too)
  #
  #=> test pre- and post-treatment (treatment effect) for faeces data, per genotype
  #=>	test diff in post-treatment between genotypes, per material
  ################################
  #Estimate richness / alpha diversity
  #=> phyloseq calculates several indices at a time, which is great for lazy people...
  otu.alpha_div <- estimate_richness(my.data);
  #Paste together diversity estimates with all sample data
  curr.data <- merge(sample.data, otu.alpha_div, by="row.names"); rownames(curr.data)=curr.data$Row.names;
  curr.data$Time <- factor(curr.data$Time, levels = c("pre", "post"));
  #Export current alpha div data as flat table (tsv)
  write.table(curr.data, file = paste(PARAM$folder$results, "data.alpha_div.", c, ".tsv", sep = ""), sep="\t", col.names=NA);
  ################################
  #Iterate through alpha diversity indices
  for (alpha in c("ACE", "Simpson", "InvSimpson", "Shannon")) {
    ################################
    #First group (not immune-suppressed), faeces
    #=> effect of time (paired tests)
    #=> effect of DSS treatment
    ################################
    curr.subset <- which(curr.data$Experiment %in% c(6,8) & curr.data$Material == "faeces");
    curr.frame <- curr.data[curr.subset,];
    #Plot by interaction(Parasite, Drug)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_6_8.per_time.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Drug)", y=alpha, fill = "Time")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #Plot by interaction(Parasite, Time)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_6_8.per_dss.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Time)", y=alpha, fill = "Drug")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #ANOVA
    curr.aov <- aov(curr.frame[, alpha] ~ Time*Parasite*Drug, data=curr.frame);
    capture.output(summary(curr.aov), file = paste(PARAM$folder$results, "anova.alpha_div.group_6_8.per_time.", c, ".", alpha, ".txt", sep=""));
    #Concrete tests
    #=> t-tests
    tmp <- matrix(nrow=0, ncol=5); colnames(tmp) <- c("Factor_1", "Factor_2", "Variable", "p.value", "statistic"); curr.tests <- data.frame(tmp);
    i <- 1;
    #=> paired tests between time points per drug and parasite treatment
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("water", "DSS")) {
      #Sub-select current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Drug==fac.2, c("Animal_ID", "Time", alpha)];
      #Cast into wide format for paired testing
      cast.frame <- dcast(tmp.frame, Animal_ID ~ Time, value.var=alpha);
      #Perform paired t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Time", t(unlist(wilcox.test(cast.frame$pre, cast.frame$post, paired=T)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #=> between drug treatments per parasite treatment and time point
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("pre", "post")) {
      #Subselect current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Time==fac.2, c("Drug", alpha)];
      #Perform unpaired two-sample t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Drug", t(unlist(wilcox.test(tmp.frame[, alpha]~tmp.frame$Drug)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #Export
    write.table(curr.tests, file=paste(PARAM$folder$results, "wilcox_test.alpha_div.group_6_8.", c, ".", alpha, ".tsv", sep=""), quote=F, sep="\t", row.names=F);
    ################################
    #First group (not immune-suppressed), material effect
    #=> only looking at post-treatment
    #=> comparison between materials (paired)
    #=> effect of DSS treatment
    ################################
    curr.subset <- which(curr.data$Experiment %in% c(6,8) & curr.data$Time == "post");
    curr.frame <- curr.data[curr.subset,];
    #Plot by interaction(Parasite, Drug)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_6_8.per_material.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Drug)", y=alpha, fill = "Material")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #Plot by interaction(Parasite, Material)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_6_8.per_dss_material.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Material)", y=alpha, fill = "Drug")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #ANOVA
    curr.aov <- aov(curr.frame[, alpha] ~ Material*Parasite*Drug, data=curr.frame);
    capture.output(summary(curr.aov), file = paste(PARAM$folder$results, "anova.alpha_div.group_6_8.per_material.", c, ".", alpha, ".txt", sep=""));
    #Concrete tests
    #=> t-tests
    tmp <- matrix(nrow=0, ncol=5); colnames(tmp) <- c("Factor_1", "Factor_2", "Variable", "p.value", "statistic"); curr.tests <- data.frame(tmp);
    i <- 1;
    #=> paired tests between materials per drug and parasite treatment
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("water", "DSS")) {
      #Sub-select current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Drug==fac.2, c("Animal_ID", "Material", alpha)];
      #Cast into wide format for paired testing
      cast.frame <- dcast(tmp.frame, Animal_ID ~ Material, value.var=alpha);
      #Perform paired t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Material", t(unlist(wilcox.test(cast.frame$faeces, cast.frame[, "caecal content"], paired=T)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #=> between drug treatments per parasite treatment and material
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("faeces", "caecal content")) {
      #Subselect current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Material==fac.2, c("Drug", alpha)];
      #Perform unpaired two-sample t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Drug", t(unlist(wilcox.test(tmp.frame[, alpha]~tmp.frame$Drug)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #Export
    write.table(curr.tests, file=paste(PARAM$folder$results, "wilcox_test.alpha_div.group_6_8.per_material.", c, ".", alpha, ".tsv", sep=""), quote=F, sep="\t", row.names=F);
    ################################
    #Second group (immune-suppressed), faeces
    #=> effects of time (paired tests)
    #=> effect of immune suppression
    ################################
    curr.subset <- which(curr.data$Experiment %in% c(9,10) & curr.data$Material == "faeces" & curr.data$Drug=="DSS");
    curr.frame <- curr.data[curr.subset,];
    #Plot by interaction(Parasite, Immune_State)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_9_10.per_time.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Immune_State)", y=alpha, fill = "Time")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #Plot by interaction(Parasite, Time)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_9_10.per_immune_state.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Time)", y=alpha, fill = "Immune_State")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #ANOVA
    curr.aov <- aov(curr.frame[, alpha] ~ Time*Parasite*Immune_State, data=curr.frame);
    capture.output(summary(curr.aov), file = paste(PARAM$folder$results, "anova.alpha_div.group_9_10.per_time.", c, ".", alpha, ".txt", sep=""));
    #Concrete tests
    #=> t-tests
    tmp <- matrix(nrow=0, ncol=5); colnames(tmp) <- c("Factor_1", "Factor_2", "Variable", "p.value", "statistic"); curr.tests <- data.frame(tmp);
    i <- 1;
    #=> paired tests between time points per immune_state and parasite treatment
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("control", "IS")) {
      #Sub-select current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Immune_State==fac.2, c("Animal_ID", "Time", alpha)];
      #Cast into wide format for paired testing
      cast.frame <- dcast(tmp.frame, Animal_ID ~ Time, value.var=alpha);
      #Perform paired t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Time", t(unlist(wilcox.test(cast.frame$pre, cast.frame$post, paired=T)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #=> between immune states per parasite treatment and time point
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("pre", "post")) {
      #Subselect current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Time==fac.2, c("Immune_State", alpha)];
      #Perform unpaired two-sample t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "Immune State", t(unlist(wilcox.test(tmp.frame[, alpha]~tmp.frame$Immune_State)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #Export
    write.table(curr.tests, file=paste(PARAM$folder$results, "wilcox_test.alpha_div.group_9_10.per_time.", c, ".", alpha, ".tsv", sep=""), quote=F, sep="\t", row.names=F);
    ################################
    #Second group (immune-suppressed), material effect
    #=> effects of material (paired tests) => three-way
    #=> effect of immune suppression
    ################################
    curr.subset <- which(curr.data$Experiment %in% c(9,10) & curr.data$Time == "post" & curr.data$Drug=="DSS");
    curr.frame <- curr.data[curr.subset,];
    #Plot by interaction(Parasite, Immune_State)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_9_10.per_material.", c, ".", alpha, ".pdf", sep = ""), width=10, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Immune_State)", y=alpha, fill = "Material")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #Plot by interaction(Parasite, Material)
    pdf(file = paste(PARAM$folder$results, "boxplot.alpha_div.group_9_10.per_immune_state_material.", c, ".", alpha, ".pdf", sep = ""), width=20, height=10);
    curr.plot <- ggplot(curr.frame, aes_string(x="interaction(Parasite, Material)", y=alpha, fill = "Immune_State")) + 
      geom_boxplot(alpha = 0.7, na.rm = T, outlier.colour = NA) +
      geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), size=3);
    print(curr.plot);
    dev.off();
    #ANOVA
    curr.aov <- aov(curr.frame[, alpha] ~ Material*Parasite*Immune_State, data=curr.frame);
    capture.output(summary(curr.aov), file = paste(PARAM$folder$results, "anova.alpha_div.group_9_10.per_material.", c, ".", alpha, ".txt", sep=""));
    #Concrete tests
    #=> t-tests
    tmp <- matrix(nrow=0, ncol=6); colnames(tmp) <- c("Factor_1", "Factor_2", "Variable_1", "Variable_2", "p.value", "statistic"); curr.tests <- data.frame(tmp);
    i <- 1;
    #=> paired tests between time points per immune_state and parasite treatment
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("control", "IS")) {for (var.1 in c("faeces", "caecal content", "adherent bacteria")) {for (var.2 in c("caecal content", "adherent bacteria")){
      if (var.1 == var.2) {next()}
      #Sub-select current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Immune_State==fac.2 & curr.frame$Material %in% c(var.1, var.2), c("Animal_ID", "Material", alpha)];
      #Cast into wide format for paired testing
      cast.frame <- dcast(tmp.frame, Animal_ID ~ Material, value.var=alpha);
      #Perform paired t-test
      curr.tests[i, ] <- c(fac.1, fac.2, var.1, var.2, t(unlist(wilcox.test(cast.frame[, var.1], cast.frame[, var.2], paired=T)[c("p.value", "statistic")])))
      i <- i+1;
    }}}}
    #=> between drug treatments per parasite treatment and time point
    for (fac.1 in c("vehicle", "TSO")) {for (fac.2 in c("faeces", "caecal content", "adherent bacteria")) {
      #Subselect current data
      tmp.frame <- curr.frame[curr.frame$Parasite==fac.1 & curr.frame$Material==fac.2, c("Immune_State", alpha)];
      #Perform unpaired two-sample t-test
      curr.tests[i, ] <- c(fac.1, fac.2, "control", "IS", t(unlist(wilcox.test(tmp.frame[, alpha]~tmp.frame$Immune_State)[c("p.value", "statistic")])))
      i <- i+1;
    }}
    #Export
    write.table(curr.tests, file=paste(PARAM$folder$results, "wilcox_test.alpha_div.group_9_10.per_material.", c, ".", alpha, ".tsv", sep=""), quote=F, sep="\t", row.names=F);
    ################################
    
    
    ################################
    #Test for correlations of alpha diversity with clinical parameters
    ################################
    #Preallocate tests
    #=> correlations with clinical parameters (CP)
    alpha_clin.tests <- list(); t <- 1;
    if ("test" == "test") {
      #=> group 6&8: diversity in faeces pre-treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_6_8";
      alpha_clin.tests[[t]]$name <- "faeces_pre";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "pre")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of faecal diversity pre-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Drug";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "pre treatment, faeces");
      t <- t + 1;
      #=> group 6&8: diversity in faeces post-treatement vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_6_8";
      alpha_clin.tests[[t]]$name <- "faeces_post";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of faecal diversity post-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Drug";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "post treatment, faeces");
      t <- t + 1;
      #=> group 6&8: diversity in caecum post-treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_6_8";
      alpha_clin.tests[[t]]$name <- "caecum_post";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of caecal diversity post-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Drug";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "post treatment, caecum");
      t <- t + 1;
      #=> group 6&8: change in diversity in faeces upon drug treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_6_8";
      alpha_clin.tests[[t]]$name <- "faeces_delta";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces")];
      tmp.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      new.frame <- tmp.frame[tmp.frame$Animal_ID %in% intersect(tmp.frame$Animal_ID[tmp.frame$Time == "post"], tmp.frame$Animal_ID[tmp.frame$Time == "pre"]) & tmp.frame$Time == "post", ];
      new.frame[[alpha]] <- unlist(by(tmp.frame, tmp.frame$Animal_ID, function(X) {X[X$Time == "post", alpha] - X[X$Time == "pre", alpha]}));
      alpha_clin.tests[[t]]$curr.frame <- new.frame;
      alpha_clin.tests[[t]]$curr.title <- "Correlation of the change in faecal diversity with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Drug";
      alpha_clin.tests[[t]]$curr.x_label <- paste("Change in ", alpha);
      t <- t + 1;
      #=> group 9&10: diversity in faeces pre-treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_9_10";
      alpha_clin.tests[[t]]$name <- "faeces_pre";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "pre")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of faecal diversity pre-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Immune_State";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "pre treatment, faeces");
      t <- t + 1;
      #=> group 9&10: diversity in faeces post-treatement vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_9_10";
      alpha_clin.tests[[t]]$name <- "faeces_post";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of faecal diversity post-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Immune_State";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "post treatment, faeces");
      t <- t + 1;
      #=> group 9&10: diversity in caecum post-treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_9_10";
      alpha_clin.tests[[t]]$name <- "caecum_post";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of caecal diversity post-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Immune_State";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "post treatment, caecum");
      t <- t + 1;
      #=> group 9&10: diversity in adherent bacteria post-treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_9_10";
      alpha_clin.tests[[t]]$name <- "adherent_post";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post")];
      alpha_clin.tests[[t]]$curr.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      alpha_clin.tests[[t]]$curr.title <- "Correlation of adherent bacteria diversity post-treatment with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Immune_State";
      alpha_clin.tests[[t]]$curr.x_label <- paste(alpha, "post treatment, adherent bacteria");
      t <- t + 1;
      #=> group 9&10: change in diversity in faeces upon drug treatment vs CP
      alpha_clin.tests[[t]] <- list();
      alpha_clin.tests[[t]]$group <- "group_9_10";
      alpha_clin.tests[[t]]$name <- "faeces_delta";
      alpha_clin.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces")];
      tmp.frame <- curr.data[alpha_clin.tests[[t]]$curr.subset, ];
      new.frame <- tmp.frame[tmp.frame$Animal_ID %in% intersect(tmp.frame$Animal_ID[tmp.frame$Time == "post"], tmp.frame$Animal_ID[tmp.frame$Time == "pre"]) & tmp.frame$Time == "post", ];
      new.frame[[alpha]] <- unlist(by(tmp.frame, tmp.frame$Animal_ID, function(X) {X[X$Time == "post", alpha] - X[X$Time == "pre", alpha]}));
      alpha_clin.tests[[t]]$curr.frame <- new.frame;
      alpha_clin.tests[[t]]$curr.title <- "Correlation of the change in faecal diversity with ";
      alpha_clin.tests[[t]]$curr.color <- "Parasite";
      alpha_clin.tests[[t]]$curr.shape <- "Immune_State";
      alpha_clin.tests[[t]]$curr.x_label <- paste("Change in ", alpha);
      t <- t + 1;
    }
    
    #Iterate through clinical parameters
    for (cp in PARAM$clinical.parameters) {
      #Iterate through tests
      for (t in seq(1, length(alpha_clin.tests))) {
        #Start with a clean sheet
        rm(curr.frame, curr.subset, name, group);
        #Get current test parameters
        attach(alpha_clin.tests[[t]]);
        #Get interaction terms for current variable parameters
        curr.frame$Interaction <- interaction(curr.frame[[curr.color]], curr.frame[[curr.shape]]);
        curr.interactions <- levels(curr.frame$Interaction);
        #Get current correlations per interaction term	
        curr.cor <- curr.sig <- numeric(length=4); names(curr.cor) <- names(curr.sig) <- curr.interactions; cc <- 1;
        for (ia in curr.interactions) {
          tmp.test <- cor.test(curr.frame[[alpha]][curr.frame$Interaction == ia], curr.frame[[cp]][curr.frame$Interaction == ia], method = "spearman");
          curr.cor[cc] <- unname(tmp.test$estimate); curr.sig[cc] <- unname(tmp.test$p.value);
          cc <- cc+1;
        }
        #Plot current correlations
        pdf(file=paste(PARAM$folder$results, paste("alpha_div.clinical_parameters", c, alpha, cp, group, name, "pdf", sep="."), sep=""), width=10, height=10);
        curr.plot <- ggplot(curr.frame, aes_string(x=alpha, y=cp, shape="Interaction", color="Interaction")) +
          scale_color_brewer(type="div", palette="PuOr") +
          geom_point(size=5) +
          geom_smooth(method=lm, se=F, fullrange=F) +
          xlab(curr.x_label) +
          ggtitle(paste(paste(curr.title, cp), paste(sapply(curr.interactions, function(x) {paste(x, curr.cor[x], curr.sig[x])}), collapse="\n"), sep="\n"));
        print(curr.plot);
        dev.off();
        #Tidy up
        detach(alpha_clin.tests[[t]]);
      }
    }
  }
  ################################
  
  
  ################################
  #Beta diversity
  #=> weighted UniFrac (see http://en.wikipedia.org/wiki/UniFrac & http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317376/)
  #=> Bray-Curtis dissimilarity (see http://en.wikipedia.org/wiki/Bray???Curtis_dissimilarity)
  #
  #=> plot PCoA ordinations
  #=> at each stage, assess significance of separation in PERMANOVA tests (vegan::adonis)
  ################################
  #Calculate pairwise sample distances
  #Register parallel backend (for UniFrac calculation)
  registerDoParallel(cores = PARAM$use.cores);
  #Calculate pairwise sample distances as weighted UniFrac & Bray-Curtis
  file.dist <- paste(PARAM$folder$results, "data.beta_div.distances.", c, ".RData", sep = "");
  #Check if a pre-calculated distance file is available
  if (file.exists(file.dist)) {load(file.dist)} else {
    data$sample$dist <- list();
    data$sample$dist$unifrac <- UniFrac(my.data, weighted = T, normalized = T, parallel = T, fast = T);
    data$sample$dist$bray_curtis <- distance(my.data, method = "bray");
    #Save current distances
    save(data, file=file.dist);
  }
  ################################
  #Preallocate tests to perform
  beta_div.tests <- list(); t <- 1;
  if ("test" == "test") {
    ################################
    #First group (exp 6 & 8)
    ################################
    #Make one global overview plot and test contribution of all factors (room, TSO treatment, DSS treatment, material)
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "all_factors";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8))];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Global ordination of all data, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Room + Parasite*Drug*Time + Material";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS & water treatment (post vs pre) on TSO(+/-) rabbits, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Drug*Time";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in TSO-treated rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.TSO.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "TSO")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS & water treatment (post vs pre) on TSO-treated rabbits, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Drug*Time";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in non-TSO-treated rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.non_TSO.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "vehicle")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS & water treatment (post vs pre) on TSO-treated rabbits, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Drug*Time";
    t <- t + 1;
    #Treatment effect in faeces post treatment: does DSS treatment (+/-) have an effect on microbiota composition in TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.post";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug)
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS(+/-) and TSO(+/-) in faeces post treatment, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Drug";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Drug";
    t <- t + 1;
    #Treatment effect in caecum post treatment: does DSS treatment (+/-) have an effect on microbiota composition in TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_6_8";
    beta_div.tests[[t]]$name <- "treatment_effect.caecum.post";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Drug)
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS(+/-) and TSO(+/-) in caecum post treatment, experiments 6 & 8";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Drug";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Drug";
    t <- t + 1;
    ################################
    
    
    ################################
    #Second group (exp 9 & 10)
    ################################
    #Make one global overview plot and test contribution of all factors (room, TSO treatment, Immune_State, material)
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "all_factors";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10))];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Global ordination of all data, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Room + Parasite*Immune_State*Drug*Time + Material";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in IS(+/-) TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment (post vs pre) on IS(+/-) TSO(+/-) rabbits, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Immune_State*Time";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in TSO-treated IS(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.TSO.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment (post vs pre) on IS(+/-) TSO rabbits, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Immune_State*Time";
    t <- t + 1;
    #Treatment effect in faeces: does DSS treatment (pre vs post) have an effect on microbiota composition in non-TSO-treated IS(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.TSO.post_vs_pre";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment (post vs pre) on IS(+/-) non-TSO rabbits, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Time";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Immune_State*Time";
    t <- t + 1;
    #Treatment effect in faeces post treatment: does DSS treatment have an effect on microbiota composition in IS(+/-) TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.faeces.post";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Time == "post")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment on IS(+/-) TSO(+/-) rabbits in faeces, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Immune_State";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Immune_State";
    t <- t + 1;
    #Treatment effect in caecum post treatment: does DSS treatment have an effect on microbiota composition in IS(+/-) TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.caecum.post";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Drug == "DSS" & sample.data$Time == "post")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment on IS(+/-) TSO(+/-) rabbits in caecum, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Immune_State";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Immune_State";
    t <- t + 1;
    #Treatment effect in adherent bacteria post treatment: does DSS treatment have an effect on microbiota composition in IS(+/-) TSO(+/-) rabbits?
    beta_div.tests[[t]] <- list();
    beta_div.tests[[t]]$group <- "group_9_10";
    beta_div.tests[[t]]$name <- "treatment_effect.adherent_bacteria.post";
    beta_div.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Drug == "DSS" & sample.data$Time == "post")];
    beta_div.tests[[t]]$curr.physeq <- prune_samples(beta_div.tests[[t]]$curr.subset, tmp.physeq);
    tmp.sample_data <- as(sample_data(beta_div.tests[[t]]$curr.physeq), "data.frame");
    tmp.sample_data$Interaction <- interaction(tmp.sample_data$Parasite, tmp.sample_data$Immune_State);
    beta_div.tests[[t]]$curr.data <- tmp.sample_data;
    beta_div.tests[[t]]$curr.title <- "Effect of DSS treatment on IS(+/-) TSO(+/-) rabbits in adherent bacteria, experiments 9 & 10";
    beta_div.tests[[t]]$curr.color <- "Interaction";
    beta_div.tests[[t]]$curr.shape <- "Immune_State";
    beta_div.tests[[t]]$permanova.formula <- "curr.dist ~ Parasite*Immune_State";
    t <- t + 1;
  }
  ################################
  #Iterate through tests
  for (t in seq(1, length(beta_div.tests))) {
    #Get current test parameters
    attach(beta_div.tests[[t]]);
    #Remove possbily attached stuff (just to be sure)
    rm(curr.data, curr.physeq, name, group);
    for (dist.fun in c("unifrac", "bray_curtis")) {
      curr.dist <- as.dist(as.matrix(data$sample$dist[[dist.fun]])[curr.subset, curr.subset]);
      #Plot ordination
      curr.pcoa <- ordinate(curr.physeq, method = "PCoA", distance=curr.dist);
      #Prepare data for plotting
      plot.data <- data.frame(curr.data, Axis.1=curr.pcoa$vectors[, "Axis.1"], Axis.2=curr.pcoa$vectors[, "Axis.2"], Group=interaction(curr.data[[curr.color]], curr.data[[curr.shape]]));
      plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"], digits=1);
      plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
      plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
      curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
      #Make plot
      pdf(file=paste(PARAM$folder$results, paste("ordination", c, group, name, dist.fun, "pdf", sep="."), sep=""), width=20, height=20);
      curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", color=curr.color, shape=curr.shape)) +
        geom_point(size=5) +
        geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
        geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
        geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
        ggtitle(curr.title) +
        xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
        ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep=""));
      print(curr.plot);
      dev.off();
      #PERMANOVA
      tmp.aov <- adonis(as.formula(permanova.formula), data=curr.data);
      capture.output(tmp.aov, file = paste(PARAM$folder$results, paste("permanova", c, group, name, dist.fun, "txt", sep="."), sep=""));
    }
    #Tidy up
    detach(beta_div.tests[[t]]);
  }
  #Be sure to detach, even if loop crashes for some reason...
  if ("beta_div.tests[[t]]" %in% search()) {detach(beta_div.tests[[t]])}
  ################################
  
  
  ################################
  #Significantly associated OTUs
  #=> EdgeR-based analysis
  #=> Assess which OTUs are significantly associated to different treatments and/or states
  #################################
  #Pre-remove all OTUs with trivially low variance
  tmp.var <- apply(otu_table(tmp.physeq), 1, var);
  keep.OTUs = names(which(tmp.var > PARAM$variance.threshold));
  pre_filt.physeq <- phyloseq(otu_table(tmp.physeq)[keep.OTUs,], sample_data(tmp.physeq), tax_table(tmp.physeq)[keep.OTUs,]);
  ################################
  #Preallocate tests to perform
  assoc.tests <- list(); t <- 1;
  if ("test" == "test") {
    ################################
    #First group (exp 6 & 8)
    ################################
    #Paired tests
    #=> which OTUs are associated to time (ageing) in non-TSO, non-DSS rabbits (full controls)?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "paired.faeces.non_TSO.non_DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "vehicle" & sample.data$Drug == "water")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in TSO, non-DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "paired.faeces.TSO.non_DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "TSO" & sample.data$Drug == "water")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in TSO, DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "paired.faeces.TSO.DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "TSO" & sample.data$Drug == "DSS")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in non-TSO, DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "paired.faeces.non_TSO.DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Parasite == "vehicle" & sample.data$Drug == "DSS")];
    t <- t + 1;
    ################################
    #"Classic" (unpaired) tests
    #post vs pre
    #=> which OTUs are associated to non-TSO, non-DSS treatment in faeces (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.non_TSO.non_DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "pre")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "vehicle" & sample.data$Drug == "water")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to TSO, non-DSS treatment in faeces (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.TSO.non_DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "pre")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "TSO" & sample.data$Drug == "water")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to TSO, DSS treatment in faeces (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.TSO.DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "pre")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "TSO" & sample.data$Drug == "DSS")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to non-TSO, DSS treatment in faeces (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.non_TSO.DSS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "pre")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "vehicle" & sample.data$Drug == "DSS")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    ################################
    #Conditions post treatment, faeces
    #=> which OTUs are associated to DSS in non-TSO rabbits, in faeces post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.non_TSO.post.DSS_vs_water";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "vehicle")];
    assoc.tests[[t]]$group_by <- "Drug";
    assoc.tests[[t]]$baseline <- "water";
    assoc.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #=> which OTUs are associated to DSS in TSO rabbits, in faeces post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.TSO.post.DSS_vs_water";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "TSO")];
    assoc.tests[[t]]$group_by <- "Drug";
    assoc.tests[[t]]$baseline <- "water";
    assoc.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #=> which OTUs are associated to TSO in non-DSS rabbits, in faeces post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.non_DSS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "water")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    #=> which OTUs are associated to TSO in DSS rabbits, in faeces post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    ################################
    #Conditions post treatment, caecum
    #=> which OTUs are associated to DSS in non-TSO rabbits, in caecum post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.caecum.non_TSO.post.DSS_vs_water";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "vehicle")];
    assoc.tests[[t]]$group_by <- "Drug";
    assoc.tests[[t]]$baseline <- "water";
    assoc.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #=> which OTUs are associated to DSS in TSO rabbits, in caecum post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.caecum.TSO.post.DSS_vs_water";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "TSO")];
    assoc.tests[[t]]$group_by <- "Drug";
    assoc.tests[[t]]$baseline <- "water";
    assoc.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #=> which OTUs are associated to TSO in non-DSS rabbits, in caecum post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.caecum.non_DSS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "water")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    #=> which OTUs are associated to TSO in DSS rabbits, in caecum post treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_6_8";
    assoc.tests[[t]]$name <- "classic.caecum.DSS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "DSS")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    ################################
    
    
    ################################
    #Group 9 & 10
    ################################
    #Paired tests
    #=> which OTUs are associated to time (ageing) in non-TSO, non-IS DSS rabbits (full controls)?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "paired.faeces.DSS.non_TSO.non_IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle" & sample.data$Immune_State == "control")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in TSO, non-IS DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "paired.faeces.DSS.TSO.non_IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO" & sample.data$Immune_State == "control")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in TSO, IS DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "paired.faeces.DSS.TSO.IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO" & sample.data$Immune_State == "IS")];
    t <- t + 1;
    #=> which OTUs are associated to time (ageing) in non-TSO, IS DSS rabbits?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "paired";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "paired.faeces.DSS.non_TSO.IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle" & sample.data$Immune_State == "IS")];
    t <- t + 1;
    ################################
    #"Classic" (unpaired) tests
    #post vs pre
    #=> which OTUs are associated to non-TSO, non-IS treatment in faeces & DSS (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.non_TSO.non_IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "pre" & sample.data$Drug == "DSS")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle" & sample.data$Immune_State == "control")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to TSO, non-IS treatment in faeces & DSS (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.TSO.non_IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "pre" & sample.data$Drug == "DSS")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO" & sample.data$Immune_State == "control")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to TSO, IS treatment in faeces & DSS (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.TSO.IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "pre" & sample.data$Drug == "DSS")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO" & sample.data$Immune_State == "IS")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    #=> which OTUs are associated to non-TSO, IS treatment in faeces & DSS (unpaired against combined baseline)
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.non_TSO.IS.post_vs_pre";
    assoc.tests[[t]]$curr.subset <- c(
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "pre" & sample.data$Drug == "DSS")],
      rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle" & sample.data$Immune_State == "IS")]
    );
    assoc.tests[[t]]$group_by <- "Time";
    assoc.tests[[t]]$baseline <- "pre";
    assoc.tests[[t]]$response <- "post";
    t <- t + 1;
    ################################
    #Conditions post treatment, faeces
    #=> which OTUs are associated to IS in non-TSO rabbits, in faeces post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.non_TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to IS in TSO rabbits, in faeces post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to TSO in non-IS rabbits, in faeces post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.non_IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "control")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    #=> which OTUs are associated to TSO in IS rabbits, in faeces post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.faeces.DSS.IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "IS")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    ################################
    #Conditions post treatment, caecum
    #=> which OTUs are associated to IS in non-TSO rabbits, in caecum post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.caecum.DSS.non_TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to IS in TSO rabbits, in caecum post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.caecum.DSS.TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to TSO in non-IS rabbits, in caecum post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.caecum.DSS.non_IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "control")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    #=> which OTUs are associated to TSO in IS rabbits, in caecum post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.caecum.DSS.IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "IS")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    ################################
    #Conditions post treatment, adherent bacteria
    #=> which OTUs are associated to IS in non-TSO rabbits, in adherent bacteria post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.adherent.DSS.non_TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "vehicle")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to IS in TSO rabbits, in adherent bacteria post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.adherent.DSS.TSO.post.IS_vs_control";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Parasite == "TSO")];
    assoc.tests[[t]]$group_by <- "Immune_State";
    assoc.tests[[t]]$baseline <- "control";
    assoc.tests[[t]]$response <- "IS";
    t <- t + 1;
    #=> which OTUs are associated to TSO in non-IS rabbits, in adherent bacteria post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.adherent.DSS.non_IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "control")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
    #=> which OTUs are associated to TSO in IS rabbits, in adherent bacteria post DSS treatment?
    assoc.tests[[t]] <- list();
    assoc.tests[[t]]$type <- "classic";
    assoc.tests[[t]]$group <- "group_9_10";
    assoc.tests[[t]]$name <- "classic.adherent.DSS.IS.post.TSO_vs_vehicle";
    assoc.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post" & sample.data$Drug == "DSS" & sample.data$Immune_State == "IS")];
    assoc.tests[[t]]$group_by <- "Parasite";
    assoc.tests[[t]]$baseline <- "vehicle";
    assoc.tests[[t]]$response <- "TSO";
    t <- t + 1;
  }
  
  ################################
  #Iterate through tests
  for (t in seq(1, length(assoc.tests))) {
    #Get current data
    attach(assoc.tests[[t]]);
    
    #Tidy up
    rm(curr.subset, curr.physeq, name, group);
    
    #Switch between "paired" (time-dependent) and "unpaired" cases
    if (type == "paired") {
      #Subset current data
      curr.physeq <- prune_samples(curr.subset, pre_filt.physeq);
      #Get variables
      tmp.animal_id <- get_variable(curr.physeq, "Animal_ID");
      tmp.time <- get_variable(curr.physeq, "Time");
      #Focus on animals for which data is available at both time points
      keep.animals = names(which(tapply(tmp.time, tmp.animal_id, function(x){length(unique(x))}) == 2));
      paired.physeq <- subset_samples(curr.physeq, Animal_ID %in% keep.animals);
      animal_id <- get_variable(paired.physeq, "Animal_ID");
      time <- get_variable(paired.physeq, "Time");
      #Prefilter to remove OTUs with trivially low variance
      curr.var <- apply(otu_table(paired.physeq), 1, var);
      keep.OTUs = names(which(curr.var > PARAM$variance.threshold));
      pruned.physeq <- phyloseq(otu_table(paired.physeq)[keep.OTUs,], sample_data(paired.physeq), tax_table(paired.physeq)[keep.OTUs,]);
      #Add one to protect against overflow, log(0) issues
      x <- as(otu_table(pruned.physeq), "matrix") + 1L;
      taxonomy <- data.frame(as(tax_table(pruned.physeq), "matrix"));
      #Coerce into EdgeR-understandable DGEList object
      x <- DGEList(counts=x, group=time, genes=taxonomy, remove.zeros=T);
      #Get design matrix
      design <- model.matrix(~ animal_id + time);
      #Estimate dispoersion
      x <- calcNormFactors(x, method="RLE");
      x <- estimateGLMCommonDisp(x, design);
      x <- estimateGLMTrendedDisp(x, design);
      x <- estimateGLMTagwiseDisp(x, design);
      #Fit GLM & detect differentially abundant OTUs
      fit <- glmFit(x, design);
      lrt <- glmLRT(fit);
      #Extract results
      tt <- topTags(lrt, n=nrow(x), adjust.method="BH", sort.by="PValue");
    } else {
      #Subset current data
      curr.physeq <- prune_samples(curr.subset, pre_filt.physeq);
      #Prefilter to remove OTUs with trivially low variance
      curr.var <- apply(otu_table(curr.physeq), 1, var);
      keep.OTUs = names(which(curr.var > PARAM$variance.threshold));
      pruned.physeq <- phyloseq(otu_table(curr.physeq)[keep.OTUs,], sample_data(curr.physeq), tax_table(curr.physeq)[keep.OTUs,]);
      #Coerce into DGEList object
      x <- phyloseq_to_edgeR(pruned.physeq, group=group_by);
      #Perform exact test
      et <- exactTest(x, pair=c(baseline, response));
      #Extract results
      tt <- topTags(et, n=nrow(x), adjust.method="BH", sort.by="PValue");
    }
    
    #Make volcano plot
    #Prepare data
    volcano.data <- data.frame(PValue=-log10(tt$table$PValue), LogFC=tt$table$logFC, Significance=as.factor(tt$table$FDR < PARAM$alpha.cutoff));
    #Plot
    pdf(file=paste(PARAM$folder$results, paste("volcano", c, group, name, "pdf", sep="."), sep=""), width=10, height=10);
    curr.plot <- ggplot(volcano.data, aes(x=LogFC, y=PValue, colour=Significance)) + xlab("log2 fold change") + ylab("-log10 p-value");
    print(curr.plot + geom_point(size = 5));
    dev.off();
    #Export all topTags
    curr.ot <- otu_table(pruned.physeq);
    curr.otu_sizes = rowSums(curr.ot);
    curr.inc_ot <- curr.ot; curr.inc_ot[curr.inc_ot > 0] <- 1;
    curr.observed_in_samples <- rowSums(curr.inc_ot);
    curr.sample_count <- ncol(curr.ot)
    curr.res <- data.frame(tt$table, Global_OTU_Size=data$otu$size[rownames(tt$table)], Local_OTU_Size=curr.otu_sizes[rownames(tt$table)], Observed_in_Samples=curr.observed_in_samples[rownames(tt$table)], Local_Sample_Count=curr.sample_count, Significance=volcano.data$Significance);
    write.table(curr.res, file=paste(PARAM$folder$results, paste("assoc_table", c, group, name, "tsv", sep="."), sep=""), sep="\t", col.names=NA);
    #Make taxonomy-based plot
    theme_set(theme_bw());
    scale_fill_discrete <- function(palname = "Set1", ...) {
      scale_fill_brewer(palette = palname, ...)
    }
    sigtab = curr.res[(curr.res$FDR < PARAM$alpha.cutoff), ];
    sigtab$genus <- factor(sigtab$genus, levels=c(levels(sigtab$genus), "Unclassified")); sigtab$genus[is.na(sigtab$genus)] <- "Unclassified";
    sigtab$phylum <- factor(sigtab$phylum, levels=c(levels(sigtab$phylum), "Unclassified")); sigtab$phylum[is.na(sigtab$phylum)] <- "Unclassified";
    plot.tab = tapply(sigtab$logFC, sigtab$genus, function(x) max(x))
    plot.tab = sort(plot.tab, TRUE);
    sigtab$genus = factor(as.character(sigtab$genus), levels = names(plot.tab));
    if (nrow(sigtab) > 0) {
      pdf(file=paste(PARAM$folder$results, paste("otu_association.taxonomy", c, group, name, "pdf", sep="."), sep=""), width=10, height=10);
      curr.plot <- ggplot(sigtab, aes(x = genus, y = logFC, color = phylum, fill=phylum)) + geom_point(size=5, position=position_jitterdodge(dodge.width=0.9)) + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
      print (curr.plot);
      dev.off();
    }
    
    #Store current topTags
    assoc.tests[[t]]$tt <- curr.res;
    
    #Tidy up
    detach(assoc.tests[[t]]);
  }
  #Be sure to detach, even if loop crashes for some reason...
  if ("assoc.tests[[t]]" %in% search()) {detach(assoc.tests[[t]])}
  ################################
  
  
  ################################
  #Correlation of OTU abundances with clinical parameters
  #################################
  #Get filtered, normalized (relative) OTU table
  norm.ot <- as(t(t(otu_table(pre_filt.physeq)) / colSums(otu_table(pre_filt.physeq))), "matrix");
  #Preallocate tests for correlations with clinical parameters
  correlation.tests <- list(); t <- 1;
  if ("test" == "test") {
    ################################
    #First group (exp 6 & 8)
    ################################
    #Correlations in faeces post treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "unpaired";
    correlation.tests[[t]]$group <- "group_6_8";
    correlation.tests[[t]]$name <- "unpaired.faeces.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post")];
    t <- t + 1;
    #Correlations in caecum post treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "unpaired";
    correlation.tests[[t]]$group <- "group_6_8";
    correlation.tests[[t]]$name <- "unpaired.caecum.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post")];
    t <- t + 1;
    #Correlations with change in faecal abundance post vs pre treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "paired";
    correlation.tests[[t]]$group <- "group_6_8";
    correlation.tests[[t]]$name <- "paired.faeces.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces")];
    t <- t + 1;
    ################################
    #Second group (exp 9 & 10)
    ################################
    #Correlations in faeces post treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "unpaired";
    correlation.tests[[t]]$group <- "group_9_10";
    correlation.tests[[t]]$name <- "unpaired.faeces.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post")];
    t <- t + 1;
    #Correlations in caecum post treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "unpaired";
    correlation.tests[[t]]$group <- "group_9_10";
    correlation.tests[[t]]$name <- "unpaired.caecum.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post")];
    t <- t + 1;
    #Correlations in adherent bacteria post treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "unpaired";
    correlation.tests[[t]]$group <- "group_9_10";
    correlation.tests[[t]]$name <- "unpaired.adherent_bacteria.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "adherent bacteria" & sample.data$Time == "post")];
    t <- t + 1;
    #Correlations with change in faecal abundance post vs pre treatment
    correlation.tests[[t]] <- list();
    correlation.tests[[t]]$type <- "paired";
    correlation.tests[[t]]$group <- "group_9_10";
    correlation.tests[[t]]$name <- "paired.faeces.post";
    correlation.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces")];
    t <- t + 1;
  }
  
  ################################
  #Iterate through tests
  for (t in seq(1, length(correlation.tests))) {
    attach(correlation.tests[[t]]);
    #Tidy up
    rm(curr.subset, curr.physeq, name, group);
    #Subset current data
    curr.samples <- which(rownames(sample.data) %in% curr.subset);
    curr.ot <- norm.ot[, curr.subset];
    curr.otus <- rownames(curr.ot);
    #Preallocate results data.frame
    curr.res <- data.frame(otu.data[curr.otus, ], Global_OTU_Size=data$otu$size[curr.otus]);
    #Iterate through clinical parameters
    for (cp in PARAM$clinical.parameters) {
      #Switch between "paired" and "unpaired" tests
      if (type == "paired") {
        #Get variables
        tmp.animal_id <- sample.data$Animal_ID[curr.samples];
        tmp.time <- sample.data$Time[curr.samples];
        #Focus on animals for which data is available at both time points
        keep.animals = names(which(tapply(tmp.time, tmp.animal_id, function(x){length(unique(x))}) == 2));
        #Get samples pre and post treatment per animal (pairs)
        paired.samples <- data.frame(
          pre=rownames(sample.data)[sample.data$Animal_ID %in% keep.animals & sample.data$Time == "pre" & sample.data$Material == "faeces"],
          post=rownames(sample.data)[sample.data$Animal_ID %in% keep.animals & sample.data$Time == "post" & sample.data$Material == "faeces"],
          row.names=keep.animals
        );
        #Extract relevant values for current clinical parameter
        curr.cp <- sample.data[[cp]][sample.data$Animal_ID %in% keep.animals & sample.data$Time == "post" & sample.data$Material == "faeces"];
        #Get table of log-transformed changes OTU abundances
        #=> log2(post / pre);
        #Add small value to protect from log(0) overflow
        tmp.ot <- curr.ot + 10^-10;
        tmp.paired <- mclapply(rownames(tmp.ot), function(otu) {log2(tmp.ot[otu, as.character(paired.samples$post)] / tmp.ot[otu, as.character(paired.samples$pre)])}, mc.cores=PARAM$use.cores);
        #Compute correlations
        tmp.cor <- unlist(mclapply(tmp.paired, function(X) {cor(X, curr.cp, method="spearman", use="complete.obs")}, mc.cores=PARAM$use.cores));
        
      } else {
        #Compute correlations
        tmp.cor <- unlist(mclapply(rownames(curr.ot), function(otu) {cor(curr.ot[otu,], sample.data[[cp]][curr.samples], method="spearman", use="complete.obs")}, mc.cores=PARAM$use.cores));
      }
      names(tmp.cor) <- rownames(curr.ot);
      #Store to results data.frame
      curr.res[[cp]] <- tmp.cor;
      
      #Plot correlated OTUs, sorted by taxonomy#Make taxonomy-based plot
      theme_set(theme_bw());
      scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
      }
      sigtab = curr.res[!is.na(curr.res[[cp]]) & (curr.res[[cp]] < PARAM$correlation.cutoff[1] | curr.res[[cp]] > PARAM$correlation.cutoff[2]), ];
      sigtab$genus <- factor(sigtab$genus, levels=c(levels(sigtab$genus), "Unclassified")); sigtab$genus[is.na(sigtab$genus)] <- "Unclassified";
      sigtab$phylum <- factor(sigtab$phylum, levels=c(levels(sigtab$phylum), "Unclassified")); sigtab$phylum[is.na(sigtab$phylum)] <- "Unclassified";
      plot.tab = tapply(sigtab[[cp]], sigtab$genus, function(x) max(x))
      plot.tab = sort(plot.tab, TRUE);
      sigtab$genus = factor(as.character(sigtab$genus), levels = names(plot.tab));
      if (nrow(sigtab) > 0) {
        pdf(file=paste(PARAM$folder$results, paste("otu_correlation.taxonomy", c, cp, group, name, "pdf", sep="."), sep=""), width=10, height=10);
        curr.plot <- ggplot(sigtab, aes_string(x = "genus", y = cp, color = "phylum", fill="phylum")) + geom_point(size=5, position=position_jitterdodge(dodge.width=0.9)) + 
          theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
        print (curr.plot);
        dev.off();
      }
      
      #Make "correlation volcano"
      #Prepare data
      volcano.data <- data.frame(Correlation=curr.res[[cp]], LogGlobalSize=log10(curr.res$Global_OTU_Size), Significance=as.factor(!is.na(curr.res[[cp]]) & (curr.res[[cp]] < PARAM$correlation.cutoff[1] | curr.res[[cp]] > PARAM$correlation.cutoff[2])));
      #Plot
      pdf(file=paste(PARAM$folder$results, paste("otu_correlation.volcano", c, cp, group, name, "pdf", sep="."), sep=""), width=10, height=10);
      curr.plot <- ggplot(volcano.data, aes(x=Correlation, y=LogGlobalSize, colour=Significance)) + xlab("Correlation of OTU Abundance with Clinical Parameter") + ylab("log10(Global OTU Size)");
      print(curr.plot + geom_point(size = 5, alpha = .5));
      dev.off();
    }
    
    #Export results data.frame
    write.table(curr.res, file=paste(PARAM$folder$results, paste("otu_correlation.table", c, cp, group, name, "tsv", sep="."), sep=""), sep="\t", col.names=NA);
    
    #Tidy up
    detach(correlation.tests[[t]]);
  }
  #Be sure to detach, even if loop crashes for some reason...
  if ("correlation.tests[[t]]" %in% search()) {detach(correlation.tests[[t]])}
  ################################
  
  
  ################################
  #Collect data for an OTU heatmap
  #=> significantly associated OTUs across selected conditions
  #=> correlations with clinical parameters in faeces and caecum
  ################################
  #Preallocate tests and results collector
  collect.tests <- list(); t <- 1;
  collect.otus <- character();
  raw.collect.otus <- character();
  collect.data <- matrix(nrow=nrow(norm.ot), ncol=16);
  rownames(collect.data) <- rownames(norm.ot);
  tmp.colnames <- character(length=16);
  if ("test" == "test") {
    ################################
    #Associations (fold changes)
    ################################
    #Caecum, exp. 6&8, (TSO-, IS-, DSS+) vs (TSO-, IS-, DSS-)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.6_8.non_TSO.non_IS.post.DSS_vs_water";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "vehicle")];
    collect.tests[[t]]$group_by <- "Drug";
    collect.tests[[t]]$baseline <- "water";
    collect.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #Caecum, exp. 6&8, (TSO+, IS-, DSS+) vs (TSO-, IS-, DSS-)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.6_8.TSO.non_IS.post.DSS_vs_water";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "TSO")];
    collect.tests[[t]]$group_by <- "Drug";
    collect.tests[[t]]$baseline <- "water";
    collect.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #Caecum, exp. 9&10, (TSO-, IS+, DSS+) vs (TSO-, IS-, DSS+)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.9_10.non_TSO.DSS.post.IS_vs_IC";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "vehicle" & sample.data$Drug == "DSS")];
    collect.tests[[t]]$group_by <- "Immune_State";
    collect.tests[[t]]$baseline <- "control";
    collect.tests[[t]]$response <- "IS";
    t <- t + 1;
    #Caecum, exp. 9&10, (TSO+, IS+, DSS+) vs (TSO-, IS-, DSS+)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.9_10.TSO.DSS.post.IS_vs_IC";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "caecal content" & sample.data$Time == "post" & sample.data$Parasite == "TSO" & sample.data$Drug == "DSS")];
    collect.tests[[t]]$group_by <- "Immune_State";
    collect.tests[[t]]$baseline <- "control";
    collect.tests[[t]]$response <- "IS";
    t <- t + 1;
    #Faeces, exp. 6&8, (TSO-, IS-, DSS+) vs (TSO-, IS-, DSS-)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.6_8.non_TSO.non_IS.post.DSS_vs_water";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "vehicle")];
    collect.tests[[t]]$group_by <- "Drug";
    collect.tests[[t]]$baseline <- "water";
    collect.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #Faeces, exp. 6&8, (TSO+, IS-, DSS+) vs (TSO-, IS-, DSS-)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.6_8.TSO.non_IS.post.DSS_vs_water";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(6,8) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "TSO")];
    collect.tests[[t]]$group_by <- "Drug";
    collect.tests[[t]]$baseline <- "water";
    collect.tests[[t]]$response <- "DSS";
    t <- t + 1;
    #Faeces, exp. 9&10, (TSO-, IS+, DSS+) vs (TSO-, IS-, DSS+)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.9_10.non_TSO.DSS.post.IS_vs_IC";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "vehicle" & sample.data$Drug == "DSS")];
    collect.tests[[t]]$group_by <- "Immune_State";
    collect.tests[[t]]$baseline <- "control";
    collect.tests[[t]]$response <- "IS";
    t <- t + 1;
    #Faeces, exp. 9&10, (TSO+, IS+, DSS+) vs (TSO-, IS-, DSS+)
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "association";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.9_10.TSO.DSS.post.IS_vs_IC";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Experiment %in% c(9,10) & sample.data$Material == "faeces" & sample.data$Time == "post" & sample.data$Parasite == "TSO" & sample.data$Drug == "DSS")];
    collect.tests[[t]]$group_by <- "Immune_State";
    collect.tests[[t]]$baseline <- "control";
    collect.tests[[t]]$response <- "IS";
    t <- t + 1;
    ################################
    #Correlations
    ################################
    #Caecum, exp. 6&8&9&10, correlation with weight change
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.post.weight_change";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "caecal content" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Weight_Change";
    t <- t + 1;
    #Caecum, exp. 6&8&9&10, correlation with disease activity index
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.post.DAI";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "caecal content" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Disease_Activity_Index";
    t <- t + 1;
    #Caecum, exp. 6&8&9&10, correlation with Histo_Total
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.post.histo";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "caecal content" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Histo_Total";
    t <- t + 1;
    #Caecum, exp. 6&8&9&10, correlation with Eosino
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "caecum.post.eosino";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "caecal content" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Eosino";
    t <- t + 1;
    #Faeces, exp. 6&8&9&10, correlation with weight change
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.post.weight_change";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "faeces" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Weight_Change";
    t <- t + 1;
    #Faeces, exp. 6&8&9&10, correlation with disease activity index
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.post.DAI";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "faeces" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Disease_Activity_Index";
    t <- t + 1;
    #Faeces, exp. 6&8&9&10, correlation with Histo_Total
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.post.histo";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "faeces" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Histo_Total";
    t <- t + 1;
    #Faeces, exp. 6&8&9&10, correlation with Eosino
    collect.tests[[t]] <- list();
    collect.tests[[t]]$type <- "correlation";
    collect.tests[[t]]$name <- tmp.colnames[t] <- "faeces.post.eosino";
    collect.tests[[t]]$curr.subset <- rownames(sample.data)[which(sample.data$Material == "faeces" & sample.data$Time == "post")];
    collect.tests[[t]]$cp <- "Eosino";
    t <- t + 1;
  }
  colnames(collect.data) <- tmp.colnames;
  ################################
  #Iterate through tests and collect data
  for (t in seq(1, length(collect.tests))) {
    attach(collect.tests[[t]]);
    if (type == "association") {
      #Subset current data
      curr.physeq <- prune_samples(curr.subset, pre_filt.physeq);
      #Coerce into DGEList object
      x <- phyloseq_to_edgeR(curr.physeq, group=group_by);
      #Perform exact test
      et <- exactTest(x, pair=c(baseline, response));
      #Extract results
      tt <- topTags(et, n=nrow(x), adjust.method="BH", sort.by="PValue");
      #Process all topTags
      collect.tests[[t]]$res <- tt$table$logFC;
      collect.data[, t] <- tt$table[rownames(collect.data), "logFC"];
      candidate.otus <- rownames(tt$table)[tt$table$FDR < PARAM$alpha.cutoff];
      raw.collect.otus <- c(raw.collect.otus, candidate.otus);
      if (length(candidate.otus) > 20) {
        collect.otus <- c(collect.otus, candidate.otus[1:20]);
        writeLines(paste("Dropped", as.character(length(candidate.otus) - 20), "significant OTUs in", name));
      } else {
        collect.otus <- c(collect.otus, candidate.otus);
      }
    } else {
      #Subset current data
      curr.samples <- which(rownames(sample.data) %in% curr.subset);
      curr.ot <- norm.ot[, curr.subset];
      curr.otus <- rownames(curr.ot);
      #Compute correlations
      tmp.cor <- unlist(mclapply(rownames(curr.ot), function(otu) {cor(curr.ot[otu,], sample.data[[cp]][curr.samples], method="spearman", use="complete.obs")}, mc.cores=PARAM$use.cores));
      names(tmp.cor) <- rownames(curr.ot);
      #Process and store
      collect.tests[[t]]$res <- tmp.cor;
      collect.data[, t] <- tmp.cor;
      ordered.cor <- sort(abs(tmp.cor), decreasing=T);
      candidate.otus <- names(ordered.cor)[ordered.cor > PARAM$correlation.cutoff[2]];
      raw.collect.otus <- c(raw.collect.otus, candidate.otus);
      if (length(candidate.otus) > 20) {
        collect.otus <- c(collect.otus, candidate.otus[1:20]);
        writeLines(paste("Dropped", as.character(length(candidate.otus) - 20), "significant OTUs in", name));
      } else {
        collect.otus <- c(collect.otus, candidate.otus);
      }
    }
    #Tidy up
    detach(collect.tests[[t]]);
  }
  #Be sure to detach, even if loop crashes for some reason...
  if ("collect.tests[[t]]" %in% search()) {detach(collect.tests[[t]])}
  ################################
  #Summarize & sort results
  unique.otus <- unique(collect.otus);
  n.sig_hits <- length(unique.otus);
  sig.data <- collect.data[unique.otus, ];
  sig.otu_data <- otu.data[unique.otus, ];
  sig.size <- data$otu$size[unique.otus] / sum(data$otu$size);
  sig.clust <- hclust(Dist(sig.data[,1:8], method="pearson"));
  sig.order <- sig.clust$order;
  sig.export <- data.frame(Name=rownames(sig.data)[sig.order], sig.otu_data[sig.order, ], Size=sig.size[sig.order], sig.data[sig.order, ]);
  sig.export$Name <- factor(sig.export$Name, levels=rownames(sig.data)[sig.order]);
  write.table(sig.export, file=paste(PARAM$folder$results, paste("significant_otus", c, "tsv", sep="."), sep=""), sep="\t", col.names=NA);
  #Plot heatmaps (which will afterwards be combined into one figure)
  #=> taxonomy at order level
  curr.melt <- melt(sig.export, measure.vars="order");
  pdf(file=paste(PARAM$folder$results, paste("heatmap.taxonomy", c, "pdf", sep="."), sep=""), width=4, height=8);
  curr.plot <- ggplot(curr.melt, aes(x=variable, y=Name, fill=value)) +
    geom_tile() +
    scale_x_discrete(breaks=NULL) + 
    theme(axis.text.y=element_text(size=6), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    coord_fixed(ratio=0.3);
  print(curr.plot);
  dev.off();
  #=> OTU sizes
  curr.melt <- melt(sig.export, measure.vars="Size");
  pdf(file=paste(PARAM$folder$results, paste("heatmap.sizes", c, "pdf", sep="."), sep=""), width=4, height=8);
  curr.plot <- ggplot(curr.melt, aes(x=variable, y=Name, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(limits=c(10^-5, 2*10^-3), midpoint=min(curr.melt$value), high = "black") +
    scale_x_discrete(breaks=NULL) + 
    theme(axis.text.y=element_text(size=6), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    coord_fixed(ratio=0.3);
  print(curr.plot);
  dev.off();
  #=> OTU associations (fold changes)
  curr.melt <- melt(sig.export, measure.vars=tmp.colnames[1:8]);
  pdf(file=paste(PARAM$folder$results, paste("heatmap.fold_changes", c, "pdf", sep="."), sep=""), width=8, height=8);
  curr.plot <- ggplot(curr.melt, aes(x=variable, y=Name, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(breaks=c(-6,-4,-2,0,2,4,6), colours=rev(brewer.pal(11, "RdBu"))) +
    scale_x_discrete(breaks=NULL) +
    theme(axis.text.y=element_text(size=6), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    coord_fixed(ratio=0.3);
  print(curr.plot);
  dev.off();
  #=> OTU correlations
  curr.melt <- melt(sig.export, measure.vars=tmp.colnames[9:16]);
  pdf(file=paste(PARAM$folder$results, paste("heatmap.correlations", c, "pdf", sep="."), sep=""), width=8, height=8);
  curr.plot <- ggplot(curr.melt, aes(x=variable, y=Name, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(breaks=c(-1,-0.5,0,0.5,1), colours=rev(brewer.pal(11, "PiYG"))) +
    scale_x_discrete(breaks=NULL) +
    theme(axis.text.y=element_text(size=6), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    coord_fixed(ratio=0.3);
  print(curr.plot);
  dev.off();
}









############ PRELIM CODE ###############

#Must be run using the latest version of R (ver. 4.3.2)
#installing packages
#install.packages("reshape2")
library(reshape2)

#install.packages("ggplot2")
library(ggplot2)

#install.packages("ggrepel")
library(ggrepel)

#install.packages("amap")
library(amap)

#install.packages("ggridges")
library(ggridges)

#BiocManager::install("STRINGdb")
library(STRINGdb)

#BiocManager::install("clusterProfiler")
#samples are human so downloading the human database
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

#install.packages("devEMF")
library(devEMF)


#setting working directory filepath  
setwd("C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")

#loading files
em = read.table("EM.csv ", header=TRUE, row.names=1, sep= "\t")
de_mvp = read.table("DE_Senes_MtD_vs_prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mvs = read.table("DE_Senes_Mtd_vs_Senes.csv", header=TRUE, row.names=1, sep= "\t")
de_svp = read.table("DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
sample_sheet = read.table("sample_sheet.csv", header=TRUE, row.names=1, sep= "\t")
human_background = read.table("Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")


######## MY THEME #########

my_theme = theme(
  plot.title = element_text(size=25),
  axis.text.x = element_text(size=15),
  axis.text.y = element_text(size=15),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  legend.text = element_text(size=15),
  legend.position = "top"
)


########### DENSITY PLOT FUNCTION ############

density_exp_plot = function(expression_table, columns, file_path)
{
  em.m = melt(expression_table)
  ggp = ggplot(em.m, aes(x = log10(value+0.01), fill=variable)) +
  geom_density(alpha =0.75) + 
  labs(title="Density Plot", x="Log10 Expression", y="Density") + 
  facet_wrap(~variable, ncol=columns) + my_theme +
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"),
        strip.text = element_text(size = 8)) +
  #below code groups the replicates into sets of 3 and assigns a color for each group 
  scale_fill_manual(values = rep(c("blue3", "darkorange1", "cornflowerblue"), each = 3))
  
  emf(file.path(file_path, "density_plot.emf"), height=4, width=8)
  print(ggp)
  dev.off()

}


####### FUNCTION TO PARSE DE TABLES #######
format_de_table = function(de_file_path)
{
  #reading in the file 
  de = read.table(de_file_path, header=TRUE, row.names=1, sep= "\t")
  #merging the DE table with the EM_annotation table made earlier 
  de_em_annotated = merge(de,em_annotation, by.x=0, by.y=1)
  #omitting any NA values 
  de_em_annotated = na.omit(de_em_annotated)
  #ordering the table based on p.adj values 
  sorted_order = order(de_em_annotated[,"p.adj"], decreasing=FALSE)
  de_em_annotated = de_em_annotated[sorted_order,]
  #column of mean expression values 
  means = rowMeans(de_em_annotated[,c(2:10)])
  de_em_annotated$mean_expression = means
  #new column: flagging significance 
  de_em_annotated$sig = as.factor(de_em_annotated$p.adj < 0.05 & abs(de_em_annotated$log2fold) > 1.0)
  #new column: log10 of p.adj
  de_em_annotated$mlog10padj = -log10(de_em_annotated$p.adj)
  #making the bold row the gene SYMBOLS - continuity across all created tables 
  rownames(de_em_annotated) = de_em_annotated[,"SYMBOL"]
  #columns_for_removal = 11:19
  #de_em_annotated = de_em_annotated[,-columns_for_removal]
  return(de_em_annotated)
}



######## FUNCTION FOR SIG DE TABLE #########

de_sig = function(de_master_table, p_threshold, fold_threshold)
{
  master_sig = subset(de_master_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)
  return(master_sig)
}


############ PCA PLOT FUNCTION ###############

# PCA PLOT
plot_general_pca = function(em_scaled, groups, variables, table_name, file_path)
{
  em_matrix = as.matrix(sapply(em_scaled,as.numeric)) 
  #transpose matrix 
  pca = prcomp(t(em_matrix))
  #extract component data?? 
  pca_coordinates = data.frame(pca$x)
  #percent values 
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  ggp = ggplot(pca_coordinates, aes(x = PC1, y = PC2, colour = groups)) + geom_point() + 
    scale_color_manual(values=c("blue3", "darkorange1", "cornflowerblue"), name='Sample Group') + 
    geom_label_repel(aes(label = variables), show.legend=FALSE, position=position_jitter(width=1,height=1)) + xlab(x_axis_label) + ylab(y_axis_label) + 
    my_theme + labs(title= table_name)
  
  emf(file_path, height=6, width=8)
  print(ggp)
  dev.off()
}


plot_pca = function(em_scaled, groups, variables, plot_title, file_path)
{
  em_matrix = as.matrix(sapply(em_scaled,as.numeric)) 
  #transpose matrix 
  pca = prcomp(t(em_matrix))
  #extract component data?? 
  pca_coordinates = data.frame(pca$x)
  #percent values 
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  ggp = ggplot(pca_coordinates, aes(x = PC1, y = PC2, colour = groups)) + geom_point() + 
    scale_color_manual(values=c("blue3", "darkorange1"), name='Sample Type') + 
    geom_text_repel(aes(label = variables), show.legend=FALSE) + xlab(x_axis_label) + ylab(y_axis_label) + 
    my_theme + labs(title= plot_title)
  
  emf(file_path, height=6, width=8)
  print(ggp)
  dev.off()
}


############ VOLCANO PLOT FUNCTION ############

plot_volcano = function(de_table, p_threshold, fold_threshold, plot_title, file_name)
{
  de_table$mlog10p = -log10(de_table$p)
  sort_order = order(de_table[,"p"], decreasing=FALSE)
  de_table = de_table[sort_order,]
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)
  de_table_sig_up = subset(de_table_sig, p.adj < 0.05 & log2fold > 1)
  de_table_sig_down = subset(de_table_sig, p.adj < 0.05 & log2fold < -1)
  de_table_sig_up_top5 = de_table_sig_up[1:5,]
  de_table_sig_down_top5 = de_table_sig_down[1:5,]
  
  ggp = ggplot(de_table, aes(x=log2fold, y=mlog10p, colour=sig)) + 
    geom_point(aes(colour="a"))+
    geom_point(data= de_table_sig_down, aes(colour="b")) + 
    geom_point(data= de_table_sig_up, aes(colour="c")) +
    labs(title= plot_title, x="Log2Fold Change", y="-Log10 p") + 
    theme_bw() + geom_vline(xintercept= -1) + 
    geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) + 
    geom_vline(xintercept=-1, linetype="dashed", color = "grey", size=0.5) +
    geom_hline(yintercept=-log10(0.05),linetype="dashed", color = "grey", size=0.5) +
    geom_label_repel(data=de_table_sig_up_top5, aes(label=rownames(de_table_sig_up_top5), colour='c'), show.legend=FALSE, max.overlaps = Inf, position=position_jitter(width=1,height=1)) + 
    geom_label_repel(data=de_table_sig_down_top5, aes(label=rownames(de_table_sig_down_top5), colour='b'), show.legend=FALSE, max.overlaps = Inf, position=position_jitter(width=1,height=1)) +
    scale_colour_manual(values = c("black", "blue", "red"), labels=c("No change","Down","Up"), name='Gene Significance') + 
    my_theme
  
  emf(file_name, height=6, width=8)
  print(ggp)
  dev.off()
}


############# MA PLOT FUNCTION ################

plot_MA = function(table, p_threshold, fold_threshold, plot_title, file_name)
{
  sort_order = order(table[,"p"], decreasing=FALSE)
  table = table[sort_order,]
  table_sig = subset(table, p.adj < p_threshold & abs(log2fold) > fold_threshold)
  table_sig_up = subset(table_sig, p.adj < 0.05 & log2fold > 1)
  table_sig_down = subset(table_sig, p.adj < 0.05 & log2fold < -1)
  
  ggp = ggplot(table, aes(x=log10(mean_expression), y=log2fold)) + 
    geom_point(aes(colour="a")) +
    geom_point(data= table_sig_down, aes(colour="b")) + 
    geom_point(data= table_sig_up, aes(colour="c")) +
    labs(title= plot_title, x="Mean Expression (log10)", y="Log2Fold Change") + theme_bw() + 
    geom_hline(yintercept=1, linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") +
    geom_hline(yintercept=-1, linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") +
    scale_colour_manual(values = c("black", "blue", "red"), labels=c("No change","Down","Up"), name= "Gene Significance") + 
    my_theme
  
  emf(file_name, height=6, width=8)
  print(ggp)
  dev.off()
}



############### BOX PLOT FUNCTIONS ###############
#box plot for one gene
box_plot_one_gene = function(gene_name, em_symbols, ss, plot_title, file_name)
{
  gene_data = em_symbols[gene_name,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD"))
  names(gene_data) = c("expression","sample_group")
  
  
  ggp_box = ggplot(gene_data,aes(x=sample_group,y=expression,fill=sample_group)) +
  geom_boxplot(size = 1, outlier.size = 0, alpha = 0.5) + 
  geom_jitter() +
  labs(title= plot_title, x="Sample Group", y="Expression") + theme_bw() + 
  scale_fill_manual(values = c("blue3", "darkorange2", "cornflowerblue"),name= "Sample Type") + 
  my_theme

  emf(file_name, height=6, width=8)
  print(ggp_box)
  dev.off()
}



#box plot for a vector of several genes
facet_box_plot = function(genes,em_scaled, ss, plot_name, file_name)
{
  gene_data = em_scaled[genes,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  gene_data.m = melt(gene_data, id.vars="sample_group")
  
  ggp_box = ggplot(gene_data.m,aes(x=sample_group,y=value,fill = sample_group, color = sample_group)) +
    geom_boxplot(size = 1, outlier.size = 0, alpha = 0.5) + geom_jitter(show.legend=FALSE) +
    facet_wrap(~variable,ncol=9) + 
    labs(title= plot_name, y="Expression") + theme_bw() + 
    scale_fill_manual(values = c("blue3", "darkorange2", "cornflowerblue"),name= "Sample Type") +
    scale_color_manual(values = c("blue3", "darkorange1", "cornflowerblue")) +  my_theme + 
  theme(strip.background=element_rect(colour="black",fill="white"), strip.text = element_text(size = 10), axis.text.x = element_blank()) + 
    guides(color = guide_legend(title = "Sample Group"), fill = "none") 
  
  emf(file_name, height=6, width=8)
  print(ggp_box)
  dev.off()
  
}



########### MULTI GENE BOX PLOT FUNCTION ###########
#for two sample groups 
#orientation indicates the values that go on the x and y axis. "normal" means x axis has the variable (i.e. gene name)
#flipped means the gene names are on the y axis. this helps for orientation for presentation 
multigene_boxplot <- function(em_scaled, gene_list, ss, plot_title, file_name)
{
  gene_data = em_scaled[gene_list,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  
  # Melt the table
  gene_data.m = melt(gene_data, id.vars = "sample_group")
  
  # Make the boxplot
  ggp = ggplot(gene_data.m, aes(x = variable, y = value, fill = sample_group, color = sample_group)) +
    geom_boxplot(size = 1, outlier.size = 0, alpha = 0.5) +
    labs(title = plot_title, x = "Genes", y = "Expression") +
    scale_fill_manual(values = c("blue3", "darkorange1", "cornflowerblue"), name = "Sample Group") +
    scale_color_manual(values = c("blue3", "darkorange1", "cornflowerblue")) +  theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.2)) + 
    guides(color = guide_legend(title = "Sample Group"), fill = "none") + 
    my_theme 
    
  
  emf(file_name, height=6, width=15)
  print(ggp)
  dev.off()
}


############### HEAT MAP FUNCTION ###############

#NOTE: THIS DOES NOT INCLUDE THE RUG YET OR THE X AXIS CLUSTERED
heat_map = function(em_scaled_sig_table, y_cluster = TRUE, x_cluster = TRUE, rug = TRUE, file_name)
{
  #clustering y axis
  if (y_cluster == TRUE) {
  hm.matrix = as.matrix(em_scaled_sig_table)
  y.Dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.Dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order=order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  }
  #clustering x axis 
  if (x_cluster == TRUE) {
  hm.matrix = as.matrix(em_scaled_sig_table)
  x.Dist = Dist(hm.matrix, method="spearman")
  x.cluster = hclust(y.Dist, method="average")
  x.dd = as.dendrogram(y.cluster)
  x.dd.reorder = reorder(y.dd,0,FUN="average")
  x.order=order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[x.order,]
  }
  
  hm.matrix = melt(hm.matrix_clustered)
  
  colours = c("blue","azure2", "red")
  
  ggp = ggplot(hm.matrix, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100)) + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(),legend.spacing.x = unit(0.25, 'cm'))

  if (rug == TRUE){
    
    groups_data = as.matrix(as.numeric(as.factor(ss$SAMPLE_GROUP)))
    groups_data = melt(groups_data)
    
    ggp = ggp + geom_rug(groups_data, aes(x = Var2, y = Var1, fill = value)) + 
      geom_tile() + scale_fill_gradientn(colours = colours) + geom_tile(linetype="blank") + 
      labs(x = "", y = "") + 
      theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y =element_blank(), axis.ticks=element_blank())
  }
  
  emf(file_name, height=6, width=8)
  print(ggp)
  dev.off()
}


########### ORA AND GSE FUNCTIONS #############

#for all sig genes, as well as upregulated and downregulated genes 
#desired ora analysis will determine if the test is run for all sig, up sig, or down siora_results = function(data, desired_ora_analysis, folder_name, file_path)
ora_results = function(data, desired_ora_analysis, folder_name, file_path)
{
  
  #creating the folder where the created plots will go 
  output_folder = file.path(file_path,folder_name)
  dir.create(output_folder)
  
  if (desired_ora_analysis == "all"){
  
  #converting the IDs
    sig_genes_entrez_all = bitr(data, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
    #running the ORA analysis based on the desired test 
    ora_results_all = enrichGO(gene = sig_genes_entrez_all$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.001, qvalueCutoff = 0.10)
    
    #creating ORA plots for data visualization 
    #results for all sig genes 
    ora_barplot_all = barplot(ora_results_all, showCategory=10)
    ora_dotplot_all = dotplot(ora_results_all, showCategory=10)
    ora_goplot_all = goplot(ora_results_all, showCategory = 10)
    ora_cnetplot_all = cnetplot(ora_results_all, categorySize= "pvalue")
    
    #pulling the gene-set, p and description information from the pathway analysis result for all sig genes 
    gene_sets = ora_results_all$geneID
    description = ora_results_all$Description
    p.adj = ora_results_all$p.adjust
    #data frame 
    ora_results_all_table = data.frame(cbind(description, gene_sets, p.adj))
    #list of the most enriched genes 
    enriched_gene_set = as.character(ora_results_all_table[1,2])
    candidate_genes_all = unlist(strsplit(enriched_gene_set, "/"))
    
    #prep to run string
    candidate_genes_table = data.frame(candidate_genes_all)
    names(candidate_genes_table) = "gene"
  
    # loading the database
    string_db = STRINGdb$new(version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")
    # mapping genes to the database
    string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
    
    #plotting results 
    emf(file.path(output_folder, "string_results_all_ora.emf"))
    print(string_db$plot_network(string_mapped))
    dev.off()
    
    #creating a boxplot 
    genes_for_plot = candidate_genes_all[1:20]
    multigene_boxplot(em_scaled, genes_for_plot, sample_sheet, "Most Enriched Pathway Gene Set", file.path(output_folder, "Boxplot_Enriched_all_genes.emf"))
    #creating a heatmap
    em_candidate = em_scaled[candidate_genes_all,]
    heat_map(em_candidate, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, file.path(output_folder, "Enriched_genes_all_heatmap.emf"))
  
  
    emf(file.path(output_folder, "ora_barplot_all.emf"))
    print(ora_barplot_all)
    dev.off()
    
    emf(file.path(output_folder, "ora_dotplot_all.emf"))
    print(ora_dotplot_all)
    dev.off()
    
    emf(file.path(output_folder, "ora_goplot_all.emf"))
    print(ora_goplot_all)
    dev.off()
    
    emf(file.path(output_folder, "ora_cnetplot_all.emf"))
    print(ora_cnetplot_all)
    dev.off()
  
  }
  
  if (desired_ora_analysis == "up"){
    
  sig_genes_entrez_up = bitr(data, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results_up = enrichGO(gene = sig_genes_entrez_up$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.001, qvalueCutoff = 0.10)
  
  #results for up regulated genes
  ora_barplot_up = barplot(ora_results_up, showCategory=10)
  ora_dotplot_up = dotplot(ora_results_up, showCategory=10)
  ora_goplot_up = goplot(ora_results_up, showCategory = 10)
  ora_cnetplot_up = cnetplot(ora_results_up, categorySize= "pvalue")
  
  #pulling the gene-set, p and description information from the pathway analysis result for all sig genes 
  gene_sets = ora_results_up$geneID
  description = ora_results_up$Description
  p.adj = ora_results_up$p.adjust
  #data frame 
  ora_results_up_table = data.frame(cbind(description, gene_sets, p.adj))
  #list of the most enriched genes 
  enriched_gene_set = as.character(ora_results_up_table [1,2])
  candidate_genes_up = unlist(strsplit(enriched_gene_set, "/"))
  
  #prep to run string
  candidate_genes_table = data.frame(candidate_genes_up)
  names(candidate_genes_table) = "gene"
  
  # loading the database
  string_db = STRINGdb$new(version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")
  # mapping genes to the database
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
  
  #plotting results 
  emf(file.path(output_folder, "string_results_up_ora.emf"))
  print(string_db$plot_network(string_mapped))
  dev.off()
  
  #creating a boxplot 
  multigene_boxplot(em_scaled, candidate_genes_up, sample_sheet, "Most Enriched Pathway - Upregulated Genes", file.path(output_folder, "Boxplot_Enriched_up_genes.emf"))
  #creating a heatmap
  em_candidate = em_scaled[candidate_genes_up,]
  heat_map(em_candidate, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, file.path(output_folder,"Enriched_genes_up_heatmap.emf"))
  

  emf(file.path(output_folder, "ora_barplot_up.emf"))
  print(ora_barplot_up)
  dev.off()
  
  emf(file.path(output_folder, "ora_dotplot_up.emf"))
  print(ora_dotplot_up)
  dev.off()
  
  emf(file.path(output_folder, "ora_goplot_up.emf"))
  print(ora_goplot_up)
  dev.off()
  
  emf(file.path(output_folder, "ora_cnetplot_up.emf"))
  print(ora_cnetplot_up)
  dev.off()
  }
  
  if (desired_ora_analysis == "down"){
  sig_genes_entrez_down = bitr(data, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results_down = enrichGO(gene = sig_genes_entrez_down$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.001, qvalueCutoff = 0.10)

  #results for down regulated genes 
  ora_barplot_down = barplot(ora_results_down, showCategory=10)
  ora_dotplot_down = dotplot(ora_results_down, showCategory=10)
  ora_goplot_down = goplot(ora_results_down, showCategory = 10)
  ora_cnetplot_down = cnetplot(ora_results_down, categorySize= "pvalue")
  
  #pulling the gene-set, p and description information from the pathway analysis result for all sig genes 
  gene_sets = ora_results_down$geneID
  description = ora_results_down$Description
  p.adj = ora_results_down$p.adjust
  #data frame 
  ora_results_down_table = data.frame(cbind(description, gene_sets, p.adj))
  #list of the most enriched genes 
  enriched_gene_set = as.character(ora_results_down_table [1,2])
  candidate_genes_down = unlist(strsplit(enriched_gene_set, "/"))
  
  #prep to run string
  candidate_genes_table = data.frame(candidate_genes_down)
  names(candidate_genes_table) = "gene"
  
  # loading the database
  string_db = STRINGdb$new(version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")
  # mapping genes to the database
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
  
  #plotting results 
  emf(file.path(output_folder, "string_results_down_ora.emf"))
  print(string_db$plot_network(string_mapped))
  dev.off()
  
  #creating a boxplot 
  multigene_boxplot(em_scaled, candidate_genes_down, sample_sheet, "Enriched Gene Set", file.path(output_folder,"Boxplot_Enriched_genes_down.emf"))
  #creating a heatmap
  em_candidate = em_scaled[candidate_genes_down,]
  heat_map(em_candidate, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, file.path(output_folder,"Enriched_genes_down_heatmap.emf"))
  

  emf(file.path(output_folder, "ora_barplot_down.emf"))
  print(ora_barplot_down)
  dev.off()
  
  emf(file.path(output_folder, "ora_dotplot_down.emf"))
  print(ora_dotplot_down)
  dev.off()
  
  emf(file.path(output_folder, "ora_goplot_down.emf"))
  print(ora_goplot_down)
  dev.off()
  
  emf(file.path(output_folder, "ora_cnetplot_down.emf"))
  print(ora_cnetplot_down)
  dev.off()
  }
}


gse_results = function(table, em_table, folder_name, file_path)
{
  #grabbing the log2fold from the input table 
  gsea_input = table$log2fold
  #adding gene names 
  names(gsea_input) = row.names(table)
  #omitting any NA values 
  gsea_input = na.omit(gsea_input)
  #sort in decreasing order 
  gsea_input = sort(gsea_input, decreasing = TRUE)
  #running the GSE analysis 
  gse_results = gseGO(geneList=gsea_input,ont ="BP",keyType = "SYMBOL",nPerm = 10000,minGSSize = 3,
  maxGSSize = 800,pvalueCutoff = 0.05,verbose = TRUE,OrgDb = org.Hs.eg.db,pAdjustMethod = "none")

  output_folder = file.path(file_path,folder_name)
  dir.create(output_folder)
  
  #plot results
  ggp = ridgeplot(gse_results)
  #export as emf file 
  emf(file.path(output_folder, "Results_GSE.emf"))
  print(ggp)
  dev.off()
  
  #pulling info from the gse results 
  gene_description = gse_results$Description
  p.adj = gse_results$p.adjust
  nes_score = gse_results$NES
  genes = gse_results$core_enrichment
  
  #data frame 
  gse_results_table = data.frame(cbind(gene_description, genes, p.adj, nes_score))
  #list of the most enriched genes 
  enriched_gene_set = as.character(gse_results_table [c(1:10),2])
  candidate_genes_gse = unlist(strsplit(enriched_gene_set, "/"))
  
  #creating a boxplot 
  multigene_boxplot(em_table, table, sample_groups, "Top Most Enriched Gene Set", file.path(output_folder,"GSE_Enriched_genes_boxplot.emf"))
  
  #creating a heatmap
  em_candidate = em_table[table,]
  heat_map(em_candidate, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, file.path(output_folder,"GSE_Enriched_genes_heatmap.emf"))
  
  candidate_genes_table = data.frame(table)
  names(candidate_genes_table) = "gene"
  
  # loading the database
  string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")
  # mapping genes to the database
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
  #plotting results 
  emf(file.path(output_folder, "string_results_gse.emf"))
  print(string_db$plot_network(string_mapped))
  dev.off()
  
}


######## META GENE BOX PLOT FUNCTION ###########
# This is used for the sig gene data 
metagene_boxplot = function(em_signature, ss, plot_title, file_name)
{
  metagene = data.frame(colMeans(em_signature))
  metagene$sample_group = ss$SAMPLE_GROUP
  
  # Melt the table
  metagene.m = melt(metagene, id.vars = "sample_group")
  
  # Make the boxplot
  ggp = ggplot(metagene.m, aes(x = variable, y = value, fill = sample_group, color = sample_group)) +
    geom_boxplot(size = 1, outlier.size = 0, alpha = 0.5) +
    labs(title = plot_title, x = "Genes", y = "Expression") + theme_bw() + my_theme + 
    scale_fill_manual(values = c("blue3", "darkorange1", "cornflowerblue"), name = "Sample Group") +
    scale_color_manual(values = c("blue3", "darkorange1", "cornflowerblue")) +
    guides(fill = guide_legend(title = "Sample Group"), color = "none") + theme(strip.text.x=element_blank(), axis.text.x = element_blank())

  
  emf(file_name, height=6, width=8)
  print(ggp)
  dev.off()
}



################### MAIN CODE #####################


############ FORMATTING EM DATA
#creating a master table (particularly interested in joining em with the symbols column from the human_background table)
#merging based on gene IDs
em_annotation = merge(em, human_background, by.x=0, by.y=0)
#fix row names
row.names(em_annotation) = em_annotation[,"SYMBOL"]
#fixing the name of the column for gene ID
names(em_annotation)[1] = "GENE_ID"

#edit the sample sheet table so the below code is easier (for me)
sample_sheet$SAMPLE = rownames(sample_sheet)
#make em_symbols table using sample sheet table 
em_symbols = em_annotation[ , as.vector(sample_sheet$SAMPLE)]

#calling the density function 
ggp_density = density_exp_plot(em_symbols, 9, "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ggp_density 
#based on the results of the density plots, the data appears to be generally good quality


#creating a master table for each DE test using the format de table function 
master_mvp = format_de_table("DE_Senes_MtD_vs_prolif.csv")
master_mvs = format_de_table("DE_Senes_Mtd_vs_Senes.csv")
master_svp = format_de_table("DE_Senes_vs_Prolif.csv")


#making a list of the significant genes
sig_master_mvp = de_sig(master_mvp, 0.001, 2)
sig_master_mvs = de_sig(master_mvs, 0.001, 2)
sig_master_svp = de_sig(master_svp, 0.001, 2)


######## MORE DATA PREP 
#making vectors of the significant genes 
sig_genes_mvp = rownames(sig_master_mvp)
df_sig_master_mvp = data.frame(SYMBOL = sig_genes_mvp, stringsAsFactors = FALSE)

sig_genes_mvs = rownames(sig_master_mvs)
df_sig_master_mvs = data.frame(SYMBOL = sig_genes_mvs, stringsAsFactors = FALSE)

sig_genes_svp = rownames(sig_master_svp)
df_sig_master_svp = data.frame(SYMBOL = sig_genes_svp, stringsAsFactors = FALSE)

#now creating an expression table of just the significant genes
em_symbols_sig_mvp = em_symbols[sig_genes_mvp,]
em_symbols_sig_mvs = em_symbols[sig_genes_mvp,]
em_symbols_sig_svp = em_symbols[sig_genes_svp,]

#prepping to make PCA plots 
#scaled data is important for PCA plots(and heatmaps)
em_scaled = na.omit(data.frame(t(scale(t(em_symbols)))))

em_scaled_mvp = subset(em_scaled, select = -c(4:6))
em_scaled_mvs = subset(em_scaled, select = -c(1:3))
em_scaled_svp = subset(em_scaled, select = -c(7:9))

#writing tables for scaled em tables with just sig genes 
em_scaled_sig_mvp = em_scaled[sig_genes_mvp, ]
em_scaled_sig_mvs = em_scaled[sig_genes_mvs, ]
em_scaled_sig_svp = em_scaled[sig_genes_svp, ]

em_scaled_sig_mvp = em_scaled_sig_mvp[,-c(4:6)]
em_scaled_sig_mvs = em_scaled_sig_mvs[,-c(1:3)]
em_scaled_sig_svp = em_scaled_sig_svp[,-c(7:9)]

############ CREATING FULL MASTER TABLES ############
de_mvp$sig = factor(de_mvp$p.adj < 0.001 & abs(de_mvp$log2fold) > 2)
de_mvs$sig = factor(de_mvs$p.adj < 0.001 & abs(de_mvs$log2fold) > 2)
de_svp$sig = factor(de_svp$p.adj < 0.001 & abs(de_svp$log2fold) > 2)

names(de_svp) = c("log2fold.svp", "p.svp", "p.adj.svp", "sig.svp")
#merging the DE tables 
master = merge(de_mvp, de_mvs, by.x=0, by.y=0, suffixes=c(".mvp",".mvs"))
master = merge(master, de_svp, by.x=1, by.y=0)

#fixing the name of the column for gene ID
names(master)[1] = "GENE_ID"

master = merge(master, em_annotation, by.x=1, by.y=1)
row.names(master) = master[,"SYMBOL"]
#now have a full master table of all the DE test results 

#general PCA plot comparing all the expression values between the 9 replicates
ggp_pca = plot_general_pca(em_scaled, sample_sheet$SAMPLE_GROUP, sample_sheet$SAMPLE, "Senes MtD vs Prolif vs Senes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PCA_plot.emf")
ggp_pca
#Conclusion: there is VERY distinct clustering of the expression values for each group...almost suspiciously distinct


#isolating the relevant groups for the PCA plots 
mvp_groups = sample_sheet[-c(4:6),]
mvs_groups = sample_sheet[-c(1:3),]
svp_groups = sample_sheet[-c(7:9),]

#making a PCA plot for comparing each group 
ggp_pca2 = plot_pca(em_scaled_mvp, mvp_groups$SAMPLE_GROUP, mvp_groups$SAMPLE, "Senes MtD vs Prolif", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PCA_plot_mvp1.emf")
ggp_pca2
ggp_pca3 = plot_pca(em_scaled_mvs, mvs_groups$SAMPLE_GROUP, mvs_groups$SAMPLE, "Senes MtD vs Senes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PCA_plot_mvs.emf")
ggp_pca3
ggp_pca4 = plot_pca(em_scaled_svp, svp_groups$SAMPLE_GROUP, svp_groups$SAMPLE, "Senes vs Prolif", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PCA_plot_svp.emf")
ggp_pca4


#this will fix the infinite values in my plot for p 
smallest_value_mvp = min(master_mvp$p[master_mvp$p > 0])
master_mvp$p[master_mvp$p == 0] = smallest_value_mvp

smallest_value_mvs = min(master_mvs$p[master_mvs$p > 0])
master_mvs$p[master_mvs$p == 0] = smallest_value_mvs

smallest_value_svp = min(master_svp$p[master_svp$p > 0])
master_svp$p[master_svp$p == 0] = smallest_value_svp

#now for p.adj (if i wanted to do max, just replace min with max)
smallest_value = min(master_mvp$p.adj[master_mvp$p.adj > 0])
master_mvp$p.adj[master_mvp$p.adj == 0] = smallest_value

smallest_value_mvs = min(master_mvs$p.adj[master_mvs$p.adj > 0])
master_mvs$p.adj[master_mvs$p.adj == 0] = smallest_value_mvs

smallest_value_svp = min(master_svp$p.adj[master_svp$p.adj > 0])
master_svp$p.adj[master_svp$p.adj == 0] = smallest_value_svp

#creating a volcano plot for each DE table 
ggp_volcano1 = plot_volcano(master_mvp, 0.001, 2, "Senes MtD vs Prolif","C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\volcano_plot_mvp1.emf")
ggp_volcano1

ggp_volcano2 = plot_volcano(master_mvs, 0.001, 2, "Senes MtD vs Senes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\volcano_plot_mvs.emf")
ggp_volcano2

ggp_volcano3 = plot_volcano(master_svp, 0.001, 2,"Senes vs Prolif", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\volcano_plot_svp.emf")
ggp_volcano3

#creating a MA plot for each table 
ggp_MA1 = plot_MA(master_mvp, 0.001, 2, "Senes MtD vs Prolif", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\MA_plot_mvp.emf")
ggp_MA1

ggp_MA2 = plot_MA(master_mvs, 0.001, 2, "Senes MtD vs Senes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\MA_plot_mvs.emf")
ggp_MA2

ggp_MA3 = plot_MA(master_svp, 0.001, 2, "Senes vs Prolif", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\MA_plot_svp.emf")
ggp_MA3


#picking genes to view and test box plot 
box_plot_one_gene("NIPAL3", em_scaled, sample_sheet, "NIPAL3","C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\Gene_Box_Plot.emf") 
box_plot_one_gene("COL1A1", em_scaled, sample_sheet, "COL1A1", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\COL1A1_Box_Plot.emf") 
box_plot_one_gene("PGK1", em_scaled, sample_sheet, "PGK1", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PGK1_Box_Plot.emf") 
box_plot_one_gene("SLC2A3", em_scaled, sample_sheet, "SLC2A3", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\SLC2A3_Box_Plot.emf")
box_plot_one_gene("TPX2", em_scaled, sample_sheet, "TPX2", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\TPX_2_Box_Plot.emf")
#downregulated genes in mVs
box_plot_one_gene("CEMIP", em_scaled, sample_sheet, "CEMPIP", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\CEMIP_Box_Plot.emf")
box_plot_one_gene("RIPOR3", em_scaled, sample_sheet, "RIPOR3", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\RIPOR3_Box_Plot.emf")
box_plot_one_gene("MAN1C1", em_scaled, sample_sheet, "MAN1C1", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\MAN1C1_Box_Plot.emf")
box_plot_one_gene("ACTA2", em_scaled, sample_sheet, "ACTA2", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\ACTA2_Box_Plot.emf")
#upregulated genes in mvs
box_plot_one_gene("NDRG1", em_scaled, sample_sheet, "NDRG1", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\NDRG1_Box_Plot.emf")
box_plot_one_gene("ENO2", em_scaled, sample_sheet, "ENO2", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\ENO2_Box_Plot.emf")
box_plot_one_gene("PFKFB4", em_scaled, sample_sheet, "PFKFB4", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\PFKFB4_Box_Plot.emf")

#mvp: SLC2A3 (or TPX2)
box_plot_one_gene("SLC2A3", em_scaled_sig_mvp, mvp_groups, "SLC2A3", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\mvp_SLC2A3_Box_Plot.emf")
box_plot_one_gene("TPX2", em_scaled_sig_mvp, mvp_groups, "SLC2A3", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\mvp_TPX2_Box_Plot.emf")
#mvs: RIPOR3 (or SLC2A3)
box_plot_one_gene("SLC2A3", em_scaled, sample_sheet, "SLC2A3", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\SLC2A3_Box_Plot.emf")
#svp: TPX2 (or RIPOR3)
box_plot_one_gene("RIPOR3", em_scaled_sig_svp, svp_groups,"SLC2A3","C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\RIPOR3_svp_Box_Plot.emf") 


#pulling the top-most significant genes and turning into a vector for input into the facet plot function 
mvp_genes = em_scaled_sig_mvp[1:8,]
mvp_genes = rownames(mvp_genes)

mvs_genes = em_scaled_sig_mvs[1:8,]
mvs_genes = rownames(mvs_genes)

svp_genes = em_scaled_sig_svp[1:8,]
svp_genes = rownames(svp_genes)

#this is how i'd do it for one DE table 
facet_box_plot(mvp_genes, em_scaled_sig_mvp, mvp_groups,"Top Significant Differentially Expressed Genes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\mvp_facet_Plot.emf") 
#this is comparing the mvp sig genes across all three groups 
facet_box_plot(mvp_genes, em_scaled, sample_sheet,"Top Significant Differentially Expressed Genes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\mvp_facet_Plot.emf")
#mvs top sig genes
facet_box_plot(mvs_genes, em_scaled, sample_sheet,"Top Significant Differentially Expressed Genes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\mvs_facet_Plot.emf") 
#svp top sig genes 
facet_box_plot(svp_genes, em_scaled, sample_sheet,"Top Significant Differentially Expressed Genes", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\svp_facet_Plot.emf") 

#taking a subset of sig genes (to make heat map more visually appealing/less messy)
mvp_sig_subset = em_scaled_sig_mvp[1:100,]
mvs_sig_subset = em_scaled_sig_mvs[1:100,]
svp_sig_subset = em_scaled_sig_svp[1:100,]


#making heat maps for the subset of sig genes from each DE table 
heat_map_plot = heat_map(mvp_sig_subset, y_cluster = TRUE, x_cluster = TRUE, rug = FALSE, "mvp_heat.emf")
print(heat_map_plot)

heat_map_plot = heat_map(mvs_sig_subset, y_cluster = TRUE, x_cluster = TRUE, rug = FALSE, "mvs_heat.emf")
print(heat_map_plot)

heat_map_plot = heat_map(svp_sig_subset, y_cluster = TRUE, x_cluster = TRUE, rug = FALSE, "svp_heat.emf")
print(heat_map_plot)


## MAKING A HEAT MAP FOR ALL SIG GENES FROM ALL TESTS (ALL 3 GROUPS, ALL 9 REPLICATES)
sig_all = subset(master, sig.mvp == TRUE | sig.mvs == TRUE | sig.svp == TRUE)$"SYMBOL"
sig_all_em_scaled = em_scaled[sig_all,]

heat_map(sig_all_em_scaled, y_cluster = TRUE, x_cluster = TRUE, rug = FALSE, "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data\\all_sig_genes_heat_map.emf")

#pulling the sig gene information needed to carry out different ORA analyses 
sig_genes_mvp = row.names(sig_master_mvp)
sig_genes_mvs = row.names(sig_master_mvs)
sig_genes_svp = row.names(sig_master_svp)

sig_genes_mvp_up = row.names(subset(sig_master_mvp, log2fold > 0))
sig_genes_mvp_down = row.names(subset(sig_master_mvp, log2fold < 0))

sig_genes_mvs_up = row.names(subset(sig_master_mvs, log2fold > 0))
sig_genes_mvs_down = row.names(subset(sig_master_mvs, log2fold < 0))

sig_genes_svp_up = row.names(subset(sig_master_svp, log2fold > 0))
sig_genes_svp_down = row.names(subset(sig_master_svp, log2fold < 0))

#ORA results for all significant genes 
ora_results(sig_genes_mvp, "all", "ora_mvp_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_mvs, "all", "ora_mvs_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_svp, "all", "ora_svp_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")

#ORA results for just upregulated significant genes 
ora_results(sig_genes_mvp_up, "up", "ora_mvp_up", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_mvs_up, "up", "ora_mvs_up", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_svp_up, "up", "ora_svp_up", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")

#ORA results for just downregulated significant genes 
ora_results(sig_genes_mvp_down, "down", "ora_mvp_down", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_mvs_down, "down", "ora_mvs_down", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(sig_genes_svp_down, "down", "ora_svp_down", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")

#gse results for each DE table 
gse_results(master_mvp, em_scaled, "gse_mvp_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
gse_results(master_mvs, em_scaled, "gse_mvs_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
gse_results(master_svp, em_scaled, "gse_svp_results", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")


## SIGNATURES 
#pulling genes from the Master table that apply to the scenarios outlined 
#genes that are downregulated in senescent cells but upregulated in mtd senescent 
#genes to do with muscle and bone density (see mice study on BAX deletion) also genes involved in mitochondrial depletion??
signature_1 = row.names(subset(master, (sig.mvs == TRUE & log2fold.mvs > 0) &
                                 (sig.svp == TRUE & log2fold.svp < 0)))
#genes that are upregulated in senescent cells but downregulated in mtd senescent 
#genes involved in inflammation, cytokines, kinase inhibitors 
signature_2 = row.names(subset(master, (sig.mvs == TRUE & log2fold.mvs < 0) &
                                 (sig.svp == TRUE & log2fold.svp > 0)))
#genes that are upregulated in MvS and SvP
#genes likely involved in oxidative stress
#log2fold + in mvs, then that means the mtd senes was upregulated 
signature_3 = row.names(subset(master, (sig.mvs == TRUE & log2fold.mvs > 0) &
                                 (sig.svp == TRUE & log2fold.svp > 0)))
#genes that are downregulated in MvS and SvP
#PLEK2: promotes proliferation and migration of a certain type of lung cancer 
signature_4 = row.names(subset(master, (sig.mvs == TRUE & log2fold.mvs < 0) &
                                 (sig.svp == TRUE & log2fold.svp < 0)))
#genes that are upregulated in prolif cells and senescent cells 
#genes that are downregulated in prolif cells and senescent cells 

#pulling the EM_scaled data for each signature to use for making the heat maps below 
em_sig1 = em_scaled[signature_1,]
em_sig2 = em_scaled[signature_2,]
em_sig3 = em_scaled[signature_3,]
em_sig4 = em_scaled[signature_4,]

#making heat maps of the signature genes 
heat_map(em_sig1, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, "signature_1_heatmap.emf")
heat_map(em_sig2, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, "signature_2_heatmap.emf")
heat_map(em_sig3, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, "signature_3_heatmap.emf")
heat_map(em_sig4, y_cluster = TRUE, x_cluster = FALSE, rug = FALSE, "signature_4_heatmap.emf")

#ora analysis for each signature 
#ora_results(signature_1, "all", "signature_1_ora", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
ora_results(signature_2, "all", "signature_2_ora", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
#ora_results(signature_3, "all", "signature_3_ora", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")
#ora_results(signature_4, "all", "signature_4_ora", "C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOL5379-Project Data")

"NOTE: Signatures 1, 3, and 4 did not yield results in the ORA analysis. This is likely 
because both signatures had a very small gene set and the genes within each signature were quite diverse
in function, which likely made it difficult for the database to match the genes to a common pathway.
Time permitting, more signatures would be generated to pinpoint other gene sets that have the potential 
for producing ORA results"

#metagene plots for each signature 
metagene_boxplot(em_sig1, sample_sheet, "Signature 1 Metagenes", "metagene_boxplot_sig1.emf")
metagene_boxplot(em_sig2, sample_sheet, "Signature 2 Metagenes", "metagene_boxplot_sig2.emf")
metagene_boxplot(em_sig3, sample_sheet, "Signature 3 Metagenes", "metagene_boxplot_sig3.emf")
metagene_boxplot(em_sig4, sample_sheet, "Signature 4 Metagenes", "metagene_boxplot_sig4.emf")

#taking a subset of EM_sig1 genes to plot
em_sig1_subset = em_sig1[c(1:20),]
signature_1_subset = signature_1[1:20]

#box plot of all the genes in each signature 
multigene_boxplot(em_sig1_subset, signature_1_subset, sample_sheet, "Signature 1", "signature1_genes_boxplot.emf")
multigene_boxplot(em_sig2, signature_2, sample_sheet, "Signature 2", "signature2_genes_boxplot.emf")
multigene_boxplot(em_sig3, signature_3, sample_sheet, "Signature 3", "signature3_genes_boxplot.emf")
multigene_boxplot(em_sig4, signature_4, sample_sheet, "Signature 4", "signature4_genes_boxplot.emf")

########## SPEARMAN CORRELATION 
#pulling the sig genes for each DE comparison 
sig_mvp = subset(master, sig.mvp == TRUE)
sig_mvs = subset(master, sig.mvs == TRUE)
sig_svp = subset(master, sig.svp == TRUE)
#pulling the common significant genes, which is useful for layering the correlation plot to color based on significance 
sig_mvp_mvs = subset(master, sig.mvp == TRUE & sig.mvs == TRUE)
sig_mvs_svp = subset(master, sig.mvs == TRUE & sig.svp == TRUE)
sig_mvp_svp = subset(master, sig.mvp == TRUE & sig.svp == TRUE)

#Run Spearman correlation test and plot results 
#did not make into a function because it was logical to only make one correlation plot for this 
correlation = cor.test(master$log2fold.mvs, master$log2fold.svp,method = "spearman")
ggp = ggplot(master, aes(x=log2fold.mvs, y=log2fold.svp)) + theme_bw() + 
  geom_point(aes(colour="a")) + 
  geom_point(data = sig_mvs, aes(colour="b")) +
  geom_point(data= sig_svp, aes(colour="c")) + 
  geom_point(data= sig_mvs_svp, aes(x=log2fold.mvs, y=log2fold.svp, colour="d")) + my_theme + 
  scale_colour_manual(values = c("black", "blue3", "darkorange1","cornflowerblue"), 
  labels=c("not sig", "Sig SE_MtD V SE only", "Sig SE V PRO only", "sig both"), name="Gene Significance") +
  labs(title= paste("Correlation Coefficient:",correlation$estimate), x="Senes MtD vs Senes log2fold", y="Senes vs Prolif log2fold")

emf("test_correlation_mvs_svp.emf", height=6, width=9)
print(ggp)
dev.off()


############ SAVING TABLES FUNCTION ############
write_to_file = function(table, file_path)
{
  write.table(table, file= file_path, sep="\t")
}

write_to_file(em_symbols, "em_symbols.csv")
write_to_file(em_scaled, "em_scaled.csv")

write_to_file(master_mvp, "master_de_mtd_vs_prolif.csv")
write_to_file(master_mvs, "master_de_mtd_vs_senes.csv")
write_to_file(master_svp, "master_de_senes_vs_prolif.csv")



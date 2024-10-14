#install packages
install.packages("ggplot2")
library(ggplot2)

install.packages("factoextra")
library("factoextra")

#set working directory 
setwd('C:\\Users\\2915700s\\OneDrive\\UoG Masters\\BIOl5371\\Project')

#load data into R 
data_annotate = read.table('Annotations.csv', header=TRUE, row.names = 1, sep="\t")
DE_GOUTvsHC = read.table('De_Gout_vs_HC.csv', header=TRUE, row.names = 1, sep='\t')
DE_SAvsHC = read.table('De_SA_vs_HC.csv', header=TRUE, row.names = 1, sep='\t')
Expression = read.table('Expression_Table.csv', header=TRUE, row.names = 1, sep='\t')
sample_info = read.table('Sample_Information.csv', header=TRUE, row.names = 1, sep='\t')


############################################################################################################################################################


#########Task 1 - assessment of clinical information############


#change the category of variables for clinical information 
str(sample_info)
sample_info$SEX=as.factor(sample_info$SEX)
sample_info$GROUP=as.factor(sample_info$GROUP)

#summary statistics for clinical information 
summary(sample_info)
IQR = 66.75 - 49.75
#response to this information catalogued in report 

healthy_info = subset(sample_info, GROUP == 'HEALTHY')
GOUT_info = subset(sample_info, GROUP == 'GOUT')
SA_info = subset(sample_info, GROUP == 'SA')

summary(healthy_info)
summary(GOUT_info)
summary(SA_info)


##to compare p-values between the groups and clinical information## 

ggp = ggplot(sample_info, aes(x=AGE, fill=GROUP)) + geom_histogram() + labs(title = "Histogram of Age Across All Groups", x = "AGE", y = "Number of Samples") 
ggp
#histogram of age distribution reveals the data is skewed and there is an outlier in the group 

ggp = ggplot(sample_info, aes(x=NEUTROPHILS, fill=GROUP)) + geom_histogram() + labs(title = "Histogram of Neutrophils Across All Groups", x = "Neutrophils", y = "Number of Samples") 
ggp

ggp = ggplot(sample_info, aes(x=MONOCYTES, fill=GROUP)) + geom_histogram() + labs(title = "Histogram of Monocytes Across All Groups", x = "Monocytes", y = "Number of Samples") 
ggp


ggp = ggplot(sample_info, aes(x=GROUP, y=AGE, fill=GROUP)) + geom_boxplot() + labs(title = "Boxplot of Age Across All Groups", x = "GROUP", y = "AGE") 
ggp
ggp
#box plot confirms outliers and indicates some of the data is skewed - however not a good assessment for normality of data distribution
#box plot also shows there is not spread/variance 

#Shapiro-Wilk test to evaluate the distribution of data (i.e. if it is normal)

shapiro.test(sample_info$AGE[sample_info$GROUP=='HEALTHY'])
#p value 0.332
shapiro.test(sample_info$MONOCYTES[sample_info$GROUP=='HEALTHY'])
#p value 0.886
shapiro.test(sample_info$NEUTROPHILS[sample_info$GROUP=='HEALTHY'])
#p value 0.862

shapiro.test(sample_info$AGE[sample_info$GROUP=='GOUT'])
#p value 0.369
shapiro.test(sample_info$MONOCYTES[sample_info$GROUP=='GOUT'])
#p value 0.124
shapiro.test(sample_info$NEUTROPHILS[sample_info$GROUP=='GOUT'])
#p value 0.281

shapiro.test(sample_info$AGE[sample_info$GROUP=='SA'])
#p value 0.0924
shapiro.test(sample_info$MONOCYTES[sample_info$GROUP=='SA'])
#p value 0.710
shapiro.test(sample_info$NEUTROPHILS[sample_info$GROUP=='SA'])
#p value 0.897

#while the plots indicate that there is heavily skewed data, the tests for normal distribution do not indicate the data is not normal 
   #it is important to acknowledge the limitations of these tests and to use visualization as a preferred method 
   #furthermore, the sample size of each group is quite small, which is known to heavily impact the results of these tests

#test for variance 
bartlett.test(MONOCYTES~GROUP, sample_info)
#p value 0.01117
bartlett.test(NEUTROPHILS~GROUP, sample_info)
#p value 0.0206
bartlett.test(AGE~GROUP, sample_info)
#p value 0.03163

#p value < 0.05 for all clinical information indicates the variance is not the same for all groups

#Kruskal-Wallis test  
kruskal.test(AGE ~ GROUP, data=sample_info)
#p value 0.7706 
kruskal.test(SEX ~ GROUP, data=sample_info)
#p value 0.34

#kruskal wallis test is used to compare several groups. P value indicates no significant difference 



############################################################################################################################################################



########TASK 2 - significantly differential genes between HCvGOUT and HCvSA#######


#merge tables to compile gene information 
expression_annotated = merge(Expression, data_annotate, by.x=0,by.y=0)
SAvsHC_annotated = merge(DE_SAvsHC, data_annotate, by.x=0, by.y=0)
GOUTvsHC_annotated = merge(DE_GOUTvsHC, data_annotate, by.x=0, by.y=0)

#rename rows - not using Gene Symbol due to duplicates 
row.names(expression_annotated) = expression_annotated[,1]
row.names(SAvsHC_annotated) = SAvsHC_annotated[,1]
row.names(GOUTvsHC_annotated) = GOUTvsHC_annotated[,1]

#rename columns 
names(expression_annotated)[1] = 'Gene_ID'
names(SAvsHC_annotated)[1] = 'Gene_ID'
names(GOUTvsHC_annotated)[1] = 'Gene_ID'

#get the significant genes 
HC_GOUT_sig = subset(GOUTvsHC_annotated, p.adj < 0.05)
HC_SA_sig = subset(SAvsHC_annotated, p.adj < 0.05)


#histogram of all the differential genes from HC_Gout and HC_SA 
ggp = ggplot(HC_GOUT_sig, aes(x=log2fold)) + geom_histogram()
ggp
ggp = ggplot(HC_SA_sig, aes(x=log2fold)) + geom_histogram()
ggp
ggp = ggplot(HC_GOUT_sig,aes(x=log2fold)) +geom_histogram(colour="blue") + geom_histogram(data=HC_SA_sig,alpha = 0.5)
ggp


#table for sig up and down regulated genes 
HC_GOUT_UP = subset(HC_GOUT_sig, p.adj < 0.05 & log2fold >1)
HC_GOUT_DOWN = subset(HC_GOUT_sig, p.adj < 0.05 & log2fold < -1)

HC_SA_UP = subset(HC_SA_sig, p.adj < 0.05 & log2fold >1)
HC_SAT_DOWN = subset(HC_SA_sig, p.adj < 0.05 & log2fold < -1)

#make into a volcano plot 
ggp = ggplot(HC_GOUT_sig, aes(x=log2fold, y=-log10(p.adj))) + geom_point(colour = "black") + geom_point(data = HC_GOUT_UP, colour = "red") + geom_point(data = HC_GOUT_DOWN, colour = "blue") + labs(title="HC_GOUT Significant Genes") 
ggp

ggp = ggplot(HC_SA_sig, aes(x=log2fold, y=-log10(p.adj))) + geom_point(colour = "black") + geom_point(data = HC_SA_UP, colour = "red") + geom_point(data = HC_SAT_DOWN, colour = "blue") + labs(title="HC_SA Significant Genes")
ggp

#volcano plot of both differential tables together 
ggp = ggplot(HC_GOUT_sig, aes(x=log2fold, y=-log10(p.adj))) + geom_point(colour = "black") + geom_point(data=HC_SA_sig, colour="red") + labs(title="HC_GOUT (Black) and HC_SA (Red)  Significant Genes") + geom_vline(xintercept=-1,linetype="dashed") + geom_vline(xintercept=1,linetype="dashed") + geom_hline(yintercept=-log10(0.05),linetype="dashed") 
ggp


#merge to get the number of common sig genes from the two comparative groups 
HC_GOUT_sig_gene_IDs = HC_GOUT_sig$Gene_ID
HC_SA_sig_gene_IDs = HC_SA_sig$Gene_ID
common_sig_genes = intersect(HC_GOUT_sig_gene_IDs, HC_SA_sig_gene_IDs)
#indicates there are 3661 common differential genes 

#looking to see how much the data was cut down by parsing for commonality 
3661/6956
3661/5552 


############################################################################################################################################################



########Task 3: Select Genes and Assess with Clinical Measurements##########



#merge to get the common sig genes from the two comparative groups 
sig_genes = merge(HC_GOUT_sig, HC_SA_sig, by=1)
row.names(sig_genes) = sig_genes[,1]

#expression table with only significant genes
sig_genes_ID_only = subset(sig_genes, select = -c(2:17))
sig_gene_expression = merge(sig_genes_ID_only, Expression, by=0)
row.names(sig_gene_expression) = sig_gene_expression[,1]
sig_gene_expression = subset(sig_gene_expression, select = -c(1:2))

sig_expression_data = merge(sig_genes_ID_only, Expression,by=0)
row.names(sig_expression_data) = sig_expression_data[,1]
sig_expression_data = subset(sig_expression_data, select = -c(1:2))


#trying to identify any potential patterns or clustering using PCA 

#PCA --> check if similarities in differential expression of these genes aligns with clinical information 
PCA = prcomp(t(sig_gene_expression))
pca_coordinates = data.frame(PCA$x)
summary(PCA)

ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_info$GROUP)) + geom_point() + labs(title="PCA Plot of Significant Differentially Expressed Genes in Healthy, Gout, and SA Groups")
ggp
ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_info$AGE)) + geom_point() + labs(title="PCA Plot of Significant Differentially Expressed Genes Across Age")
ggp
ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_info$SEX)) + geom_point()
ggp
ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_info$MONOCYTES)) + geom_point()
ggp
ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_info$NEUTROPHILS)) + geom_point()
ggp

PCA2 = prcomp(t(Expression))
pca_data = data.frame(PCA2$x)
summary(PCA2)


fviz_eig(PCA2, addlabels = TRUE)
67.9 + 13.5


ggp = ggplot(pca_data, aes(x=PC1, y=PC2, colour = sample_info$GROUP)) + geom_point() + labs(title="PCA Plot of Gene Expression in Healthy, Gout, and SA Groups")
ggp
ggp = ggplot(pca_data, aes(x=PC1, y=PC2, colour = sample_info$AGE)) + geom_point() + labs(title="PCA Plot of Significant Differentially Expressed Genes Across Age")
ggp
ggp = ggplot(pca_data, aes(x=PC1, y=PC2, colour = sample_info$SEX)) + geom_point()
ggp
ggp = ggplot(pca_data, aes(x=PC1, y=PC2, colour = sample_info$MONOCYTES)) + geom_point()
ggp
ggp = ggplot(pca_data, aes(x=PC1, y=PC2, colour = sample_info$NEUTROPHILS)) + geom_point()
ggp


#pick my significant genes based on p.adj --> one is the most significant for HC_GOUT and the other is most significant for HC_SA
#HC_GOUT gene: ENSG00000108950, FAM20A, at p.adj=3.322x10^-13
#HC_SA gene: ENSG00000170439, METTL7B, at p.adj=5.98x10^-18

#pull the expression levels of genes for all 3 groups 
FAM20A_expression = data.frame(t(sig_expression_data["ENSG00000108950",]))
METTL7B_expression = data.frame(t(sig_expression_data["ENSG00000170439",]))

#compile clinical data 
FAM20A_expression$GROUP = sample_info$GROUP
FAM20A_expression$AGE = sample_info$AGE
FAM20A_expression$SEX = sample_info$SEX
FAM20A_expression$NEUTROPHILS = sample_info$NEUTROPHILS
FAM20A_expression$MONOCYTES = sample_info$MONOCYTES

METTL7B_expression$GROUP = sample_info$GROUP
METTL7B_expression$AGE = sample_info$AGE
METTL7B_expression$SEX = sample_info$SEX
METTL7B_expression$NEUTROPHILS = sample_info$NEUTROPHILS
METTL7B_expression$MONOCYTES = sample_info$MONOCYTES

FAM20A_healthy_data = subset(FAM20A_expression, GROUP == 'HEALTHY')
FAM20A_gout_data = subset(FAM20A_expression, GROUP == 'GOUT')
FAM20A_sa_data = subset(FAM20A_expression, GROUP == 'SA')

METTL7B_healthy_data = subset(METTL7B_expression, GROUP == 'HEALTHY')
METTL7B_gout_data = subset(METTL7B_expression, GROUP == 'GOUT')
METTL7B_sa_data = subset(METTL7B_expression, GROUP == 'SA')

#plot data for visualization 
ggp = ggplot(FAM20A_expression, aes(x=ENSG00000108950, fill=GROUP)) + geom_density()
ggp
ggp = ggplot(METTL7B_expression, aes(x=ENSG00000170439, fill=GROUP)) + geom_density()
ggp
#both distributions of expression are skewed to the right, with METTL plot showing a wider spread of expression levels for SA samples 
#sample size is small and distribution does not appear to be normal in either gene -- opting for non parametric test 

#statistical testing for the two genes: FAM20A vs Gout
cor(FAM20A_gout_data$ENSG00000108950, FAM20A_gout_data$AGE, method = c("spearman"))
#-0.07717793
cor(FAM20A_gout_data$ENSG00000108950, FAM20A_gout_data$NEUTROPHILS, method = c("spearman"))
#0.5384615
cor(FAM20A_gout_data$ENSG00000108950, FAM20A_gout_data$MONOCYTES, method = c("spearman"))
#-0.1032967

#statistical testing for the two genes: FAM20A vs SA 
cor(FAM20A_sa_data$ENSG00000108950, FAM20A_sa_data$AGE, method = c("spearman"))
#0.3142857
cor(FAM20A_sa_data$ENSG00000108950, FAM20A_sa_data$NEUTROPHILS, method = c("spearman"))
#0.4417582
cor(FAM20A_sa_data$ENSG00000108950, FAM20A_sa_data$MONOCYTES, method = c("spearman"))
#0.3098901

#statistical testing for the two genes: FAM20A vs Healthy 
cor(FAM20A_healthy_data$ENSG00000108950, FAM20A_healthy_data$AGE, method = c("spearman"))
#-0.3588148
cor(FAM20A_healthy_data$ENSG00000108950, FAM20A_healthy_data$NEUTROPHILS, method = c("spearman"))
#0.578022
cor(FAM20A_healthy_data$ENSG00000108950, FAM20A_healthy_data$MONOCYTES, method = c("spearman"))
#-0.432967

#statistical testing for the two genes: METTL7B vs Gout 
cor(METTL7B_gout_data$ENSG00000170439, METTL7B_gout_data$AGE, method = c("spearman"))
#-0.105844
cor(METTL7B_gout_data$ENSG00000170439, METTL7B_gout_data$NEUTROPHILS, method = c("spearman"))
#0.7846154
cor(METTL7B_gout_data$ENSG00000170439, METTL7B_gout_data$MONOCYTES, method = c("spearman"))
#-0.06813187

#statistical testing for the two genes: METTL7B vs. SA 
cor(METTL7B_sa_data$ENSG00000170439, METTL7B_sa_data$AGE, method = c("spearman"))
#0.3846154
cor(METTL7B_sa_data$ENSG00000170439, METTL7B_sa_data$NEUTROPHILS, method = c("spearman"))
#0.6967033
cor(METTL7B_sa_data$ENSG00000170439, METTL7B_sa_data$MONOCYTES, method = c("spearman"))
#-0.07692308

#statistical testing for the two genes: METTL7B vs Healthy 
cor(METTL7B_healthy_data$ENSG00000170439, METTL7B_healthy_data$AGE, method = c("spearman"))
#-0.1130845
cor(METTL7B_healthy_data$ENSG00000170439, METTL7B_healthy_data$NEUTROPHILS, method = c("spearman"))
#0.2068208
cor(METTL7B_healthy_data$ENSG00000170439, METTL7B_healthy_data$MONOCYTES, method = c("spearman"))
#-0.2024204

#plots 
ggp = ggplot(FAM20A_expression, aes(x=AGE, y=ENSG00000108950, colour=GROUP)) + geom_point() + labs(title="Relationship between Age and FAM20A Gene Expression")
ggp
ggp = ggplot(METTL7B_expression, aes(x=AGE, y=ENSG00000170439, colour=GROUP)) + geom_point() + labs(title="Relationship between Age and METTYL7B Gene Expression")
ggp
ggp = ggplot(FAM20A_expression, aes(x=NEUTROPHILS, y=ENSG00000108950, colour=GROUP)) + geom_point() + labs(title="Relationship between Neutrophils and FAM20A Gene Expression")
ggp
ggp = ggplot(METTL7B_expression, aes(x=NEUTROPHILS, y=ENSG00000170439, colour=GROUP)) + geom_point() + labs(title="Relationship between Neutrophils and METTYL7B Gene Expression")
ggp
ggp = ggplot(FAM20A_expression, aes(x=MONOCYTES, y=ENSG00000108950, colour=GROUP)) + geom_point() + labs(title="Relationship between Monocytes and FAM20A Gene Expression")
ggp
ggp = ggplot(METTL7B_expression, aes(x=MONOCYTES, y=ENSG00000170439, colour=GROUP)) + geom_point() + labs(title="Relationship between Monocytes and METTYL7B Gene Expression")
ggp
#plot point confirms that no clear relationship can be derived from these genes 

 

############################################################################################################################################################



###########Task 4: Significantly Different Genes Between Gout and SA##############


#separate tables based on group: HC, GOUT, SA 
HC_expression = subset(Expression, select = c(1:14))
GOUT_expression = subset(Expression, select = c(15:28))
SA_expression = subset(Expression, select = c(29:42))

#plot samples from both groups 
ggp = ggplot(GOUT_expression, aes(x=Gout_1)) + geom_histogram()
ggp

ggp = ggplot(SA_expression, aes(x=SA_1)) + geom_histogram()
ggp


#table to store the results 
de_GOUTvsSA = as.data.frame(matrix(0,ncol = 3, nrow = nrow(Expression)))
names(de_GOUTvsSA) = c("Log2Fold","p","p.adj")
row.names(de_GOUTvsSA) = row.names(Expression)

#create loop to get values for every gene 
for(row in 1:nrow(Expression))
{
  gene_data_GOUT = as.numeric(GOUT_expression[row,])
  gene_data_SA = as.numeric(SA_expression[row,])
  mean_GOUT = mean(gene_data_GOUT)
  mean_SA = mean(gene_data_SA)
  log2fold = log2(mean_SA) - log2(mean_GOUT)
  p = t.test(gene_data_SA,gene_data_GOUT)
  p = p$p.value
  
  de_GOUTvsSA[row,"Log2Fold"] = log2fold
  de_GOUTvsSA[row,"p"] = p
}

#add p.adj values 
pvalues = de_GOUTvsSA$p
p.adjusted = p.adjust(pvalues, method="bonferroni")
de_GOUTvsSA$p.adj = p.adjusted

p.adjusted = p.adjust(pvalues, method="BH")
de_GOUTvsSA$p.adjBH = p.adjusted


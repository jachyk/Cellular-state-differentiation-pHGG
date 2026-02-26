# Bulk RNA-seq analysis using OpenPBTA data for
# Signature scoring, regression analysis, subtype comparison 
# and statistical testing.

# Data is derived from PedCBioportal;
# for further information see: https://github.com/d3b-center/OpenPedCan-analysis


#----- Load packages -----#

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(openxlsx)
library(PMCMRplus)


options(scipen = 15)
options(stringsAsFactors = FALSE)


#----- Define input/output directories -----#

expression_file <- "XXXX/.txt"
clinical_file   <- "XXXX/.txt"

savedir <- "/Output/"


# subtype color palette
color_subtype1 <- c(
  "Diffuse midline glioma" = "#264653",
  "Diffuse hemispheric glioma" = "#e9c46a",
  "High-grade glioma" = "#d62828",
  "Infant-type hemispheric glioma" = "#2a9d8f"
)


#----- Load and format bulk expression and clinical metadata -----#

expr_data     <- read.delim(expression_file, header = TRUE, row.names = 1)
clinical_data <- read.delim(clinical_file, header = TRUE)

head(expr_data)
head(clinical_data)


# format clinical metadata
clinical_data2 <- clinical_data[5:nrow(clinical_data), ]
colnames(clinical_data2) <- clinical_data[4,]
unique(colnames(clinical_data2))

# format expression matrix sample IDs
expr_data2 <- expr_data
colN <- colnames(expr_data2)

for (i in seq_along(colN)) {
  if (startsWith(colN[i], "X")) {
    colN[i] <- substring(colN[i], 2)
  }
  colN[i] <- gsub("\\.", "-", colN[i])
}

colnames(expr_data2) <- colN


#----- Filter dataset for tumor entities of interest -----#

Toinclude <- c(
  "High-grade glioma, IDH-wildtype and H3-wildtype",
  "Diffuse midline glioma, H3 K28-altered",
  "Infant-type hemispheric glioma, ALK-altered",
  "Infant-type hemispheric glioma, NTRK-altered",
  "Diffuse hemispheric glioma, H3 G35-mutant",
  "Diffuse midline glioma, H3 K28-altered (EGFR mutant)",
  "Infant-type hemispheric glioma, ROS1-altered"
)

clinical_d_filtered <- clinical_data2[ clinical_data2$CANCER_TYPE_DETAILED %in% Toinclude,]

# filter expression matrix
metadata_sample_ids <- clinical_d_filtered$MATCHED_NORMAL_SAMPLE_ID
filtered_expr_data2 <-
  expr_data2[, colnames(expr_data2) %in% metadata_sample_ids]
counts <- filtered_expr_data2
Toinclude2 <- colnames(filtered_expr_data2)
clinical_d_filtered2 <-
  clinical_d_filtered[
    clinical_d_filtered$MATCHED_NORMAL_SAMPLE_ID %in% Toinclude2,]

# since some sample IDs are duplicated, only keep representative entries
clinical_d_filtered3 <-
  clinical_d_filtered2 %>%
  distinct(MATCHED_NORMAL_SAMPLE_ID, .keep_all = TRUE)

colnames(clinical_d_filtered3)[21] <- "Sample_ID"

#----- Cohort composition -----#

sum(clinical_d_filtered3$CANCER_GROUP == "Infant-type hemispheric glioma")
sum(clinical_d_filtered3$CANCER_GROUP == "Diffuse midline glioma")
sum(clinical_d_filtered3$CANCER_GROUP == "High-grade glioma")
sum(clinical_d_filtered3$CANCER_GROUP == "Diffuse hemispheric glioma")


#----- normalize expression -----# 
counts <- filtered_expr_data2
counts_log<- log2(counts + 1)



#----- Regression analysis between signatures -----#
# as shown in Figure 4
name1 <- "OPC"
signature1 <- c("OLIG1","OLIG2","SOX10","GPR17","TNR","CSPG4","MBP","PLP1","MBP")
signature_expression1 <- colMeans(counts_log[signature1, ])

name2 <- "INP"
signature2 <- c("DCX","CD24","SOX4","DLX5","DLX6","DLX2","DLX1","NRXN3")
signature_expression2 <- colMeans(counts_log[signature2, ])

data <- data.frame(
  Signature1 = signature_expression1,
  Signature2 = signature_expression2
)

ggplot(data, aes(x = Signature1, y = Signature2)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = paste0("Relation of ", name1, " vs ", name2),
    x = paste0(name1, " Expression"),
    y = paste0(name2, " Expression")
  )

#statistical testing 
model <- lm(Signature1 ~ Signature2, data = data)
summary(model)

cor(data$Signature1, data$Signature2, method = "pearson")


#----- Regression analysis between individual genes -----#

Gene1 <- "NTRK2"
Gene2 <- "SPRY4"

single_expression1 <- t(counts_log[Gene1, ])
single_expression2 <- t(counts_log[Gene2, ])

data <- data.frame(Gene1 = single_expression1,
                   Gene2 = single_expression2)

data$Sample_ID <- rownames(data)

mdata <- merge(clinical_d_filtered3, data, by = "Sample_ID")

ggplot(data = mdata, aes(x = Gene1, y = Gene2, col = ))+ 
  geom_point()+
  geom_smooth(method = "lm", color = "red") +
  labs(title= paste0("Relation of ",Gene1, " Expression vs. ",Gene2," Expression"), x = paste0(Gene1," Expression"), y = paste0(Gene2," Expression"))+
  theme(axis.text.x = element_text(size = 14, face= "bold")) +
  theme(axis.text.y = element_text(size = 14, face= "bold")) +
  theme(title = element_text(size = 14, face= "bold"))+
  theme(legend.position = "none")


# Statistical testing 
# print statistics to that plot and check wether results are significant
model<- lm(formula = Gene1 ~ Gene2, data = data)
summary(model)
Coeff <- cor(data$Gene1, data$Gene2, method = "pearson")
Coeff


#----- Marker expression comparison -----#

marker <- c("MOXD1")
markers <- c("AQP4","SLC1A2","SPARCL1")

counts_log <- log2(filtered_expr_data2 + 1)

#----- calculate values -----#
# When you want to plot a group of markers
signature_expression <- colMeans(counts_log[markers,])
head(signature_expression)
data <- data.frame(Signaturescore = signature_expression)
data$Sample_ID <- names(signature_expression)

# When you only want to plot ONE marker 
signature_expression <- counts_log[marker,]
head(signature_expression)
data <- data.frame(Signaturescore = signature_expression)
data <- as.data.frame(t(data))
data$Sample_ID <- names(signature_expression)

#------ prepare plotting----#
colnames(data)[1] <- "Signaturescore"
mdata <- merge(clinical_d_filtered3, data, by = "Sample_ID")

#set order of the x-axis
X<- "Diffuse hemispheric glioma"
XX<- "Infant-type hemispheric glioma"
XXX<- "Diffuse midline glioma"
XXXX<- "High-grade glioma"
mdata$CANCER_GROUP <- factor(mdata$CANCER_GROUP, levels = c(XXX,XXXX,X,XX))

Title <- paste0(marker," expression")

paste0(Title)


tiff(filename = paste0( exp_dir,Title,".tiff"),width = 450, height = 500)
print(ggplot(mdata,aes(x = CANCER_GROUP, y = Signaturescore, fill = CANCER_GROUP)) +
        geom_boxplot()+
        ggtitle(Title)+
        scale_fill_manual(values = color_subtype1, name ="")+
        theme_classic()+
        ylab("Expression")+
        theme(axis.title.x = element_blank(),
              title = element_text(size = 20),
              axis.title.y = element_text(size = 18,vjust = 2),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 19 ),
              legend.title = element_text(size = 16),
              legend.position = "none",
              legend.text = element_text(size = 18)))
dev.off()

#----- Statistical testing -----#

results <- aov(Signaturescore ~ CANCER_GROUP, data = mdata)
summary(results)

#post-hoc testing 
TK <-TukeyHSD(results)
games_howell<- gamesHowellTest(Signaturescore ~ CANCER_GROUP, data = mdata)
welch_results <- oneway.test(
          Signaturescore ~ CANCER_GROUP,
          data = mdata,
          var.equal = FALSE)


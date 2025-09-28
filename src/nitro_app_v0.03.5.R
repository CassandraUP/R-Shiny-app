library(shiny)

if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(shinyjs)){
  install.packages("shinyjs")
  library(shinyjs)
}
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}
if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}
if(!require(shinythemes)){
  install.packages("shinythemes")
  library(shinythemes)
}
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
if(!require(lubridate)){
  install.packages("lubridate")
  library(lubridate)
}
if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require(ComplexHeatmap)){
  BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap)
}
if(!require(circlize)){
  install.packages("circlize")
  library(circlize)
}
if(!require(dendextend)){
  install.packages("dendextend")
  library(dendextend)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
if(!require(grid)){
  install.packages("grid")
  library(grid)
}
if(!require(DT)){
  install.packages("DT")
  library(DT)
}
if(!require(rclipboard)){
  install.packages("rclipboard")
  library(rclipboard)
}
if(!require(plotrix)){
  install.packages("plotrix")
  library(plotrix)
}
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(ggtext)){
  install.packages("ggtext")
  library(ggtext)
}

options(warn = -1)

### Read in data

# Orthogroup annotation - can contain duplicates
ortho_annot <- read.csv("./nitro_ortho_annot.csv", sep = ",")
ortho_annot <- ortho_annot[!is.na(ortho_annot$Orthogroup), ]

# Data to plot bar charts for all species by orthogroup with standard error bars
ortho_bar_data <- read.csv("./nitro_bar_data_ortho.csv", sep = ",")
ortho_bar_data <- ortho_bar_data[!is.na(ortho_bar_data$Orthogroup), ]
ortho_bar_data$cols <- as.character(ortho_bar_data$cols)

# All E. natalensis data for tissue specific bar charts and heatmaps
ena_data_all <- read.csv("./nitro_app_ena_data.csv")

### Format data for bar charts
ena_tissue_bar <- ena_data_all[, c(1,2,25:43)]

ena_tissue_bar$sig_TF <- is.na(ena_tissue_bar$ena_signif)
ena_tissue_bar$sig_TF <- gsub(TRUE, "", ena_tissue_bar$sig_TF)
ena_tissue_bar$sig_TF <- gsub(FALSE, "*", ena_tissue_bar$sig_TF)
ena_tissue_bar$id_lab <- paste(ena_tissue_bar$sig_TF, ena_tissue_bar$ena_id, sep = "")

ena_tissue_bar_mean <- colnames(ena_tissue_bar[4:12])
ena_tissue_bar_se <- colnames(ena_tissue_bar[13:21])
ena_tissue_bar_melt <- melt(setDT(ena_tissue_bar), measure = list(ena_tissue_bar_mean, ena_tissue_bar_se), variable.name = "Tissue", value.name = c("mean", "se"))
ena_tissue_bar_melt$Tissue <- gsub(1, "CR", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(2, "GS", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(3, "IL", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(4, "IR", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(5, "LPR", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(6, "ML", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(7, "MR", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(8, "ST", ena_tissue_bar_melt$Tissue)
ena_tissue_bar_melt$Tissue <- gsub(9, "UPR", ena_tissue_bar_melt$Tissue)

### Assign colors to tissues
ena_tissue_bar_melt$cols <- ena_tissue_bar_melt$Tissue
ena_tissue_bar_melt$cols <- gsub("CR", "#00AFAF", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("GS", "#851E76", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("IL", "#B9E9AF", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("IR", "#80D86F", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("LPR", "#94564C", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("ML", "#46B030", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("MR", "#399027", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("ST", "#C08E87", ena_tissue_bar_melt$cols)
ena_tissue_bar_melt$cols <- gsub("UPR", "#6A3E37", ena_tissue_bar_melt$cols)

### Make color palette
ena_bar_cols <- unique(ena_tissue_bar_melt[ , c(6,9)])
ena_bar_my_cols <- ena_bar_cols$cols
names(ena_bar_my_cols) <- ena_bar_cols$Tissue
ena_bar_cols_tissue <- scale_fill_manual(name = "Tissue", values = ena_bar_my_cols) 

### Format E.nat tissue data for heatmap
ena_hm_mat <- as.matrix(ena_data_all[ , 26:34])
row.names(ena_hm_mat) <- ena_data_all$ena_id
colnames(ena_hm_mat) <- c("CR", "GS", "IL", "IR", "LPR", "ML", "MR", "ST", "UPR")
ena_hm_mat <- t(scale(t(ena_hm_mat), center = T, scale = T))
ena_hm_mat[is.nan(ena_hm_mat)] = 0

tissue_annot_test <- c("CR", "GS", "IL", "IR", "LPR", "ML", "MR", "ST", "UPR")
tissue_annot <- factor(tissue_annot_test, levels=c("CR","UPR","LPR","ST","IR","MR","IL","ML","GS"))
ena_tissue_col_annot <- HeatmapAnnotation(Tissue = tissue_annot, 
                                          col = list(Tissue= c("CR" = "#00AFAF", "GS" = "#851E76", "IL" = "#B9E9AF", "IR" = "#80D86F", 
                                                               "LPR" = "#94564C", "ML" = "#46B030", "MR" = "#399027", "ST" = "#C08E87", "UPR" = "#6A3E37")))
tissue_data <- ena_data_all[ , c(1,2,25,26:34)]
tissue_data$max <- apply(ena_data_all[ , 4:12], 1, max)
tissue_data$max[tissue_data$max < 1] <- NA
tissue_data <- mutate(tissue_data, percentile_rank = ntile(tissue_data$max,100))
tissue_data$max[is.na(tissue_data$max)] <- 0
tissue_data$percentile_rank[is.na(tissue_data$percentile_rank)] <- 0

### Heatmap percentile colors
per_row_col = circlize::colorRamp2(c(0, 100), c("white", "darkgray"))

### Read in data and format for other species heatmaps

#Azolla

afi_data <- read.csv("./afi_data.csv", header = T, sep = ",")
afi_hm_mat <- afi_data[ , 3:11]
afi_hm_mat <- as.matrix(afi_hm_mat)
row.names(afi_hm_mat) <- afi_data$afi_id
colnames(afi_hm_mat) <- c("NP_CM_1", "NP_CM_2", "NP_CM_3", "NP_CP_1", "NP_CP_2", "NP_CP_3", "NM_CP_1", "NM_CP_2", "NM_CP_3")
afi_hm_mat <- t(scale(t(afi_hm_mat), center = T, scale = T))
afi_hm_mat[is.nan(afi_hm_mat)] = 0

afi_annot <- c(rep("Nitro+ Cyano-",3), rep("Nitro+ Cyano+",3),rep("Nitro- Cyano+",3))
afi_annot <- factor(afi_annot, levels=c("Nitro- Cyano+", "Nitro+ Cyano+","Nitro+ Cyano-"))
afi_annot_col <- HeatmapAnnotation(Experiment = afi_annot, col = list(Experiment = c("Nitro+ Cyano-" = "#F8F8F8", "Nitro+ Cyano+" = "#FFC2C4","Nitro- Cyano+" = "#FF595E")))

afi_sample_data <- afi_data[ , 2:11]
colnames(afi_sample_data) <- c("afi_id","NP_CM_1", "NP_CM_2", "NP_CM_3", "NP_CP_1", "NP_CP_2", "NP_CP_3", "NM_CP_1", "NM_CP_2", "NM_CP_3")
afi_sample_data$max <- apply(afi_sample_data[ , 2:10], 1, max)
afi_sample_data$max[afi_sample_data$max < 1] <- NA
afi_sample_data <- mutate(afi_sample_data, percentile_rank = ntile(afi_sample_data$max,100))
afi_sample_data$max[is.na(afi_sample_data$max)] <- 0
afi_sample_data$percentile_rank[is.na(afi_sample_data$percentile_rank)] <- 0

afi_annot <- read.csv("afi_annot.csv", header = T, sep = ",")

#Anthoceros

apu_data <- read.csv("./apu_data.csv", header = T, sep = ",")
apu_hm_mat <- apu_data[ , 3:8]
apu_hm_mat <- as.matrix(apu_hm_mat)
row.names(apu_hm_mat) <- apu_data$apu_id
colnames(apu_hm_mat) <- c("NM_CP_1", "NM_CP_2", "NM_CP_3", "NP_CM_1", "NP_CM_2", "NP_CM_3")
apu_hm_mat <- t(scale(t(apu_hm_mat), center = T, scale = T))
apu_hm_mat[is.nan(apu_hm_mat)] = 0

apu_annot <- c(rep("Nitro- Cyano+",3), rep("Nitro+ Cyano-",3))
apu_annot <- factor(apu_annot, levels=c("Nitro- Cyano+", "Nitro+ Cyano-"))
apu_annot_col <- HeatmapAnnotation(Experiment = apu_annot, col = list(Experiment = c("Nitro+ Cyano-" = "#F8F8F8", "Nitro- Cyano+" = "#FFCA3A")))

apu_sample_data <- apu_data[ , 2:8]
colnames(apu_sample_data) <- c("apu_id", "NM_CP_1", "NM_CP_2", "NM_CP_3","NP_CM_1", "NP_CM_2", "NP_CM_3")
apu_sample_data$max <- apply(apu_sample_data[ , 2:7], 1, max)
apu_sample_data$max[apu_sample_data$max < 1] <- NA
apu_sample_data <- mutate(apu_sample_data, percentile_rank = ntile(apu_sample_data$max,100))
apu_sample_data$max[is.na(apu_sample_data$max)] <- 0
apu_sample_data$percentile_rank[is.na(apu_sample_data$percentile_rank)] <- 0

apu_annot <- read.csv("apu_annot.csv", header = T, sep = ",")

#Dactylis glomerata

dgl_data <- read.csv("./dgl_data.csv", header = T, sep = ",")
dgl_hm_mat <- dgl_data[ , 3:14]
dgl_hm_mat <- as.matrix(dgl_hm_mat)
row.names(dgl_hm_mat) <- dgl_data$dgl_id
colnames(dgl_hm_mat) <- c("nodule_1", "nodule_2", "nodule_3", "nodule_4", "nodule_5", "nodule_6", "root_1", "root_2", "root_3", "root_4", "root_5", "root_6")
dgl_hm_mat <- t(scale(t(dgl_hm_mat), center = T, scale = T))
dgl_hm_mat[is.nan(dgl_hm_mat)] = 0

dgl_annot <- c(rep("Nodule",6), rep("Root",6))
dgl_annot <- factor(dgl_annot, levels=c("Nodule", "Root"))
dgl_annot_col <- HeatmapAnnotation(Tissue = dgl_annot, col = list(Tissue = c("Nodule" = "#6A4C93", "Root" = "#F8F8F8")))

dgl_sample_data <- dgl_data[ , 2:14]
colnames(dgl_sample_data) <- c("dgl_id", "nodule_1", "nodule_2", "nodule_3", "nodule_4", "nodule_5", "nodule_6", "root_1", "root_2", "root_3", "root_4", "root_5", "root_6")
dgl_sample_data$max <- apply(dgl_sample_data[ , 2:13], 1, max)
dgl_sample_data$max[dgl_sample_data$max < 1] <- NA
dgl_sample_data <- mutate(dgl_sample_data, percentile_rank = ntile(dgl_sample_data$max,100))
dgl_sample_data$max[is.na(dgl_sample_data$max)] <- 0
dgl_sample_data$percentile_rank[is.na(dgl_sample_data$percentile_rank)] <- 0

dgl_annot <- read.csv("dgl_annot.csv", header = T, sep = ",")

### Gunnera perpensa

gpe_data <- read.csv("./gpe_data.csv", header = T, sep = ",")
gpe_hm_mat <- gpe_data[ , 4:9]
gpe_hm_mat <- as.matrix(gpe_hm_mat)
row.names(gpe_hm_mat) <- gpe_data$gpe_id
colnames(gpe_hm_mat) <- c("Control 1", "Control 2", "Control 3", "Cyano+ 1", "Cyano+ 2", "Cyano+ 3")
gpe_hm_mat <- t(scale(t(gpe_hm_mat), center = T, scale = T))
gpe_hm_mat[is.nan(gpe_hm_mat)] = 0

gpe_annot <- c(rep("Control",3), rep("Cyano+",3))
gpe_annot <- factor(gpe_annot, levels=c("Cyano+", "Control"))
gpe_annot_col <- HeatmapAnnotation(Tissue = gpe_annot, col = list(Tissue = c("Cyano+" = "#1982C4", "Control" = "#F8F8F8")))

gpe_sample_data <- gpe_data[ , c(2,4:9)]
colnames(gpe_sample_data) <- c("gpe_id", "Control 1", "Control 2", "Control 3", "Cyano+ 1", "Cyano+ 2", "Cyano+ 3")
gpe_sample_data$max <- apply(gpe_sample_data[ , 2:7], 1, max)
gpe_sample_data$max[gpe_sample_data$max < 1] <- NA
gpe_sample_data <- mutate(gpe_sample_data, percentile_rank = ntile(gpe_sample_data$max,100))
gpe_sample_data$max[is.na(gpe_sample_data$max)] <- 0
gpe_sample_data$percentile_rank[is.na(gpe_sample_data$percentile_rank)] <- 0

gpe_annot <- read.csv("./gpe_annot.csv", header = T, sep = ",")

### Medicago truncatula

mtr_data <- read.csv("./mtr_data.csv", header = T, sep = ",")
mtr_hm_mat <- mtr_data[ , 4:9]
mtr_hm_mat <- as.matrix(mtr_hm_mat)
row.names(mtr_hm_mat) <- mtr_data$mtr_id
colnames(mtr_hm_mat) <- c("Root 1", "Root 2", "Root 3", "Nodule 1", "Nodule 2", "Nodule 3")
mtr_hm_mat <- t(scale(t(mtr_hm_mat), center = T, scale = T))
mtr_hm_mat[is.nan(mtr_hm_mat)] = 0

mtr_annot <- c(rep("Root",3), rep("Nodule",3))
mtr_annot <- factor(mtr_annot, levels=c("Nodule", "Root"))
mtr_annot_col <- HeatmapAnnotation(Tissue = mtr_annot, col = list(Tissue = c("Nodule" = "#6A4C93", "Root" = "#F8F8F8")))

mtr_sample_data <- mtr_data[ , c(2,4:9)]
colnames(mtr_sample_data) <- c("mtr_id", "Root 1", "Root 2", "Root 3", "Nodule 1", "Nodule 2", "Nodule 3")
mtr_sample_data$max <- apply(mtr_sample_data[ , 2:7], 1, max)
mtr_sample_data$max[mtr_sample_data$max < 1] <- NA
mtr_sample_data <- mutate(mtr_sample_data, percentile_rank = ntile(mtr_sample_data$max,100))
mtr_sample_data$max[is.na(mtr_sample_data$max)] <- 0
mtr_sample_data$percentile_rank[is.na(mtr_sample_data$percentile_rank)] <- 0

mtr_annot <- read.csv("./mtr_annot.csv", header = T, sep = ",")


### Ena annotation for convert
ena_annot <- read.csv("./ena_annot_best.csv", header = T, sep = ",")
colnames(ena_annot) <- c("ena_id", "Gu_per", "Az_fil", "An_pun", "Me_tru", "ath_id", "ath_description", 
                         "Orthogroup", "Singleton", "KO", "KEGG_definintion", "mm_bincode", "mm_definition",
                         "pfam_id", "pfam_description", "interpro_id", "interpro_description", "GO_slim")
ena_annot <- ena_annot[ ,c(1,6:18)]
ena_best_hits <- read.csv("ena_annot_best_hits_all.csv", header = T, sep = ",")
ena_best_any <- ena_best_hits %>% gather(Species, id, -En_nat, -Orthogroup)

### page_setup_ortho_bar
sidebar_content_bar_ortho <- sidebarPanel(
  textInput("text_bar_ortho", label = h5("Enter orthogroup ID"), value = "OG0000380"),
  p("Plot download options"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("bar_ortho_height", "Height:", 7, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("bar_ortho_width", "Width:", 10, min = 1, max = 100)), 
  submitButton("Update plot")
)
main_content_bar_ortho <- mainPanel(
  p("Plot a bar chart of all genes within an orthogroup by infected (filled) vs uninfected (unfilled) samples"),
  h5("Please cite ***Danielle & Cas to send appropriate references"),
  plotOutput("barPlot_ortho", height = 500),
  h5("Bold indicates significant differential expression between infected (I - ) and Uninfected (U - ) tissues"),
  downloadButton(outputId = "down_bar_ortho", label = "Download plot"),
  div(tableOutput('tbl_bar_ortho'), style = "font-size:80%")
)

bar_panel_ortho <- tabPanel(
  "Orthogroup specific bar chart - All species",
  sidebarLayout(
    sidebar_content_bar_ortho, main_content_bar_ortho
  )
)

### page_setup_ena_ortho_bar
sidebar_content_bar_ortho_ena <- sidebarPanel(
  textInput("text_bar_ortho_ena", label = h5("Enter orthogroup ID"), value = "OG0000380"),
  submitButton("Update plot")
)
main_content_bar_ortho_ena <- mainPanel(
  p("Plot a bar chart of all genes within an orthogroup across E. natalensis tissues"),
  h5("Please cite ***Danielle & Cas in prep"),
  plotOutput("barPlot_ortho_ena", height = 500),
  downloadButton(outputId = "down_bar_ortho_ena", label = "Download plot"),
  div(tableOutput('tbl_bar_ortho_ena'), style = "font-size:80%")
)

bar_panel_ortho_ena <- tabPanel(
  "Orthogroup specific bar chart - E. natalensis",
  sidebarLayout(
    sidebar_content_bar_ortho_ena, main_content_bar_ortho_ena
  )
)

### Page setup E. natalensis tissue heatmap

sidebar_content_hm_tissue <- sidebarPanel(
  shiny::textAreaInput("genes_hm_tissue", 'Input gene names separated by new line:', 
                       value = paste("Enat_0007569", "Enat_0004862", "Enat_0011750", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_ena_tissue", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_tissue_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("tissue_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_tissue", label = h5("Plot name"), value = "Enat_tissue"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_tissue_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_tissue_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_tissue <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in E. natalensis"),
  h5("Please cite ***"),
  plotOutput("hmPlot_tissue"),
  downloadButton(outputId = "hm_down_tissue", label = "Download plot"),
  downloadButton("download_tissue_clust", "Download result table"),
  div(tableOutput('hm_tbl_tissue'), style = "font-size:80%")
)


hm_panel_tissue <- tabPanel(
  "E. natalensis tissue-specific Heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_tissue, main_content_hm_tissue
  )
)

### Page setup - convert

sidebar_content_convert <- sidebarPanel(
    shiny::selectInput("convert_from", label = h5("From    "), 
                       choices = list("Any" = "Any", 
                                      "En. natalenesis" = "En_nat", 
                                      "Me. truncatula" = "Me_tru",
                                      "At. thaliana" = "At_tha"), 
                       selected = "Any"),
    shiny::selectInput("convert_to", label = h5("To     "), 
                       choices = list("En. natalenesis" = "En_nat", 
                                      "Me. truncatula" = "Me_tru",
                                      "At. thaliana" = "At_tha"), 
                       selected = "En_nat"),
    shiny::textAreaInput("convert_genes", 'Input gene names separated by new line:', 
                         value = paste("G_perpensa_evgTRINITY_DN16775_c5_g1", "Apun_evm.model.utg000031l.60.1", "Azfi_s0001.g000194", "Medtr1g026110", "AT3G17470", sep = "\n"), 
                         width = NULL),
  div(style="display:inline-block",submitButton("Update output")),
  rclipboardSetup(),
  div(style="display:inline-block",uiOutput("clip")),
  textInput("name_converter", label = h5("Name for download"), value = "list"),
  downloadButton("download_converter", "Download conversion table")
)  

main_content_convert <- mainPanel(
  p("Convert gene IDs between species"),
  tableOutput('convert_tbl')
)

convert_panel <- tabPanel(
  "Convert gene IDs",
  sidebarLayout(
    sidebar_content_convert, main_content_convert
  )
)

### Page setup Ena tissue bar charts & best hits

sidebar_content_ena_tissue <- sidebarPanel(
  textInput("text_bar_tissue_ena", label = h5("Enter E. natalensis ID"), value = "Enat_0016279"),
  submitButton("Update plot")
)

main_content_ena_tissue <- mainPanel(
  p("Plot a bar chart of a single E. natalensis gene across tissues"),
  h5("Please cite ***Danielle & Cas in prep"),
  plotOutput("barPlot_tissue", height = 500),
  downloadButton(outputId = "down_bar_tissue", label = "Download plot"),
  div(tableOutput('tbl_bar_tissue'), style = "font-size:80%")
)

ena_panel_tissue <- tabPanel(
  "Tissue bar chart - E. natalensis",
  sidebarLayout(
    sidebar_content_ena_tissue, main_content_ena_tissue
  )
)

### Page setup A. filiculoides sample heatmap

sidebar_content_hm_afi <- sidebarPanel(
  shiny::textAreaInput("genes_hm_afi", 'Input gene names separated by new line:', 
                       value = paste("Azfi_s0641.g080692", "Azfi_s0020.g015478", "Azfi_s0008.g011616", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_afi", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_afi_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("afi_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_afi", label = h5("Plot name"), value = "Azfi_hm"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_afi_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_afi_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_afi <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in A. filiculoides"),
  h5("Please cite ***"),
  shiny::plotOutput("hmPlot_afi"),
  downloadButton(outputId = "hm_down_afi", label = "Download plot"),
  downloadButton("download_afi_clust", "Download result table"),
  div(shiny::tableOutput('hm_tbl_afi'), style = "font-size:80%")
)


hm_panel_afi <- tabPanel(
  "A. filiculoides heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_afi, main_content_hm_afi
  )
)


### Page setup A. punctatus sample heatmap

sidebar_content_hm_apu <- sidebarPanel(
  shiny::textAreaInput("genes_hm_apu", 'Input gene names separated by new line:', 
                       value = paste("Apun_evm.model.utg000023l.1510.1", "Apun_evm.model.utg000135l.173.1", "Apun_evm.model.utg000173l.125.1", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_apu", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_apu_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("apu_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_apu", label = h5("Plot name"), value = "Anpu_hm"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_apu_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_apu_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_apu <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in A. punctatus"),
  h5("Please cite ***"),
  shiny::plotOutput("hmPlot_apu"),
  downloadButton(outputId = "hm_down_apu", label = "Download plot"),
  downloadButton("download_apu_clust", "Download result table"),
  div(shiny::tableOutput('hm_tbl_apu'), style = "font-size:80%")
)


hm_panel_apu <- tabPanel(
  "A. punctatus heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_apu, main_content_hm_apu
  )
)

### Page setup D. glomerata sample heatmap

sidebar_content_hm_dgl <- sidebarPanel(
  shiny::textAreaInput("genes_hm_dgl", 'Input gene names separated by new line:', 
                       value = paste("DgTrNR04516_a1_i1", "DgTrNR00454_a1_i1", "DgTrNR05876_a1_i2", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_dgl", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_dgl_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("dgl_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_dgl", label = h5("Plot name"), value = "Dagl_hm"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_dgl_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_dgl_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_dgl <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in D. glomerata"),
  h5("Please cite ***"),
  shiny::plotOutput("hmPlot_dgl"),
  downloadButton(outputId = "hm_down_dgl", label = "Download plot"),
  downloadButton("download_dgl_clust", "Download result table"),
  div(shiny::tableOutput('hm_tbl_dgl'), style = "font-size:80%")
)


hm_panel_dgl <- tabPanel(
  "D. glomerata heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_dgl, main_content_hm_dgl
  )
)
### Page setup G. perpensa sample heatmap

sidebar_content_hm_gpe <- sidebarPanel(
  shiny::textAreaInput("genes_hm_gpe", 'Input gene names separated by new line:', 
                       value = paste("G_perpensa_evgTRINITY_DN21508_c5_g7", "G_perpensa_evgd5483TRINITY_DN19016_c0_g1", "G_perpensa_evgTRINITY_DN393_c0_g1", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_gpe", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_gpe_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("gpe_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_gpe", label = h5("Plot name"), value = "Gupe_hm"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_gpe_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_gpe_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_gpe <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in G. perpensa"),
  h5("Please cite ***"),
  shiny::plotOutput("hmPlot_gpe"),
  downloadButton(outputId = "hm_down_gpe", label = "Download plot"),
  downloadButton("download_gpe_clust", "Download result table"),
  div(shiny::tableOutput('hm_tbl_gpe'), style = "font-size:80%")
)


hm_panel_gpe <- tabPanel(
  "G. perpensa heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_gpe, main_content_hm_gpe
  )
)

### Page setup M. truncatula sample heatmap

sidebar_content_hm_mtr <- sidebarPanel(
  shiny::textAreaInput("genes_hm_mtr", 'Input gene names separated by new line:', 
                       value = paste("MtrunA17Chr3g0100391", "MtrunA17Chr4g0070011", "MtrunA17Chr1g0157221", sep = "\n"), 
                       width = NULL),
  submitButton("Update plot"),
  shiny::sliderInput("k_mtr", "K-means clustering", min = 2, 
                     max = 10, value = 2),
  div(style="display: inline-block;vertical-align:top;",shiny::numericInput("hm_mtr_percentile", label = h5("Min % rank"),1, min = 1, 
                                                                            max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 10px;",HTML("<br>")),
  div(style="display: inline-block;vertical-align:top;", shiny::selectInput("mtr_hm_plot_row_names_show", label = h5("Show row names"), 
                                                                            choices = list("Yes" = TRUE, "No" = FALSE), 
                                                                            selected = TRUE)),
  p("Plot download options"),
  textInput("hm_title_mtr", label = h5("Plot name"), value = "Metr_hm"),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_mtr_height", "Height:", 5, min = 1, max = 100)),
  div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("hm_mtr_width", "Width:", 7, min = 1, max = 100))
)  

main_content_hm_mtr <- mainPanel(
  p("Plot a k-means clustered heatmap of tissue-specific gene expression in M. truncatula"),
  h5("Please cite ***"),
  shiny::plotOutput("hmPlot_mtr"),
  downloadButton(outputId = "hm_down_mtr", label = "Download plot"),
  downloadButton("download_mtr_clust", "Download result table"),
  div(shiny::tableOutput('hm_tbl_mtr'), style = "font-size:80%")
)


hm_panel_mtr <- tabPanel(
  "M. truncatula heatmap clustering",
  sidebarLayout(
    sidebar_content_hm_mtr, main_content_hm_mtr
  )
)


ui <- navbarPage("Nitrogen symbiosis", theme = shinytheme("flatly"),
                 navbarMenu("E. natalensis",
                            tabPanel(title = "Bar chart", ena_panel_tissue),
                            tabPanel(title = "Heatmap", hm_panel_tissue)
                            ),
                  navbarMenu("Orthogroup",
                             tabPanel(title = "Infected vs uninfected bar chart", bar_panel_ortho),
                             tabPanel(title = "E. natalensis tissue bar chart", bar_panel_ortho_ena)
                             ),
                  navbarMenu("Convert",
                             tabPanel(title = "Convert gene IDs", convert_panel)
                             ),
                 navbarMenu("Species heatmaps",
                            tabPanel(title = "Azolla filiculoides", hm_panel_afi),
                            tabPanel(title = "Anthoceros punctatus", hm_panel_apu),
                            tabPanel(title = "Dactylis glomerata", hm_panel_dgl),
                            tabPanel(title = "Gunnera perpensa", hm_panel_gpe),
                            tabPanel(title = "Medicago truncatula", hm_panel_mtr))
                 
)

server <- function(input, output, session) {  
  vals <- reactiveValues()

### Server_bar_ortho
  
  output$barPlot_ortho <- renderPlot({
    ortho_bar_id <- toupper(input$text_bar_ortho)
    validate(
      need(ortho_bar_id %in% ortho_bar_data$Orthogroup, "Orthogroup not in dataset")
    )
    ### Set colors per species as defined in the data file column "cols" in ortho_bar_data
    bar_cols <- unique(ortho_bar_data[ , c(10,9)])
    my_cols <- bar_cols$cols
    names(my_cols) <- bar_cols$cols_id
    bar_cols_scale <- scale_fill_manual(name = "cols_id", values = my_cols) 
    bar_data <- ortho_bar_data[ortho_bar_data$Orthogroup == ortho_bar_id, ]
    bar_data <- bar_data %>%
      mutate(x.label = paste("<span style = ",
                             ";'>",
                             ifelse(is.na(signif), "", "**"),
                             id_lab,
                             ifelse(is.na(signif), "", "**"),
                             "</span>", sep = ""),
      )
    bar_data$species <- factor(bar_data$species,levels=c("An_pun","Az_fil","En_nat","Gu_per","Me_tru","Da_glo"))
    ortho_bar_gg <- ggplot(bar_data, aes(x=x.label, y=mean, fill=cols_id)) + 
      geom_bar(stat="identity", 
               position="dodge", color = "black") +
      bar_cols_scale +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.25,
                    position=position_dodge(.9)) +
      labs(x="", y = "Expression level") +
      theme_classic() +
      theme(legend.position="none") +
      facet_grid(.~species, scales = "free", space = "free") +
      theme(axis.text.x = element_markdown(angle = 90, hjust = 1)) + 
      theme(axis.text=element_text(size=12)) +
      theme(legend.position = "none") +
      ggtitle(ortho_bar_id) +
      theme(plot.title = element_text(face = "bold",
                                      margin = margin(10, 0, 10, 0),
                                      size = 14))
    vals$ortho_bar_gg_plot <- ortho_bar_gg
    print(ortho_bar_gg)
  })
  output$tbl_bar_ortho <- renderTable({
    tbl_out_bar_ortho <- unique(ortho_annot[ortho_annot$Orthogroup == toupper(input$text_bar_ortho),])
    tbl_out_bar_ortho
  })
  output$down_bar_ortho <- downloadHandler(
    filename = function(){paste(input$text_bar_ortho, 'ortho_bar.pdf', sep = '_')},
    content = function(file){
      pdf(file, width = input$bar_ortho_width, height = input$bar_ortho_height)
      print(vals$ortho_bar_gg_plot)
      dev.off()
    })

### Server_bar_ortho_ena  
  
  output$barPlot_ortho_ena <- renderPlot({
    ortho_bar_id_ena <- toupper(input$text_bar_ortho_ena)
    validate(
      need(ortho_bar_id_ena %in% ena_tissue_bar_melt$Orthogroup, "Orthogroup not in dataset")
      )
      bar_data_ena <- ena_tissue_bar_melt[ena_tissue_bar_melt$Orthogroup == ortho_bar_id_ena, ]
      bar_data_ena$Tissue <- factor(bar_data_ena$Tissue,levels=c("CR","UPR","LPR","ST","IR","MR","IL","ML","GS"))
      ortho_bar_gg_ena <- ggplot(bar_data_ena, aes(x=Tissue, y=mean, fill=Tissue)) + 
        geom_bar(stat="identity", 
                 position="dodge", color = "black") +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.25,
                      position=position_dodge(.9)) +
        labs(x="", y = "Expression level") +
        theme_classic() +
        ena_bar_cols_tissue +
        facet_wrap(.~ena_id) +
        ggtitle(ortho_bar_id_ena) +
        theme(plot.title = element_text(face = "bold",
                                        margin = margin(10, 0, 10, 0),
                                        size = 14))
        theme(axis.text=element_text(size=12))
      vals$ortho_bar_gg_plot_ena <- ortho_bar_gg_ena
      print(ortho_bar_gg_ena)
  })
    output$tbl_bar_ortho_ena <- renderTable({
      tbl_out_bar_ortho_ena <- unique(ortho_annot[ortho_annot$Orthogroup == toupper(input$text_bar_ortho_ena),])
      tbl_out_bar_ortho_ena
    })
    output$down_bar_ortho_ena <- downloadHandler(
      filename = function(){paste(input$text_bar_ortho, 'ortho_bar_ena.pdf', sep = '_')},
      content = function(file){
        pdf(file, width = 12, height = 4)
        print(vals$ortho_bar_gg_plot_ena)
        dev.off()
    
  })
    
### Server E. natalensis heatmap
    output$hmPlot_tissue <- renderPlot({
      tissue_hm_plot_row_names <- input$tissue_hm_plot_row_names_show
      input_ids_tissue <- gsub(" ", "\n",input$genes_hm_tissue)
      input_ids_tissue <- gsub(", ", "\n",input$genes_hm_tissue)
      input_ids_split_tissue <- strsplit(input_ids_tissue, "\n")
      input_ids_split_tissue <- as.data.frame(input_ids_split_tissue)
      validate(
        need(nrow(input_ids_split_tissue) > 1, "Needs more than one gene ID, did you check that your input is separated by a newline or pasted from the convert clipboard?")
      )
      colnames(input_ids_split_tissue) <- c("ena_id")
      data_hm_tissue <- ena_hm_mat[row.names(ena_hm_mat) %in% input_ids_split_tissue$ena_id, ]
      percentile_filter_tissue <- tissue_data[tissue_data$percentile_rank >= input$hm_tissue_percentile, ]
      data_hm_tissue <- data_hm_tissue[row.names(data_hm_tissue) %in% percentile_filter_tissue$ena_id, ]
      data_hm_tissue <- unique.matrix(data_hm_tissue)
      val_k_ena_tissue <- input$k_ena_tissue
      validate(
        need((nrow(input_ids_split_tissue) - 1) >= val_k_ena_tissue, "Not enough genes for k-value - try a lower k or more genes"))
      hca_ena_tissue <- hclust(dist(data_hm_tissue))
      clust_ena_tissue <- cutree(hca_ena_tissue,k=val_k_ena_tissue,order_clusters_as_data = FALSE)
      clust_IDs_ena_tissue <- as.data.frame(clust_ena_tissue)
      clust_IDs_ena_tissue$ena_id <- names(clust_ena_tissue)
      dend_ena_tissue <-as.dendrogram(hca_ena_tissue)
      dend1_ena_tissue <- color_branches(dend_ena_tissue, k = val_k_ena_tissue, groupLabels = TRUE)
      row_dend_ena_tissue <- dend1_ena_tissue
      genes_hm_tissue_percentile <- tissue_data[tissue_data$ena_id %in% row.names(data_hm_tissue), ]
      genes_hm_tissue_percentile <- genes_hm_tissue_percentile[!duplicated(genes_hm_tissue_percentile$ena_id), ]
      tissue_row_annot <- rowAnnotation(Percentile = genes_hm_tissue_percentile$percentile_rank, 
                                        col = list(Percentile = per_row_col), na_col = "white", 
                                        annotation_label = "%",
                                        annotation_legend_param = list(Percentile = 
                                                                         list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
      hm_gg_tissue_plot <- Heatmap(data_hm_tissue, name = "Scaled VST", top_annotation = ena_tissue_col_annot, right_annotation = tissue_row_annot,
                                   column_title = "Tissue-specific E. natalensis heatmap",
                                   show_row_names = tissue_hm_plot_row_names, show_column_names = TRUE, column_order = c("CR", "UPR", "LPR", "ST", "IR","MR", "IL", "ML", "GS"),
                                   show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_ena_tissue, row_split = as.numeric(val_k_ena_tissue),
                                   heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
      hm_gg_tissue <- draw(hm_gg_tissue_plot, heatmap_legend_side = "bottom", annotation_legend_side = "right", merge_legend = F)
      vals$hm_gg_tissue <- hm_gg_tissue
      print(hm_gg_tissue)
    })
    output$hm_down_tissue <- downloadHandler(
      filename = function(){paste(input$hm_title_tissue, 'tissue.pdf', sep = '_')},
      content = function(file){
        pdf(file, width = input$hm_tissue_width, height = input$hm_tissue_height)
        print(vals$hm_gg_tissue)
        dev.off()
      })
    output$hm_tbl_tissue <- renderTable({
      input_ids_tissue <- gsub(" ", "\n",input$genes_hm_tissue)
      input_ids_tissue <- gsub(", ", "\n",input$genes_hm_tissue)
      input_ids_split_tissue <- strsplit(input_ids_tissue, "\n")
      input_ids_split_tissue <- as.data.frame(input_ids_split_tissue)
      validate(
        need(nrow(input_ids_split_tissue) > 1, "Needs more than one gene ID, did you check that your input is separated by a newline or pasted from the convert clipboard?")
      )
      colnames(input_ids_split_tissue) <- c("ena_id")
      data_hm_tissue <- ena_hm_mat[row.names(ena_hm_mat) %in% input_ids_split_tissue$ena_id, ]
      percentile_filter_tissue <- tissue_data[tissue_data$percentile >= input$hm_tissue_percentile, ]
      data_hm_tissue <- data_hm_tissue[row.names(data_hm_tissue) %in% percentile_filter_tissue$ena_id, ]
      val_k_ena_tissue <- input$k_ena_tissue
      hca_ena_tissue <- hclust(dist(data_hm_tissue))
      clust_ena_tissue <- cutree(hca_ena_tissue,k=val_k_ena_tissue,order_clusters_as_data = FALSE)
      clust_IDs_ena_tissue <- as.data.frame(clust_ena_tissue)
      clust_IDs_ena_tissue$ena_id <- names(clust_ena_tissue)
      dend_ena_tissue <-as.dendrogram(hca_ena_tissue)
      dend1_ena_tissue <- color_branches(dend_ena_tissue, k = val_k_ena_tissue, groupLabels = TRUE)
      row_dend_ena_tissue <- dend1_ena_tissue
      table_tissue <- ena_annot[ena_annot$ena_id %in% input_ids_split_tissue$ena_id, ]
      tissue_table <- merge(table_tissue, clust_IDs_ena_tissue, by = "ena_id")
      tissue_table <- merge(tissue_table, tissue_data[ , c(2,3,14)], by = "ena_id")
      tissue_table <- unique(tissue_table)
      tissue_table
      vals$hm_gg_clust_table_tissue <- tissue_table
    })
    output$download_tissue_clust <- downloadHandler(
      filename = function() {
        paste(input$hm_title_tissue, 'cluster_tissue.csv', sep = '_')
      },
      content = function(file) {
        write.csv(vals$hm_gg_clust_table_tissue, file, row.names = FALSE)
      }
    )
    ### Server - convert
    
    output$convert_tbl <- renderTable({  
      validate(
        need(input$convert_from != input$convert_to, "Remember to change convert from and convert to species - they cannot be the same")
      )
      if (input$convert_from == "At_tha" && input$convert_to =='En_nat'){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("At_tha") 
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "At_tha")
        convert_tbl <- convert_tbl[ ,c(3,1,2)]
      } 
      else if (input$convert_from == "At_tha" && input$convert_to =='Me_tru'){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("At_tha") 
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "At_tha")
        convert_tbl <- convert_tbl[ ,c(7,1,2)]
      } 
      else if (input$convert_from == "Me_tru" && input$convert_to =='At_tha'){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("Me_tru")
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "Me_tru")
        convert_tbl <- convert_tbl[ ,c(7,1,2)]
      }
      else if (input$convert_from == "En_nat" && input$convert_to =='At_tha'){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("En_nat")
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "En_nat")
        convert_tbl <- convert_tbl[ ,c(7,1,2)]
      }
      else if (input$convert_from == "En_nat" && input$convert_to =='Me_tru'){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("En_nat")
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "En_nat")
        convert_tbl <- convert_tbl[ ,c(6,1,2)]
      }
      else if (input$convert_from == "Me_tru" && input$convert_to =="En_nat"){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("Me_tru")
        convert_tbl <- merge(ena_best_hits, input_ids_split, by = "Me_tru")
        convert_tbl <- convert_tbl[ ,c(3,1,2)]
      }
      else if (input$convert_from == "Any" && input$convert_to =="En_nat"){
        input_ids <- input$convert_genes
        input_ids_split <- strsplit(input_ids, "\n")
        input_ids_split <- as.data.frame(input_ids_split)
        colnames(input_ids_split) <- c("id")
        convert_tbl <- merge(ena_best_any, input_ids_split, by.x = "id")
        convert_tbl <- convert_tbl[ ,c(3,1,2,4)]
      }
      vals$convert_tbl <- convert_tbl
      return(convert_tbl)
      
    })
    output$clip <- renderUI({
      values <- toString(vals$convert_tbl[,1])
      values <- gsub(", ", "\n",values)
      rclipButton("clipbtn", "Copy", values, icon("clipboard"))
    })
    output$download_converter <- downloadHandler(
      filename = function() {
        paste(input$name_converter, 'convert.csv', sep = '_')
      },
      content = function(file) {
        write.csv(vals$convert_tbl, file, row.names = FALSE)
      }
    )
    
    ### Server_bar_tissue_ena  
    
    output$barPlot_tissue <- renderPlot({
      bar_id_ena <- input$text_bar_tissue_ena
      validate(
        need(bar_id_ena %in% ena_tissue_bar_melt$ena_id, "Gene ID not in dataset")
      )
      bar_data_ena <- ena_tissue_bar_melt[ena_tissue_bar_melt$ena_id == input$text_bar_tissue_ena, ]
      bar_data_ena$Tissue <- factor(bar_data_ena$Tissue,levels=c("CR","UPR","LPR","ST","IR","MR","IL", "ML", "GS"))
      tissue_bar_gg_ena <- ggplot(bar_data_ena, aes(x=Tissue, y=mean, fill=Tissue)) + 
        geom_bar(stat="identity", 
                 position="dodge", color = "black") +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.25,
                      position=position_dodge(.9)) +
        labs(x="Tissue", y = "Expression level") +
        theme_classic() +
        ena_bar_cols_tissue +
        theme(axis.text=element_text(size=12)) +
        theme(legend.position = "none") +
        ggtitle(input$text_bar_tissue_ena) +
        theme(plot.title = element_text(face = "bold",
                                        margin = margin(10, 0, 10, 0),
                                        size = 14))
      vals$tissue_bar <- tissue_bar_gg_ena
      print(vals$tissue_bar)
    })
    output$tbl_bar_tissue <- renderTable({
      tbl_out_bar_ena <- unique(ena_annot[ena_annot$ena_id == input$text_bar_tissue_ena, ])
      tbl_out_bar_tissue <- merge(tbl_out_bar_ena, tissue_data[ , c(2,3,14)], by = "ena_id")
      tbl_out_bar_tissue
    })
    output$down_bar_tissue <- downloadHandler(
      filename = function(){paste(input$text_bar_tissue_ena, 'tissue_bar_ena.pdf', sep = '_')},
      content = function(file){
        pdf(file, width = 12, height = 4)
        print(vals$tissue_bar)
        dev.off()
      })
    
    
    ### Server "A. filiculoides heatmap
    output$hmPlot_afi <- renderPlot({
      afi_hm_plot_row_names <- input$afi_hm_plot_row_names_show
      input_ids_afi <- gsub(" ", "\n",input$genes_hm_afi)
      input_ids_afi <- gsub(", ", "\n",input$genes_hm_afi)
      input_ids_split_afi <- strsplit(input_ids_afi, "\n")
      input_ids_split_afi <- as.data.frame(input_ids_split_afi)
      validate(
        need(nrow(input_ids_split_afi) > 1, "Needs more than one gene ID")
      )
      colnames(input_ids_split_afi) <- c("afi_id")
      data_hm_afi <- afi_hm_mat[row.names(afi_hm_mat) %in% input_ids_split_afi$afi_id, ]
      percentile_filter_afi <- afi_sample_data[afi_sample_data$percentile_rank >= input$hm_afi_percentile, ]
      data_hm_afi <- data_hm_afi[row.names(data_hm_afi) %in% percentile_filter_afi$afi_id, ]
      data_hm_afi <- unique.matrix(data_hm_afi)
      val_k_afi <- input$k_afi
      validate(
        need((nrow(input_ids_split_afi) - 1) >= val_k_afi, "Not enough genes for k-value - try a lower k or more genes"))
      hca_afi <- hclust(dist(data_hm_afi))
      clust_afi <- cutree(hca_afi,k=val_k_afi,order_clusters_as_data = FALSE)
      clust_IDs_afi <- as.data.frame(clust_afi)
      clust_IDs_afi$afi_id <- names(clust_afi)
      dend_afi <-as.dendrogram(hca_afi)
      dend1_afi <- color_branches(dend_afi, k = val_k_afi, groupLabels = TRUE)
      row_dend_afi <- dend1_afi
      genes_hm_afi_percentile <- afi_sample_data[afi_sample_data$afi_id %in% row.names(data_hm_afi), ]
      genes_hm_afi_percentile <- genes_hm_afi_percentile[!duplicated(genes_hm_afi_percentile$afi_id), ]
      afi_row_annot <- rowAnnotation(Percentile = genes_hm_afi_percentile$percentile_rank, 
                                     col = list(Percentile = per_row_col), na_col = "white", 
                                     annotation_label = "%",
                                     annotation_legend_param = list(Percentile = 
                                                                      list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
      hm_gg_afi_plot <- Heatmap(data_hm_afi, name = "Scaled expression", top_annotation = afi_annot_col, right_annotation = afi_row_annot,
                                column_title = "A. filiculoides heatmap",
                                show_row_names = afi_hm_plot_row_names, show_column_names = TRUE, column_order = c("NM_CP_1", "NM_CP_2", "NM_CP_3", "NP_CP_1", "NP_CP_2","NP_CP_3","NP_CM_1", "NP_CM_2", "NP_CM_3"),
                                show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_afi, row_split = as.numeric(val_k_afi),
                                heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
      hm_gg_afi <- draw(hm_gg_afi_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
      vals$hm_gg_afi <- hm_gg_afi
      out_table_afi <- unique(afi_annot[afi_annot$afi_id %in% input_ids_split_afi$afi_id, ])
      out_table_afi <- merge(out_table_afi, afi_sample_data[ , c(1,12)], by = "afi_id")
      out_table_afi <- out_table_afi[!is.na(out_table_afi$afi_id), ]
      out_table_afi <- unique(out_table_afi)
      out_table_afi <- merge(out_table_afi, clust_IDs_afi, by = "afi_id")
      vals$hm_gg_clust_table_afi <- out_table_afi
      print(hm_gg_afi)
    })
    output$hm_down_afi <- downloadHandler(
      filename = function(){paste(input$hm_title_afi, '.pdf', sep = '_')},
      content = function(file){
        pdf(file, width = input$hm_afi_width, height = input$hm_afi_height)
        print(vals$hm_gg_afi)
        dev.off()
      })
    output$hm_tbl_afi <- renderTable({
      out_table_afi_clust <- vals$hm_gg_clust_table_afi
      out_table_afi_clust
      })
    output$download_afi_clust <- downloadHandler(
      filename = function() {
        paste(input$hm_title_afi, 'cluster_afi.csv', sep = '_')
      },
      content = function(file) {
        write.csv(vals$hm_gg_clust_table_afi, file, row.names = FALSE)
      }
    )
    ### Server A. punctatus heatmap
    output$hmPlot_apu <- renderPlot({
      apu_hm_plot_row_names <- input$apu_hm_plot_row_names_show
      input_ids_apu <- gsub(" ", "\n",input$genes_hm_apu)
      input_ids_apu <- gsub(", ", "\n",input$genes_hm_apu)
      input_ids_split_apu <- strsplit(input_ids_apu, "\n")
      input_ids_split_apu <- as.data.frame(input_ids_split_apu)
      validate(
        need(nrow(input_ids_split_apu) > 1, "Needs more than one gene ID")
      )
      colnames(input_ids_split_apu) <- c("apu_id")
      data_hm_apu <- apu_hm_mat[row.names(apu_hm_mat) %in% input_ids_split_apu$apu_id, ]
      percentile_filter_apu <- apu_sample_data[apu_sample_data$percentile_rank >= input$hm_apu_percentile, ]
      data_hm_apu <- data_hm_apu[row.names(data_hm_apu) %in% percentile_filter_apu$apu_id, ]
      data_hm_apu <- unique.matrix(data_hm_apu)
      val_k_apu <- input$k_apu
      validate(
        need((nrow(input_ids_split_apu) - 1) >= val_k_apu, "Not enough genes for k-value - try a lower k or more genes"))
      hca_apu <- hclust(dist(data_hm_apu))
      clust_apu <- cutree(hca_apu,k=val_k_apu,order_clusters_as_data = FALSE)
      clust_IDs_apu <- as.data.frame(clust_apu)
      clust_IDs_apu$apu_id <- names(clust_apu)
      dend_apu <-as.dendrogram(hca_apu)
      dend1_apu <- color_branches(dend_apu, k = val_k_apu, groupLabels = TRUE)
      row_dend_apu <- dend1_apu
      genes_hm_apu_percentile <- apu_sample_data[apu_sample_data$apu_id %in% row.names(data_hm_apu), ]
      genes_hm_apu_percentile <- genes_hm_apu_percentile[!duplicated(genes_hm_apu_percentile$apu_id), ]
      apu_row_annot <- rowAnnotation(Percentile = genes_hm_apu_percentile$percentile_rank, 
                                     col = list(Percentile = per_row_col), na_col = "white", 
                                     annotation_label = "%",
                                     annotation_legend_param = list(Percentile = 
                                                                      list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
      hm_gg_apu_plot <- Heatmap(data_hm_apu, name = "Scaled expression", top_annotation = apu_annot_col, right_annotation = apu_row_annot,
                                column_title = "A. punctatus heatmap",
                                show_row_names = apu_hm_plot_row_names, show_column_names = TRUE, column_order = c("NM_CP_1", "NM_CP_2", "NM_CP_3", "NP_CM_1", "NP_CM_2", "NP_CM_3"),
                                show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_apu, row_split = as.numeric(val_k_apu),
                                heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
      hm_gg_apu <- draw(hm_gg_apu_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
      vals$hm_gg_apu <- hm_gg_apu
      out_table_apu <- unique(apu_annot[apu_annot$apu_id %in% input_ids_split_apu$apu_id, ])
      out_table_apu <- merge(out_table_apu, apu_sample_data[ , c(1,9)], by = "apu_id")
      out_table_apu <- out_table_apu[!is.na(out_table_apu$apu_id), ]
      out_table_apu <- unique(out_table_apu)
      out_table_apu <- merge(out_table_apu, clust_IDs_apu, by = "apu_id")
      vals$hm_gg_clust_table_apu <- out_table_apu
      print(hm_gg_apu)
    })
      output$hm_down_apu <- downloadHandler(
        filename = function(){paste(input$hm_title_apu, '.pdf', sep = '_')},
        content = function(file){
          pdf(file, width = input$hm_apu_width, height = input$hm_apu_height)
          print(vals$hm_gg_apu)
          dev.off()
        })
      output$hm_tbl_apu <- renderTable({
        out_table_apu_clust <- vals$hm_gg_clust_table_apu
        out_table_apu_clust
      })
      output$download_apu_clust <- downloadHandler(
        filename = function() {
          paste(input$hm_title_apu, 'cluster_apu.csv', sep = '_')
        },
        content = function(file) {
          write.csv(vals$hm_gg_clust_table_apu, file, row.names = FALSE)
        }
      )
      ### Server D. glomerata heatmap
      output$hmPlot_dgl <- renderPlot({
        dgl_hm_plot_row_names <- input$dgl_hm_plot_row_names_show
        input_ids_dgl <- gsub(" ", "\n",input$genes_hm_dgl)
        input_ids_dgl <- gsub(", ", "\n",input$genes_hm_dgl)
        input_ids_split_dgl <- strsplit(input_ids_dgl, "\n")
        input_ids_split_dgl <- as.data.frame(input_ids_split_dgl)
        validate(
          need(nrow(input_ids_split_dgl) > 1, "Needs more than one gene ID")
        )
        colnames(input_ids_split_dgl) <- c("dgl_id")
        data_hm_dgl <- dgl_hm_mat[row.names(dgl_hm_mat) %in% input_ids_split_dgl$dgl_id, ]
        percentile_filter_dgl <- dgl_sample_data[dgl_sample_data$percentile_rank >= input$hm_dgl_percentile, ]
        data_hm_dgl <- data_hm_dgl[row.names(data_hm_dgl) %in% percentile_filter_dgl$dgl_id, ]
        data_hm_dgl <- unique.matrix(data_hm_dgl)
        val_k_dgl <- input$k_dgl
        validate(
          need((nrow(input_ids_split_dgl) - 1) >= val_k_dgl, "Not enough genes for k-value - try a lower k or more genes"))
        hca_dgl <- hclust(dist(data_hm_dgl))
        clust_dgl <- cutree(hca_dgl,k=val_k_dgl,order_clusters_as_data = FALSE)
        clust_IDs_dgl <- as.data.frame(clust_dgl)
        clust_IDs_dgl$dgl_id <- names(clust_dgl)
        dend_dgl <-as.dendrogram(hca_dgl)
        dend1_dgl <- color_branches(dend_dgl, k = val_k_dgl, groupLabels = TRUE)
        row_dend_dgl <- dend1_dgl
        genes_hm_dgl_percentile <- dgl_sample_data[dgl_sample_data$dgl_id %in% row.names(data_hm_dgl), ]
        genes_hm_dgl_percentile <- genes_hm_dgl_percentile[!duplicated(genes_hm_dgl_percentile$dgl_id), ]
        dgl_row_annot <- rowAnnotation(Percentile = genes_hm_dgl_percentile$percentile_rank, 
                                       col = list(Percentile = per_row_col), na_col = "white", 
                                       annotation_label = "%",
                                       annotation_legend_param = list(Percentile = 
                                                                        list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
        hm_gg_dgl_plot <- Heatmap(data_hm_dgl, name = "Scaled expression", top_annotation = dgl_annot_col, right_annotation = dgl_row_annot,
                                  column_title = "D. glomerata heatmap",
                                  show_row_names = dgl_hm_plot_row_names, show_column_names = TRUE, column_order = c("nodule_1", "nodule_2", "nodule_3", "nodule_4", "nodule_5", "nodule_6", "root_1", "root_2", "root_3", "root_4", "root_5", "root_6"),
                                  show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_dgl, row_split = as.numeric(val_k_dgl),
                                  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
        hm_gg_dgl <- draw(hm_gg_dgl_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
        vals$hm_gg_dgl <- hm_gg_dgl
        out_table_dgl <- unique(dgl_annot[dgl_annot$dgl_id %in% input_ids_split_dgl$dgl_id, ])
        out_table_dgl <- merge(out_table_dgl, dgl_sample_data[ , c(1,15)], by = "dgl_id")
        out_table_dgl <- out_table_dgl[!is.na(out_table_dgl$dgl_id), ]
        out_table_dgl <- unique(out_table_dgl)
        out_table_dgl <- merge(out_table_dgl, clust_IDs_dgl, by = "dgl_id")
        vals$hm_gg_clust_table_dgl <- out_table_dgl
        print(hm_gg_dgl)
      })
      output$hm_down_dgl <- downloadHandler(
        filename = function(){paste(input$hm_title_dgl, '.pdf', sep = '_')},
        content = function(file){
          pdf(file, width = input$hm_dgl_width, height = input$hm_dgl_height)
          print(vals$hm_gg_dgl)
          dev.off()
        })
      output$hm_tbl_dgl <- renderTable({
        out_table_dgl_clust <- vals$hm_gg_clust_table_dgl
        out_table_dgl_clust
      })
      output$download_dgl_clust <- downloadHandler(
        filename = function() {
          paste(input$hm_title_dgl, 'cluster_dgl.csv', sep = '_')
        },
        content = function(file) {
          write.csv(vals$hm_gg_clust_table_dgl, file, row.names = FALSE)
        }
      )
      
      ### Server G. perpensa heatmap
      output$hmPlot_gpe <- renderPlot({
        gpe_hm_plot_row_names <- input$gpe_hm_plot_row_names_show
        input_ids_gpe <- gsub(" ", "\n",input$genes_hm_gpe)
        input_ids_gpe <- gsub(", ", "\n",input$genes_hm_gpe)
        input_ids_split_gpe <- strsplit(input_ids_gpe, "\n")
        input_ids_split_gpe <- as.data.frame(input_ids_split_gpe)
        validate(
          need(nrow(input_ids_split_gpe) > 1, "Needs more than one gene ID")
        )
        colnames(input_ids_split_gpe) <- c("gpe_id")
        data_hm_gpe <- gpe_hm_mat[row.names(gpe_hm_mat) %in% input_ids_split_gpe$gpe_id, ]
        percentile_filter_gpe <- gpe_sample_data[gpe_sample_data$percentile_rank >= input$hm_gpe_percentile, ]
        data_hm_gpe <- data_hm_gpe[row.names(data_hm_gpe) %in% percentile_filter_gpe$gpe_id, ]
        data_hm_gpe <- unique.matrix(data_hm_gpe)
        val_k_gpe <- input$k_gpe
        validate(
          need((nrow(input_ids_split_gpe) - 1) >= val_k_gpe, "Not enough genes for k-value - try a lower k or more genes"))
        hca_gpe <- hclust(dist(data_hm_gpe))
        clust_gpe <- cutree(hca_gpe,k=val_k_gpe,order_clusters_as_data = FALSE)
        clust_IDs_gpe <- as.data.frame(clust_gpe)
        clust_IDs_gpe$gpe_id <- names(clust_gpe)
        dend_gpe <-as.dendrogram(hca_gpe)
        dend1_gpe <- color_branches(dend_gpe, k = val_k_gpe, groupLabels = TRUE)
        row_dend_gpe <- dend1_gpe
        genes_hm_gpe_percentile <- gpe_sample_data[gpe_sample_data$gpe_id %in% row.names(data_hm_gpe), ]
        genes_hm_gpe_percentile <- genes_hm_gpe_percentile[!duplicated(genes_hm_gpe_percentile$gpe_id), ]
        gpe_row_annot <- rowAnnotation(Percentile = genes_hm_gpe_percentile$percentile_rank, 
                                       col = list(Percentile = per_row_col), na_col = "white", 
                                       annotation_label = "%",
                                       annotation_legend_param = list(Percentile = 
                                                                        list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
        hm_gg_gpe_plot <- Heatmap(data_hm_gpe, name = "Scaled expression", top_annotation = gpe_annot_col, right_annotation = gpe_row_annot,
                                  column_title = "G. perpensa heatmap",
                                  show_row_names = gpe_hm_plot_row_names, show_column_names = TRUE, column_order = c("Cyano+ 1", "Cyano+ 2", "Cyano+ 3", "Control 1", "Control 2", "Control 3"),
                                  show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_gpe, row_split = as.numeric(val_k_gpe),
                                  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
        hm_gg_gpe <- draw(hm_gg_gpe_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
        vals$hm_gg_gpe <- hm_gg_gpe
        out_table_gpe <- unique(gpe_annot[gpe_annot$gpe_id %in% input_ids_split_gpe$gpe_id, ])
        out_table_gpe <- merge(out_table_gpe, gpe_sample_data[ , c(1,9)], by = "gpe_id")
        out_table_gpe <- out_table_gpe[!is.na(out_table_gpe$gpe_id), ]
        out_table_gpe <- unique(out_table_gpe)
        out_table_gpe <- merge(out_table_gpe, clust_IDs_gpe, by = "gpe_id")
        vals$hm_gg_clust_table_gpe <- out_table_gpe
        print(hm_gg_gpe)
      })
      output$hm_down_gpe <- downloadHandler(
        filename = function(){paste(input$hm_title_gpe, '.pdf', sep = '_')},
        content = function(file){
          pdf(file, width = input$hm_gpe_width, height = input$hm_gpe_height)
          print(vals$hm_gg_gpe)
          dev.off()
        })
      output$hm_tbl_gpe <- renderTable({
        out_table_gpe_clust <- vals$hm_gg_clust_table_gpe
        out_table_gpe_clust
      })
      output$download_gpe_clust <- downloadHandler(
        filename = function() {
          paste(input$hm_title_gpe, 'cluster_gpe.csv', sep = '_')
        },
        content = function(file) {
          write.csv(vals$hm_gg_clust_table_gpe, file, row.names = FALSE)
        }
      )
      ### Server M. truncatula heatmap
      output$hmPlot_mtr <- renderPlot({
        mtr_hm_plot_row_names <- input$mtr_hm_plot_row_names_show
        input_ids_mtr <- gsub(" ", "\n",input$genes_hm_mtr)
        input_ids_mtr <- gsub(", ", "\n",input$genes_hm_mtr)
        input_ids_split_mtr <- strsplit(input_ids_mtr, "\n")
        input_ids_split_mtr <- as.data.frame(input_ids_split_mtr)
        validate(
          need(nrow(input_ids_split_mtr) > 1, "Needs more than one gene ID")
        )
        colnames(input_ids_split_mtr) <- c("mtr_id")
        data_hm_mtr <- mtr_hm_mat[row.names(mtr_hm_mat) %in% input_ids_split_mtr$mtr_id, ]
        percentile_filter_mtr <- mtr_sample_data[mtr_sample_data$percentile_rank >= input$hm_mtr_percentile, ]
        data_hm_mtr <- data_hm_mtr[row.names(data_hm_mtr) %in% percentile_filter_mtr$mtr_id, ]
        data_hm_mtr <- unique.matrix(data_hm_mtr)
        val_k_mtr <- input$k_mtr
        validate(
          need((nrow(input_ids_split_mtr) - 1) >= val_k_mtr, "Not enough genes for k-value - try a lower k or more genes"))
        hca_mtr <- hclust(dist(data_hm_mtr))
        clust_mtr <- cutree(hca_mtr,k=val_k_mtr,order_clusters_as_data = FALSE)
        clust_IDs_mtr <- as.data.frame(clust_mtr)
        clust_IDs_mtr$mtr_id <- names(clust_mtr)
        dend_mtr <-as.dendrogram(hca_mtr)
        dend1_mtr <- color_branches(dend_mtr, k = val_k_mtr, groupLabels = TRUE)
        row_dend_mtr <- dend1_mtr
        genes_hm_mtr_percentile <- mtr_sample_data[mtr_sample_data$mtr_id %in% row.names(data_hm_mtr), ]
        genes_hm_mtr_percentile <- genes_hm_mtr_percentile[!duplicated(genes_hm_mtr_percentile$mtr_id), ]
        mtr_row_annot <- rowAnnotation(Percentile = genes_hm_mtr_percentile$percentile_rank, 
                                       col = list(Percentile = per_row_col), na_col = "white", 
                                       annotation_label = "%",
                                       annotation_legend_param = list(Percentile = 
                                                                        list(direction = "horizontal", title = "Percentile", title_position = "topcenter")))
        hm_gg_mtr_plot <- Heatmap(data_hm_mtr, name = "Scaled expression", top_annotation = mtr_annot_col, right_annotation = mtr_row_annot,
                                  column_title = "M. truncatula heatmap",
                                  show_row_names = mtr_hm_plot_row_names, show_column_names = TRUE, column_order = c("Nodule 1", "Nodule 2", "Nodule 3", "Root 1", "Root 2", "Root 3"),
                                  show_row_dend = TRUE, show_column_dend = FALSE, cluster_rows = row_dend_mtr, row_split = as.numeric(val_k_mtr),
                                  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
        hm_gg_mtr <- draw(hm_gg_mtr_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
        vals$hm_gg_mtr <- hm_gg_mtr
        out_table_mtr <- unique(mtr_annot[mtr_annot$mtr_id %in% input_ids_split_mtr$mtr_id, ])
        out_table_mtr <- merge(out_table_mtr, mtr_sample_data[ , c(1,9)], by = "mtr_id")
        out_table_mtr <- out_table_mtr[!is.na(out_table_mtr$mtr_id), ]
        out_table_mtr <- unique(out_table_mtr)
        out_table_mtr <- merge(out_table_mtr, clust_IDs_mtr, by = "mtr_id")
        vals$hm_gg_clust_table_mtr <- out_table_mtr
        print(hm_gg_mtr)
      })
      output$hm_down_mtr <- downloadHandler(
        filename = function(){paste(input$hm_title_mtr, '.pdf', sep = '_')},
        content = function(file){
          pdf(file, width = input$hm_mtr_width, height = input$hm_mtr_height)
          print(vals$hm_gg_mtr)
          dev.off()
        })
      output$hm_tbl_mtr <- renderTable({
        out_table_mtr_clust <- vals$hm_gg_clust_table_mtr
        out_table_mtr_clust
      })
      output$download_mtr_clust <- downloadHandler(
        filename = function() {
          paste(input$hm_title_mtr, 'cluster_mtr.csv', sep = '_')
        },
        content = function(file) {
          write.csv(vals$hm_gg_clust_table_mtr, file, row.names = FALSE)
        }
      )
      
}

shinyApp(ui = ui, server = server)

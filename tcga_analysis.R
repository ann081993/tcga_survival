#### TCGA survival analysis
# Author: Seungchan An
library(survival)
library(ggpubr)
library(cowplot)
library(survminer)
library(progress)

# Load TCGA dataset ----
tcga_names <- strsplit("LAML 	Acute Myeloid Leukemia
ACC 	Adrenocortical Cancer
BLCA 	Bladder Urothelial Carcinoma
LGG 	Brain Lower Grade Glioma
BRCA 	Breast Invasive Carcinoma
CESC 	Cervical & Endocervical Cancer
CHOL 	Cholangiocarcinoma
COAD 	Colon Adenocarcinoma
ESCA 	Esophageal Carcinoma
GBM 	Glioblastoma Multiforme
HNSC 	Head & Neck Squamous Cell Carcinoma
KICH 	Kidney Chromophobe
KIRC 	Kidney Clear Cell Carcinoma
KIRP 	Kidney Papillary Cell Carcinoma
LIHC 	Liver Hepatocellular Carcinoma
LUAD 	Lung Adenocarcinoma
LUSC 	Lung Squamous Cell Carcinoma
DLBC 	Diffuse Large B-Cell Lymphoma
MESO 	Mesothelioma
OV 	Ovarian Serous Cystadenocarcinoma
PAAD 	Pancreatic Adenocarcinoma
PCPG 	Pheochromocytoma & Paraganglioma 
PRAD 	Prostate Adenocarcinoma
READ 	Rectum Adenocarcinoma
SARC 	Sarcoma
SKCM 	Skin Cutaneous Melanoma
STAD 	Stomach Adenocarcinoma
TGCT 	Testicular Germ Cell Tumor
THYM 	Thymoma
THCA 	Thyroid Carcinoma
UCS 	Uterine Carcinosarcoma
UCEC 	Uterine Corpus Endometrioid Carcinoma
UVM 	Uveal Melanoma", split = "\n")[[1]]
tcga_types <- trimws(sapply(strsplit(tcga_names, split = "\t"), "[[", 1))
tcga_types_full <- trimws(sapply(strsplit(tcga_names, split = "\t"), "[[", 2))

setwd("/mnt/sdd1/tcga")
tcga_types <- readLines("./tcga_types.txt")
annotation <- read.csv("./gencode.v22.annotation.csv", stringsAsFactors = FALSE)
genes <- annotation$gene

find_gene <- function(g) {
        return(genes[grep(g, genes)])
}

# Expression analysis ----
# function ge_tcga
# get gene expression level for all cancer types using a target gene as an input
ge_tcga <- function(target_gene) {
        ge_data <- NULL
        ind_gene <- which(target_gene == genes)
        pb <- progress_bar$new(format = "Progress: [:bar] :percent, Elapsed time :elapsedfull",
                               total = length(tcga_types), width = 80, clear = FALSE)
        for(t in tcga_types) {
                exp <- read.csv(paste0("./data/TCGA-", t, ".htseq_fpkm.tsv.gz"),
                                stringsAsFactors = FALSE, sep = "\t", skip = ind_gene, nrows = 1, header = FALSE)
                exp <- as.numeric(exp[, -1])
                exp <- data.frame(type = paste0(t, " (n=", length(exp), ")"),
                                  exp = exp, stringsAsFactors = FALSE)
                ge_data <- rbind(ge_data, exp)
                pb$tick()
        }
        ge <- list()
        ge$data <- ge_data
        expplot <- ggstripchart(ge_data, x = "type", y = "exp",
                                color = "darkgray",
                                xlab = "TCGA Cancer Types",
                                ylab = expression(paste(log[2], "(FPKM + 1)",sep=""))) +
                geom_hline(aes(yintercept = median(ge_data$exp)), color = alpha("blue", 0.5)) +
                labs(title = paste0(target_gene, " expression")) +
                scale_y_continuous(expand = c(0, 0)) +
                geom_boxplot(alpha = 0, outlier.color = NA, width = 0.5) +
                stat_summary(fun = median, geom = "point", shape = 18,
                             size = 2.5, color = "red") +
                rotate_x_text(45) +
                theme(plot.title = element_text(hjust = 0.5, size = 12),
                      axis.title = element_text(size = 12),
                      axis.text = element_text(size = 12))
        ge$expplot <- expplot
        return(ge)
}
ge <- ge_tcga("PDZK1IP1")
ge$expplot

# Survival analysis ----
# fuction ge_survival
# get survival fit from target gene, target cancer pair
ge_survival <- function(target_gene, target_cancer, thrsld = 0.2, return_plot = TRUE, legend_dodge = FALSE, regression = FALSE) {
        ind_gene <- which(target_gene == genes)
        samples <- read.csv(paste0("./data/TCGA-", target_cancer, ".htseq_fpkm.tsv.gz"),
                            stringsAsFactors = FALSE, sep = "\t", skip = 0, nrows = 1, header = FALSE, as.is = TRUE)
        samples <- as.character(samples[, -1])
        exp <- read.csv(paste0("./data/TCGA-", target_cancer, ".htseq_fpkm.tsv.gz"),
                        stringsAsFactors = FALSE, sep = "\t", skip = ind_gene, nrows = 1, header = FALSE)
        exp[1, 1] == annotation$Ensembl_ID[ind_gene]
        exp <- as.numeric(exp[, -1])
        
        if(mean(exp) == 0) {
                fit <- NULL
                fit$pval <- NA
                fit$high_exp <- NA
                fit$low_exp <- NA
                fit$exp_fc <- NA
                fit$exp_pval <- NA
                if(regression == TRUE) fit$day <- NA
        } else {
                cutoff <- quantile(exp, probs = c(thrsld, 1 - thrsld))
                group <- rep("", length(exp))
                group[exp >= cutoff[2]] <- paste0("High ", target_gene, " (n=", sum(exp >= cutoff[2]), ")")
                group[exp <= cutoff[1]] <- paste0("Low ", target_gene, " (n=", sum(exp <= cutoff[1]), ")")
                high_exp <- exp[exp >= cutoff[2]]
                low_exp <- exp[exp <= cutoff[1]]
        
                samples <- samples[group != ""]
                group <- group[group != ""]
        
                samples <- data.frame(samples = samples, 
                                      group = group, stringsAsFactors = FALSE)
        
                survival <- read.csv(paste0("./data/TCGA-", target_cancer, ".survival.tsv.gz"),
                                     stringsAsFactors = FALSE, sep = "\t")
                survival <- merge(x = samples, y = survival, by.x = "samples", by.y = "sample")
                survival$OS.time <- survival$OS.time / 30
                fit <- survfit(Surv(OS.time, OS) ~ group, data = survival)
                fit$pval <- surv_pvalue(fit, data = survival)$pval
                fit$high_exp <- high_exp
                fit$low_exp <- low_exp
                fit$exp_fc <- 2^(mean(high_exp) - mean(low_exp))
                fit$exp_pval <- t.test(high_exp, low_exp)$p.val
                names(fit$strata) <- gsub("group=", "", names(fit$strata))
                if(regression == TRUE) {
                        reg <- survreg(Surv(OS.time, OS) ~ group, data = survival)
                        pct <- c(0.5, 0.7)
                        ptime <- predict(reg, newdata=data.frame(group=names(table(survival$group))), type='quantile', p=pct)
                        day <- round(as.numeric(ptime) * 30)
                        names(day) <- c("high 50% surv", "low 50% surv", "high 70% surv", "low 70% surv")
                        fit$day <- day
                }
                if(return_plot == TRUE) {
                        fit$survplot <- ggsurvplot(fit, data = survival, xlab = "Follow-Up Time (Months)", ylab = "Percent Survival",
                                                   palette = c("red", "blue"), 
                                                   risk.table = FALSE, conf.int = TRUE)
                        plab <- paste0("Log-rank P ",
                                       ifelse(fit$pval < 0.001, "< 0.001",
                                              paste0("= ", round(fit$pval, 4))))
                        fit$survplot$plot <- fit$survplot$plot +
                                scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
                                scale_x_continuous(expand = c(0, 0)) +
                                labs(title = paste0(target_gene, " expression ~ overall survival in ", target_cancer),
                                     fill = plab, color = plab) +
                                #scale_color_discrete(labels = c(1,2)) +
                                guides(color = guide_legend(title.position = ifelse(legend_dodge, "top", "bottom"))) +
                                theme(aspect.ratio = 1,
                                      plot.title = element_text(hjust = 0.5, size = 12), 
                                      legend.title.align = 1,
                                      legend.title = element_text(size = 12, hjust = 1),
                                      legend.text = element_text(size = 12), 
                                      legend.position = c(0.98, ifelse(legend_dodge, 0.12, 0.93)),
                                      legend.justification = "right",
                                      legend.background = element_blank())
                        fit$survplot$plot <- ggpar(fit$survplot$plot, 
                                                   font.x = c(12),
                                                   font.y = c(12),
                                                   font.tickslab = c(12))
                }
        }
        return(fit)
}
fit <- ge_survival("ADIPOQ", "BRCA", thrsld = 0.2)
fit$survplot$plot
fit$pval

# fuction ge_surv_screen
# get p-values from survival fit for all cancer types using a target gene as an input
ge_surv_screen <- function(target_gene, thrsld = 0.2, regression = TRUE) {
        screen <- NULL
        pb <- progress_bar$new(format = "Progress: [:bar] :percent, Elapsed time :elapsedfull",
                               total = length(tcga_types), width = 80, clear = FALSE)
        for(t in tcga_types) {
                fit <- ge_survival(target_gene = target_gene, target_cancer = t,
                                   thrsld = thrsld, return_plot = FALSE, regression = TRUE)
                screen <- rbind(screen, c(fit$pval, as.integer(fit$day)))
                pb$tick()
        }
        rownames(screen) <- tcga_types
        colnames(screen) <- c("p value", "high 50% surv", "low 50% surv", "high 70% surv", "low 70% surv")
        return(screen)
}
ge_surv_screen("PDZK1IP1")
ge_surv_screen("LINC00853")

sc02 <- ge_surv_screen("NCOR1", thrsld = 0.2)
sc05 <- ge_surv_screen("NCOR1", thrsld = 0.5)

ggbarplot(sc05, x = "tcga_types", y = "pval") + rotate_x_text(45) + geom_hline(aes(yintercept = 0.05))

# GBM LGG LUAD LUSC MESO SARC KICH KIRC UVM BLCA
fit02 <- ge_survival(target_gene = "ADIPOQ", target_cancer = "THYM", thrsld = 0.2, legend_dodge = F)
fit05 <- ge_survival(target_gene = "ADIPOQ", target_cancer = "THYM", thrsld = 0.5, legend_dodge = F)
plot_grid(fit02$survplot$plot, fit05$survplot$plot, ncol = 2, nrow = 1)

#### Survial analysis - screening mode ----
genes <- c("ADIPOQ", "LEP")

pvals <- NULL
for(n in 1:length(genes)) { #
        surv <- ge_surv_screen(genes[n])
        pvals <- cbind(pvals, surv)
        cat(n, "/", length(genes), genes[n], "\n")
}

#### Download dataset
# download FPKM data
for(t in tcga_types) {
        file.name <- paste0("TCGA-", t, ".htseq_fpkm.tsv.gz")
        download.file(url = paste0("https://gdc.xenahubs.net/download/", file.name),
                      destfile = paste0("./data/", file.name))
}

# download survival data
for(t in tcga_types) {
        file.name <- paste0("TCGA-", t, ".survival.tsv.gz")
        download.file(url = paste0("https://gdc.xenahubs.net/download/", file.name),
                      destfile = paste0("./data/", file.name))
}

#!/usr/bin/env Rscript

# Make an HMM to do DGX using the depmixS4 package

# Load libraries
library(depmixS4)
library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(stringi)
library(dplyr)
library(readr)
library(stringr)

#_____________________________________________________________________________________
# Real data; X chromosome gene expression
#_____________________________________________________________________________________
# Constants
COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3"
METRICS <- "~/XIST/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"

# Create ensemblGenome object for storing Ensembl genomic annotation data
ENS <- ensemblGenome()

# Read in .gff annotation file as ensemblGenome object
read.gtf(ENS, GENCODE)

# Metadata
Metrics <- read_tsv(METRICS) # Contains tissue sample info
Phenotypes <- read_tsv(PHENOTYPES) # Contains sex info

# Read in counts
Gene_Cts <- fread(COUNTS)

# Get annotations for X chromosome
X_Annot <- extractSeqids(ENS, 'X')

# List genes on X chromosome
X_Genes <- X_Annot@ev$gtf[, c("gene_id","gene_name", "start", "end")]

# Remove numbers after decimal in gene ID to get the gene IDs, not the exon IDs
X_Genes$gene_id <- sub("\\.\\d+$", "", X_Genes$gene_id)

# Remove duplicate genes
# Only difference bw transcript and gene IDs is the decimal after the gene ID
X_Genes <- X_Genes[!duplicated(X_Genes$gene_id), ]

# Order genes sequentially by location on chromosome (i.e. by start position)
X_Genes <- X_Genes[order(X_Genes$start),]

# Drop columns in Metrics that aren't needed 
Metrics <- Metrics %>% select(Sample, Note)

# Rename column
colnames(Metrics)[2] <- "Tissue"

# Get list of individual GTEx IDs
Individual_IDs <- unique(str_extract(Metrics$Sample, "GTEX-[0-9A-Z]+"))

# Remove any missing values
which(is.na(Individual_IDs)) # 738; last item in list
Individual_IDs <- Individual_IDs[!is.na(Individual_IDs)]

# For each tissue, make a df of samples from that tissue and store in list.
Types <- unique(Metrics$Tissue)

Tissue_Lst <- list()
for (i in Types){
  Tissue_Lst[[i]] <- Metrics[Metrics$Tissue == i,] 
}

# Add column containing sample ID to each df in list
for (i in seq_along(Tissue_Lst)){
  Tissue_Lst[[i]]$ID <- str_extract(Tissue_Lst[[i]]$Sample, "GTEX-[0-9A-Z]+")
}

# Check for missing values
sapply(Tissue_Lst, function(x) sum(is.na(x))) # Cells - Leukemia cell line (CML) 

# Drop cell lines
Cells <- c("Cells - Leukemia cell line (CML)", "Cells - EBV-transformed lymphocytes", "Cells - Transformed fibroblasts")
Tissue_Lst <- Tissue_Lst[!names(Tissue_Lst) %in% Cells]

# Get list of female IDs
Fem.IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Get list of sample sex for each tissue
Sex_Lst <- lapply(Tissue_Lst, function(x) {
  x <- with(x['ID'], ifelse(ID %in% Fem.IDs, "Female", "Male"))
})

# Add column containing sample sex to each df in list 
Tissue_Lst <- mapply(cbind, Tissue_Lst, "Sex"=Sex_Lst, SIMPLIFY=F)

# For each individual, make a data frame of samples that comes from the same 
# person and store in list.
Ind_Tissues <- list()
for (i in Individual_IDs){
  Ind_Tissues[[i]] <- Metrics[Metrics$Sample %like% i, ]
}

# Get list of sample replicates
# i.e. people who have multiple samples for the same tissue type
Dup_Ind_Tissues <- lapply(Ind_Tissues, function(x) {
  x[duplicated(x[,2]), ]
})
names(Dup_Ind_Tissues) <- names(Ind_Tissues)

Sample_Replicates <- lapply(Dup_Ind_Tissues, function(x) {
  if (as.character(x['Sample']) != "character(0)") {
    as.character(x['Sample'])
  } else {}
})
Sample_Replicates <- unlist(Sample_Replicates)

# Remove sample replicates
Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[!(x$ID %in% Sample_Replicates), ]
})

# Organize sample counts by sex and tissue type; drop replicates
# For each individual, make a data frame of gene counts from samples that come 
# from the same person and store in list.
Gene_Cts <- data.frame(Gene_Cts, stringsAsFactors = F) # was both data.table and data frame
colnames(Gene_Cts) <- str_replace_all(colnames(Gene_Cts), pattern = "\\.", replacement = "-")

# Rename columns
names(Gene_Cts)[1:2] <- c("gene_id", "gene_name")

# Remove decimals in gene IDs in gene counts df
Gene_Cts$gene_id <- sub("\\.\\d+$", "", Gene_Cts$gene_id)

# Split tissue lists by sex
Split_Func <- function(x, y){
  x <- x[x[['Sex']] == y,]
  return(x)
}
f.Tissue_Lst <- Map(Split_Func, x=Tissue_Lst, y='Female')
m.Tissue_Lst <- Map(Split_Func, x=Tissue_Lst, y='Male')

# Drop tissues females/males do not have, respectively
f.remove <- c("Prostate", "Testis", "Cells - Leukemia cell line (CML)")
m.remove <- c("Ovary", "Uterus", "Vagina", "Fallopian Tube", "Cervix - Ectocervix",
              "Cervix - Endocervix", "Cells - Leukemia cell line (CML)")

f.Tissue_Lst <- f.Tissue_Lst[!names(f.Tissue_Lst) %in% f.remove]
m.Tissue_Lst <- m.Tissue_Lst[!names(m.Tissue_Lst) %in% m.remove]

# Sort count data into list of dfs by tissue 
Sort_Func <- function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
}
f.Tissue_Counts <- lapply(f.Tissue_Lst, Sort_Func)
m.Tissue_Counts <- lapply(m.Tissue_Lst, Sort_Func)

# Remove sample replicates
Drop_Replicates <- function(x){
  x <- x[, !(names(x) %in% Sample_Replicates)]
  return(x)
}
f.Tissue_Counts <- lapply(f.Tissue_Counts, Drop_Replicates)
m.Tissue_Counts <- lapply(m.Tissue_Counts, Drop_Replicates)

# Grab heart tissues from 3 females and 3 males
f.Heart <- f.Tissue_Counts[[5]][,1:5] 
m.Heart <- m.Tissue_Counts[[5]][,1:5] 

# Rename cols so I can take the difference between them
colnames(f.Heart)[3:5] <- c("Sample_1", "Sample_2", "Sample_3")
colnames(m.Heart)[3:5] <- c("Sample_1", "Sample_2", "Sample_3")

# Take the difference between the paired samples
Heart_Diff <- f.Heart[,1:2]
Heart_Diff$Delta_1 <- f.Heart$Sample_1 - m.Heart$Sample_1
Heart_Diff$Delta_2 <- f.Heart$Sample_2 - m.Heart$Sample_2
Heart_Diff$Delta_3 <- f.Heart$Sample_3 - m.Heart$Sample_3
head(Heart_Diff)

# 
#_____________________________________________________________________________________
# Simulate data 
#_____________________________________________________________________________________
# Simulate counts using a Gaussian distribution 
Gaussian_Simulate <- function(Num_Samples, Transition, Emission){
    state = 1
    observations<-0
    N <- Num_Samples
    pi_tran <- Transition 
    emission_parm <- matrix(Emission,ncol=2,byrow = T)

    for(i in 2:N){
        uni<-runif(1)
        if(runif(1)<pi_tran[1]){
            state[i]<-1
            observations[i]<-rnorm(n=1,mean = emission_parm[1,1],sd=emission_parm[1,2])
        }else if(uni<pi_tran[2]){
            state[i]<-2
            observations[i]<-rnorm(n=1,mean = emission_parm[2,1],sd=emission_parm[2,2])
        }else{
        state[i]<-3
        observations[i]<-rnorm(n=1,mean = emission_parm[3,1],sd=emission_parm[3,2])
        }
    }

    # Combine observations, states into df
    results <- data.frame(state=rep(NA,N), obs=rep(NA,N))
    results$obs <- observations
    results$state <- state
    results$state <- c("down", "nodiff", "up")[results$state]
    return(cbind(gene=1:N, results))
}


# Simulate scenario
set.seed(2)
# emission: down(mean, sd), baseline(mean, sd), up(mean, sd)
results <- Gaussian_Simulate(Num_Samples=100, Transition=c(0.1,0.3,0.4), Emission=c(-5,5,0,5,15,5))
head(results)

# Observe results
#ggplot(results, aes(x=gene, y=obs)) + geom_line()

#_____________________________________________________________________________________
# Function to make an HMM with depmix with simulated data 
#_____________________________________________________________________________________
# HMM with depmix
# haven't been able to get it to work with the trstart parameter :/
# if I explicitely set the emc parameter to FALSE, says the starting values are not feasible
fit.hmm <- function(observation, dataFrame, inStart, respStart, EM){
    mod <- depmix(observation~1, data=dataFrame, nstates=3,
                  instart=inStart, respstart=respStart)# trstart=trStart
    fit.mod <- fit(mod, emc=em.control(rand=EM))
    est.states <- posterior(fit.mod)
    head(est.states)
    tabl <- table(est.states$state, dataFrame$state)
    dataFrame$est.state.labels <- c(colnames(tabl)[which.max(tabl[1,])],
                                    colnames(tabl)[which.max(tabl[2,])],
                                    colnames(tabl)[which.max(tabl[3,])])[est.states$state]
    est.states$gene <- 1:100
    colnames(est.states)[2:4] <-c(colnames(tabl)[which.max(tabl[1,])],
                                  colnames(tabl)[which.max(tabl[2,])],
                                  colnames(tabl)[which.max(tabl[3,])])
    hmm.post.df <- melt(est.states, measure.vars=c("up", "down", "nodiff"))

    print(table(dataFrame[,c("state", "est.state.labels")]))
    # return list of 3 items; original results with est states appended, posterior info, summary stats
    return(list(hmm=dataFrame, hmm.post.df=hmm.post.df, sumstats=fit.mod))
}
hmm1 <- fit.hmm(dataFrame=results, observation=results$obs, inStart=NULL, respStart=NULL, EM=TRUE)
hmm2 <- fit.hmm(dataFrame=results, observation=results$obs, inStart=runif(3), respStart=NULL, EM=TRUE)
hmm3 <- fit.hmm(dataFrame=results, observation=results$obs, inStart=runif(3), respStart=runif(6), EM=TRUE)


#_____________________________________________________________________________________
# Function to make an HMM with depmix with real data 
#_____________________________________________________________________________________
# Make a seperate matrix for each delta vector
tmp1 <- data.frame(gene=1:length(Heart_Diff$Delta_1), obs=Heart_Diff$Delta_1)

mod <- depmix(obs~1, data=tmp1, nstates=100)
fit.mod <- fit(mod)
est.states <- posterior(fit.mod)
head(est.states)

#_____________________________________________________________________________________
# Plots to show how well the HMM fits the data and estimate the hidden states
#_____________________________________________________________________________________
mycols <- c("darkgreen", "turquoise", "darkblue")
cols <- sapply(results$state,function(x){
  if(x=="up"){mycols[1]}else if(x=="down"){mycols[2]}else{mycols[3]} 
})


plot.hmm.output <- function(model.output){
    g0 <- (ggplot(model.output[[1]], aes(x = gene, y = obs)) + 
           geom_line() +
           theme(axis.ticks = element_blank(), axis.title.y = element_blank())) %>% ggplotGrob
    g1 <- (ggplot(model.output[[1]], aes(x = gene, y = state, fill = state, col = state)) +
           geom_bar(stat = "identity", alpha = I(0.7)) +
           scale_fill_manual(values = mycols, name = "States", 
           labels = c("up", "down", "average")) +
           scale_color_manual(values = mycols, 
                              name = "States", 
                              labels = c("up", "down", "average")) +
           theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
           labs(y = "actual state")) %>% ggplotGrob
    g2 <- (ggplot(model.output[[1]], 
                  aes(x = gene,
                      y = est.state.labels,
                      fill = est.state.labels,
                      col = est.state.labels)) +
        geom_bar(stat = "identity", alpha = I(0.7)) +
        scale_fill_manual(values = mycols, name = "States", 
                          labels = c("up", "down", "average")) +
        scale_color_manual(values = mycols, name = "State", 
                           labels = c("up", "down", "average")) +
        theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
        labs(y = "estimated state")) %>% ggplotGrob
    g3 <- (ggplot(model.output$hmm.post.df, aes(x = gene, y = value, col = variable)) + 
           geom_line() +
           scale_color_manual(values = mycols, name = "State",
                              labels = c("up", "down", "average")) +
        theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
        labs(y = "posterior prob.")) %>%
        ggplotGrob
    g0$widths <- g1$widths
    return(grid.arrange(g0, g1, g2, g3, widths = 1, nrow = 4))
}
plot.hmm.output(hmm1)
plot.hmm.output(hmm2)
plot.hmm.output(hmm3)

#_____________________________________________________________________________________
# Bar plot of the resulting BIC/AIC 
#_____________________________________________________________________________________
plot_list <- c(BIC(hmm1[[3]]), BIC(hmm2[[3]]), BIC(hmm3[[3]]))
barplot(plot_list,
        main="Barplot of model BIC values", 
        xlab="models",
        names.arg=c("mod.1", "mod.2", "mod.3"),
        ylab=c("Bayesian Information Criterion"),
        ylim=c(0, max(plot_list)+100))
        

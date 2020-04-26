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

#_____________________________________________________________________________________
# Real data
#_____________________________________________________________________________________
# Constants
#COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"

# Read in counts
#Gene_Cts <- fread(COUNTS)

# Shape RNA seq data into right format for depmixS4

#_____________________________________________________________________________________
# Different way to simulate data 
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
ggplot(results, aes(x=gene, y=obs)) + geom_line()

#_____________________________________________________________________________________
# Function to make an HMM with depmix 
#_____________________________________________________________________________________
# HMM with depmix
# Changes we want to make to compare: 
fit.hmm <- function(observation, dataFrame, inStart, respStart, EM){
    mod <- depmix(observation~1, data=dataFrame, nstates=3,
                  instart=inStart, respstart=respStart)# trstart=trStart, respstart=repStart)
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
    return(list(dataFrame=dataFrame, hmm.post.df=hmm.post.df))
}
hmm1 <- fit.hmm(dataFrame=results, observation=results$obs, inStart=runif(3), respStart=runif(6), EM=TRUE)
# haven't been able to get it to work with the trstart parameter :/
# if I explicitely set the emc parameter to FALSE, says the starting values are not feasible

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




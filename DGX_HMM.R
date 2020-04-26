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
# Modify to simulate 3 states 
#_____________________________________________________________________________________
Simulate <- function(N, dice.val, count.avg, switch.val){
    up.dice <- sample(1:dice.val, N, replace=T) 
    down.dice <- sample(1:dice.val, N, replace=T)
    nodiff.dice <- sample(1:dice.val, N, replace=T)
    up.count.avg <- rpois(N, count.avg[1])
    down.count.avg <- rpois(N, count.avg[2])
    nodiff.count.avg <- rpois(N, count.avg[3])
    # States
    results <- data.frame(state=rep(NA,N), obs=rep(NA,N), dice=rep(NA,N))

    results$state[1] <- "nodiff"
    results$obs[1] <- nodiff.count.avg[1]
    results$dice[1] <- nodiff.dice[1]

    for (k in 2:N){
        if(results$state[k-1]=="up"){
             # if in the up state and rolled value > switch, stay in up
            if(results$dice[k-1] > switch.val){
                 results$state[k] <- "up"
                 results$obs[k] <- up.count.avg[k]
                 results$dice[k] <- up.dice[k]
             # if < switch, move to down state
             }else if(results$dice[k-1] < switch.val){
                 results$state[k] <- "down"
                 results$obs[k] <- down.count.avg[k]
                 results$dice[k] <- down.dice[k]
             # if = switch, move to no diff state
             }else{
                 results$state[k] <- "nodiff"
                 results$obs[k] <- nodiff.count.avg[k]
                 results$dice[k] <- nodiff.dice[k]
             }
         # if currently in down state
         }else if(results$state[k-1]=="down"){
              if(results$dice[k-1] > switch.val){
                  results$state[k] <- "up"
                  results$obs[k] <- up.count.avg[k]
                  results$dice[k] <- up.dice[k]
             # if < switch, move to down state
             }else if(results$dice[k-1] < switch.val){
                 results$state[k] <- "down"
                 results$obs[k] <- down.count.avg[k]
                 results$dice[k] <- down.dice[k]
             # if = switch, move to no diff state
             }else{
                 results$state[k] <- "nodiff"
                 results$obs[k] <- nodiff.count.avg[k]
                 results$dice[k] <- nodiff.dice[k]
             }
         # if currently in nodiff state
         }else {
             if(results$dice[k-1] > switch.val){
                  results$state[k] <- "up"
                  results$obs[k] <- up.count.avg[k]
                  results$dice[k] <- up.dice[k]
             # if < switch, move to down state
             }else if(results$dice[k-1] < switch.val){
                 results$state[k] <- "down"
                 results$obs[k] <- down.count.avg[k]
                 results$dice[k] <- down.dice[k]
             # if = switch, move to no diff state
             }else{
                 results$state[k] <- "nodiff"
                 results$obs[k] <- nodiff.count.avg[k]
                 results$dice[k] <- nodiff.dice[k]
             }
         }
    }
   return(cbind(gene=1:N, results))
}

# Simulate scenario
set.seed(2)
results <- Simulate(N=100, count.avg=c(100,1,50), dice.val=6, switch.val=3)
head(results)

# Observe results
ggplot(results, aes(x=gene, y=obs)) + geom_line()


#_____________________________________________________________________________________
# Different way to simulate data 
#_____________________________________________________________________________________
state = 1
observations<-0
N <- 100
pi_tran <- c(0.1, 0.3, 0,4)
emission_parm <- matrix(c(1,5,10,5,15,5),ncol=2,byrow = T)

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


#_____________________________________________________________________________________
# Function to make an HMM with depmix 
#_____________________________________________________________________________________
# HMM with depmix
# Changes we want to make to compare: 
fit.hmm <- function(observation, dataFrame, inStart, trStart, repStart){
    mod <- depmix(observation~1, data=dataFrame, nstates=3,
                  instart=inStart, trstart=trStart, respstart=repStart)
    fit.mod <- fit(mod)
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
hmm1 <- fit.hmm(dataFrame=results, observation=results$obs, 
                inStart=c(0.333,0.333,0.333), trStart=c(0.333,0.333,0.333),
                repStart=c(0.1,0.1,0.1,0.1,0.1,0.5))
# keep getting different errors; if I set respstart to length 3 it says it should be 6 and vice versa

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


##_____________________________________________________________________________________
## Other kinds of mixture models 
##_____________________________________________________________________________________
## Make a dataframe of randomly generated observations and label which pair they come from
## Generate observations using a poisson distribution; lambda paramater sets the mean.
#Data <- data.frame(Observation=rpois(5, lambda=10), Pair=c(1:5))
#Data
#
## Make dependent mixture model with example data; fit model; print summary
#set.seed(1)
#Model.1 <- depmix(Observation ~ 1, data = Data, nstates = 3, family=poisson())
#Fm.1 <- fit(Model.1, emc=em.control(rand=FALSE))
#Fm.1
#summary(Fm.1)
#
## Predict the states by estimating the posterior
#estStates.1 <- posterior(Fm.1)
#head(estStates.1)
#
## What happens if I randomly generate the start values with EM?
#Fm.2 <- fit(Model.1, emc=em.control(rand=TRUE))
#Fm.2
#summary(Fm.2)
#
## Predict the states by estimating the posterior
#estStates.2 <- posterior(Fm.2)
#head(estStates.2)
#
## What happens if I do not explicitely set the family parameter to poisson?
#set.seed(1)
#Model.2 <- depmix(Observation ~ 1, data = Data, nstates = 3)
#Fm.3 <- fit(Model.2, emc=em.control(rand=FALSE))
#Fm.3
#summary(Fm.3)
#
# Predict the states by estimating the posterior
#estStates.3 <- posterior(Fm.3)
#head(estStates.3)


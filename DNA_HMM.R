#!/usr/bin/env Rscript

# Make the transition matrix
states <- c("AT-rich", "GC-rich") 
AT_stateProbs <- c(0.7, 0.3)             
GC_stateProbs <- c(0.1, 0.9)             
transition_Matrix <- matrix(c(AT_stateProbs, GC_stateProbs), 2, 2, byrow = TRUE) 
rownames(transition_Matrix) <- states
colnames(transition_Matrix) <- states
transition_Matrix 

# Make the emission matrix
nucleotides <- c("A", "C", "G", "T")   
AT_nucProbs <- c(0.39, 0.1, 0.1, 0.41) 
GC_nucProbs <- c(0.1, 0.41, 0.39, 0.1) 
emission_Matrix <- matrix(c(AT_nucProbs, GC_nucProbs), 2, 4, byrow = TRUE) 
rownames(emission_Matrix) <- states
colnames(emission_Matrix) <- nucleotides
emission_Matrix

#---------------------------------------------------------------------------------------------
# Function to generate a DNA eventSequence of length seqLength 
#---------------------------------------------------------------------------------------------
makeHMMseq <- function(transitionMatrix, emissionMatrix, initialProbs, seqLength)
{
   nucleotides <- c("A", "C", "G", "T")   
   states <- c("AT-rich", "GC-rich")
   resultSequence <- character()             
   resultStates <- character()            
                                            
   firstState <- sample(states, 1, rep=TRUE, prob=initialProbs)
   probabilities <- emissionMatrix[firstState,]
   firstNucleotide <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
   resultSequence[1] <- firstNucleotide         
   resultStates[1] <- firstState             

   for (i in 2:seqLength)
   {
      prevState <- resultStates[i-1]         
      stateProbs <- transitionMatrix[prevState,]
      state  <- sample(states, 1, rep=TRUE, prob=stateProbs)
      probabilities <- emissionMatrix[state,]
      nucleotide <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
      resultSequence[i] <- nucleotide             
      resultStates[i] <- state                  
   }

   for (i in 1:length(resultSequence))
   {
      nucleotide <- resultSequence[i]
      state <- resultStates[i]
      print(paste("Position", i, ", State", state, ", Nucleotide = ", nucleotide))
   }
   return(resultSequence)
}
# Set initial probs to 50-50, since either state is assumed equally likely
initial_Probs <- c(0.5, 0.5)
DNA_Seq <- makeHMMseq(transition_Matrix, emission_Matrix, initial_Probs, 10)

#---------------------------------------------------------------------------------------------
# Function to make the Viterbi matrix
#---------------------------------------------------------------------------------------------
viterbiMatrix <- function(eventSequence, transitionMatrix, emissionMatrix)
{
   eventSequence <- toupper(eventSequence)
   numStates <- dim(transitionMatrix)[1]
   V <- matrix(NA, nrow = length(eventSequence), ncol = dim(transitionMatrix)[1])
   V[1, ] <- 0
   V[1,1] <- 1
   for (i in 2:length(eventSequence)) 
   {
      for (l in 1:numStates) 
      {
         stateNucleotideProbi <- emissionMatrix[l,eventSequence[i]]
         V[i,l] <- stateNucleotideProbi * max(V[(i-1),] * transitionMatrix[,l])
      }
  }
  return(V)
}
viterbiMatrix(DNA_Seq, transition_Matrix, emission_Matrix)

#---------------------------------------------------------------------------------------------
# Function to estimate the most likley states (AC- or GC- rich) given a DNA sequence
#---------------------------------------------------------------------------------------------
viterbiFunc <- function(eventSequence, transitionMatrix, emissionMatrix)
  {
     states <- rownames(emissionMatrix)

     V <- viterbiMatrix(eventSequence, transitionMatrix, emissionMatrix)
     mostProbablePath <- apply(V, 1, function(x) which.max(x))

     prevNucleotide <- eventSequence[1]
     prevMostProbableState <- mostProbablePath[1]
     prevMostProbableStateName <- states[prevMostProbableState]
     startIndex <- 1
     for (i in 2:length(eventSequence))
     {
        nucleotide <- eventSequence[i]
        mostProbableState <- mostProbablePath[i]
        mostProbableStateName <- states[mostProbableState]
        if (mostProbableStateName != prevMostProbableStateName)
        {
           print(paste("Positions",startIndex,"-",(i-1), "Most probable state = ", prevMostProbableStateName))
           startIndex <- i
        }
        prevNucleotide <- nucleotide
        prevMostProbableStateName <- mostProbableStateName
     }
     print(paste("Positions",startIndex,"-",i, "Most probable state = ", prevMostProbableStateName))
   }
viterbiFunc(DNA_Seq, transition_Matrix, emission_Matrix)


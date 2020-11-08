library(igraph)
library(extraDistr)

HollywoodModel <- function(NoPaths, Nv, v, alpha, theta){
  # returns a list of interaction network graphs generated from the 
  # Hollywood(alpha, theta, {v}) model.
  # **TO DO: Remove possibility of loops** 
  
  PathLengths <- rcat(NoPaths, v) # number of edges in each network
  sequence <- list() # sequence of edge sequences
  D <- c(1) # tracks how often each vertex has appeared, recorded as D[i] for vertex i
  V <- 1 # tracks how many unique vertices have appeared
  
  sequence[[1]] <- c(1, rep(0, PathLengths[1]-1)) # length of path one
  for(j in 2:(PathLengths[1])){
    t <- rcat(1, c((D - alpha), theta + alpha*V)) #new vertex
    sequence[[1]][j] <- t 
    if(t %in% 1:V){
      D[t] <- D[t] + 1
    }
    else{
      D <- c(D, 1)
      V <- V + 1
    }
  }
  
  for(i in 2:NoPaths){
    sequence[[i]] <- rep(0, PathLengths[i])
    for(j in 1:PathLengths[i]){
      t <- rcat(1, c((D - alpha), theta + alpha*V))
      sequence[[i]][j] <- t
      if(t %in% 1:V){
        D[t] <- D[t] + 1
      }
      else{
        D <- c(D, 1)
        V <- V + 1
      }
    }
  }
  
  result <- list()
  for(i in 1:NoPaths){
    A <- matrix(0, V, V)
    for(j in 1:(length(sequence[[i]]) - 1)){
      A[sequence[[i]][j], sequence[[i]][j+1] ] <- 1
      A[sequence[[i]][j+1], sequence[[i]][j] ] <- 1
    }
    result[[i]] <- graph_from_adjacency_matrix(A, mode="undirected")
  }
  
  return(result)
}
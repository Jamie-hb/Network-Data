d_Hamming <- function(g, h){
  # returns the Hamming distance between the graph g and h, defined as the number of disparate edges
  return(0.5 * sum( abs( get.adjacency(g) - get.adjacency(h) ) ) )
}

MatrixExp <- function(A){
  # returns the exponential of a matrix
  result <- matrix(0, dim(A)[1], dim(A)[1])
  for(i in 0:100){
    result <- result + ((A^i)/factorial(i))
  }
  return(result)
}

DiffusionDistance <- function(g, h){
  # returns the diffusion distance between g and h, with parameter t
  t <- 1
  Ag <- get.adjacency(g)
  Ah <- get.adjacency(h)
  Nv <- dim(Ag)[1]
  Lg <- -Ag
  Lh <- -Ah
  for(i in 1:Nv){
    Lg[i,i] <- sum(Ag[i,])
    Lh[i,i] <- sum(Ah[i,])
  }
  
  F <- MatrixExp(-t*Lg) - MatrixExp(-t*Lh)
  return(sum(F^2))
}

dSNM <- function(g, Gm, gamma, d, phi){
  # the pmf the Spherical Network Model with centroid Gm, parameter gamma, distance
  # metric g, and increasing function phi
  return( exp( -gamma * phi(d(g, Gm)) ) )
}

phi <- function(x){
  # defines phi
  return(x)
}

CentredErdosReyni <- function(Gm, p){
  # returns a graph formed by flipping each edge of Gm independently with probability p
  A <- get.adjacency(Gm)
  Nv <- dim(A)[1]
  Anew <- matrix(0, Nv, Nv)
  
  for(i in 1:(Nv-1)){
    t <- rbinom(Nv-i, rep(1, Nv-i), A[i,(i+1):Nv]*(1-p) + (1-A[i,(i+1):Nv])*p)
    Anew[i,(i+1):Nv] <- t
    Anew[(i+1):Nv,i] <- t
  }
  
  return(graph_from_adjacency_matrix(Anew, mode="undirected"))
}

CentredErdosReyniEasy <- function(Gm, p){
  A <- get.adjacency(Gm)
  Nv <- dim(A)[1]
  Anew <- matrix(0, Nv, Nv)
  
  for(i in 1:(Nv-1)){
    for(j in i:Nv){
      if(A[i,j] == 1){
        t <- rbinom(1, 1, 1-p)
        Anew[i,j] <- t
        Anew[j,i] <- t
      }
      else{
        t <- rbinom(1, 1, p)
        Anew[i,j] <- t
        Anew[j,i] <- t
      }
    }
  }
  
  return(graph_from_adjacency_matrix(Anew, mode="undirected"))
}

CEMAux <- function(x, p){
  if(x==1){
    return(p)
  }
  else{
    return(1-p)
  }
}

dCEM <- function(g, Gm, p){
  Ag <- get.adjacency(g)
  Am <- get.adjacency(Gm)
  A <- abs(Ag - Am)
  
  CEMAux <- function(x){
    if(x==1){
      return(p)
    }
    else{
      return(1-p)
    }
  }
  
  A <- sapply(A, CEMAux)
  
  for(i in 1:dim(A)[1]){
    A[i,i] <- 1
  }
  
  return(prod(A))
}

SNMMetropolisHastings <- function(NoSamples, Gm, gamma, distance, phi, p){
  # Uses a Metropolis-Hastings algorithm to sample from the Spherical Network Family model
  # with centroid Gm, parameter gamma, distance metric distance, increasing function phi.
  # Proposal Distribution: Centred Erdos-Reyni model with flip parameter p
  sample <- list()
  sample[[1]] <- Gm
  for(i in 1:NoSamples){
    CandidateGraph <- CentredErdosReyni(sample[[i]], p)
    AcceptanceProb <- dSNM(CandidateGraph, Gm, gamma, distance, phi) / dSNM(sample[[i]], Gm, gamma, distance, phi)
    # print(AcceptanceProb)
    
    t <- runif(1, 0, 1)
    if(t <= AcceptanceProb){
      sample[[i+1]] <- CandidateGraph
    }
    else{
      sample[[i+1]] <- sample[[i]]
    }
  }
  return(sample)
}
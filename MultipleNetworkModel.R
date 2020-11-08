StochasticBlockGraphWithClasses <- function(Nv, classprobs, P){
  # Returns a graph from the stochastic block model, with length(classprobs) classes, Nv edges, 
  # P(vertex in class i) = classprobs[i], P(edge between vertex of class i and vertex
  # of class j) = P[i,j]
  
  classes <- rcat(Nv, classprobs) #vector of vertex classes
  A <- matrix(0, nrow=Nv, ncol=Nv) # incidence matrix
  
  for(i in 1:(Nv - 1)){
    # for row i of A, generates a vector of 0s and 1s drawn from appropriate Bernoulli distributions, 
    # and assigns to it A[i,i+1:Nv] and A[i+1:Nv,i]
    t <- rbinom(Nv - i, rep(1, Nv - i), P[(classes[(i+1):Nv] - 1)*dim(P)[1] + classes[rep(i, Nv - i)]])
    A[i, (i+1):Nv] <- t
    A[(i+1):Nv,i] <- t
  }
  
  l <- list()
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  V(g)$class <- classes
  l[[1]] <- g
  l[[2]] <- classes
  return(l)
}

MultipleNetworkSample <- function(Nv, classprobs, B, NoSamples, P, Q){
  # Draws a network A from a stochastic block model with arguments Nv (number
  # of vertices), classprobs (probabilities of a node being in each class), and B
  # (inter class edge probabilities). Then returns NoSamples noisy realisations 
  # of A, with P[i,j] and Q[i,j] providing the false positive and false negative 
  # probabilities respectively for an edge between nodes from classes i and j
  #
  # return = list of NoSamples+1 graphs, with graphs 1:NoSamples being the noisy realisations,
  # and graph NoSamples+1 being the true graph
  
  Reality <- StochasticBlockGraphWithClasses(Nv, classprobs, B) # the actual graph of which we observe noisy realisations
  A <- get.adjacency(Reality[[1]])
  classes <- Reality[[2]]
  sample <- list()
  M <- matrix(0, nrow=Nv, ncol=Nv) # adjacency matrix of observed network
  
  for(n in 1:NoSamples){
    for(i in 1:(Nv-1)){
      t <- rbinom(Nv - i, 1, (1 - A[i,(i+1):Nv]) * P[ (classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]
                  + A[i,(i+1):Nv] * (1 - Q[(classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]))
      M[i, (i+1):Nv] <- t
      M[(i+1):Nv, i] <- t
    }
    sample[[n]] <- graph_from_adjacency_matrix(M, mode="undirected")
  }
  
  sample[[NoSamples + 1]] <- Reality[[1]]
  
  return(sample)
}

MultipleNetworkSampleBM <- function(Nv, classes, B, NoSamples, P, Q){
  # Draws a network A from a block model with arguments Nv (number
  # of vertices), classes (class membership of each node), and B
  # (inter class edge probabilities). Then returns NoSamples noisy realisations 
  # of A, with P[i,j] and Q[i,j] providing the false positive and false negative 
  # probabilities respectively for an edge between nodes from classes i and j
  
  Reality <- BlockGraph(classes, B) # the actual graph of which we observe noisy realisations
  A <- get.adjacency(Reality)
  sample <- list()
  M <- matrix(0, nrow=Nv, ncol=Nv) # adjacency matrix of observed network
  
  for(n in 1:NoSamples){
    for(i in 1:(Nv-1)){
      t <- rbinom(Nv - i, 1, (1 - A[i,(i+1):Nv]) * P[ (classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]
                  + A[i,(i+1):Nv] * (1 - Q[(classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]))
      M[i, (i+1):Nv] <- t
      M[(i+1):Nv, i] <- t
    }
    sample[[n]] <- graph_from_adjacency_matrix(M, mode="undirected")
  }
  
  sample[[NoSamples + 1]] <- Reality 
  
  return(sample)
}

MultipleNetworkFromGraph<- function(g, classes, NoSamples, P, Q){
  # Allows one to fix the "real" graph and vertex class memberships in the measurement error, rather than drawing it from 
  # a single network model. 
  A <- get.adjacency(g)
  sample <- list()
  Nv <- length(degree(g))
  M <- matrix(0, nrow=Nv, ncol=Nv) # adjacency matrix of observed network
  
  for(n in 1:NoSamples){
    for(i in 1:(Nv-1)){
      t <- rbinom(Nv - i, 1, (1 - A[i,(i+1):Nv]) * P[ (classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]
                  + A[i,(i+1):Nv] * (1 - Q[(classes[(i+1):Nv] - 1) * dim(P)[1] + classes[rep(i, Nv - i)] ]))
      M[i, (i+1):Nv] <- t
      M[(i+1):Nv, i] <- t
    }
    sample[[n]] <- graph_from_adjacency_matrix(M, mode="undirected")
  }
  
  sample[[NoSamples + 1]] <- g
  
  return(sample)
}  



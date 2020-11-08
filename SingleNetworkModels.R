library(igraph) # package for network data
library(extraDistr) # package containing commands to draw from lots of distributions
library(MASS) # contains the command mvrnorm, allowing one to draw from a multivariate normal distributions

BernouilliGraph <- function(Nv, p){ 
  # Returns a graph with Nv vertices, and each possible
  # edge present independently with probability p
  #start <- Sys.time()
  A <- matrix(0, nrow=Nv, ncol=Nv) # adjacency matrix
  for(i in 1:(Nv - 1)){ 
    t <- rbinom(Nv - i, 1, p)
    A[i,(i+1):Nv] <- t
    A[(i+1):Nv,i] <- t
  }
  
  #end <- Sys.time()
  #print(end - start)
  return(graph_from_adjacency_matrix(A, mode="undirected"))
}

BlockGraph <- function(classes, P){
  # Returns a graph from the block model with length(classes) vertices, and probability of an edge
  # between two nodes from classes i and j as P[i,j]
  # classes is a vector of the class membership of each vertex (so vertex i is in class
  # classes[i])
  
  Nv <- length(classes) # number of vertices
  A <- matrix(0, nrow=Nv, ncol=Nv) # incidence matrix
  
  for(i in 1:(Nv - 1)){
    # for row i of A, generates a vector of 0s and 1s drawn from appropriate Bernouilli distributions, 
    # and assigns to it A[i,i+1:Nv] and A[i+1:Nv,i]
    t <- rbinom(Nv - i, rep(1, Nv - i), P[(classes[(i+1):Nv] - 1)*dim(P)[1] + classes[rep(i, Nv - i)]])
    A[i, (i+1):Nv] <- t
    A[(i+1):Nv,i] <- t
  }
  
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  V(g)$class <- classes
  return(g)
}

StochasticBlockGraph <- function(Nv, classprobs, P){
  # Returns a graph from the stochastic block model, with length(classprobs) classes, Nv edges, 
  # P(vertex in class i) = classprobs[i], P(edge between vertex of class i and vertex
  # of class j) = P[i,j]
  
  classes <- rcat(Nv, classprobs) #vector of vertex classes
  A <- matrix(0, nrow=Nv, ncol=Nv) # incidence matrix
  
  for(i in 1:(Nv - 1)){
    # for row i of A, generates a vector of 0s and 1s drawn from appropriate Bernouilli distributions, 
    # and assigns to it A[i,i+1:Nv] and A[i+1:Nv,i]
    t <- rbinom(Nv - i, rep(1, Nv - i), P[(classes[(i+1):Nv] - 1)*dim(P)[1] + classes[rep(i, Nv - i)]])
    A[i, (i+1):Nv] <- t
    A[(i+1):Nv,i] <- t
  }
  
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  V(g)$class <- classes
  return(g)
}

MixedMembershipSBM <- function(Nv, K, alpha, B){
  # ***TO CHECK: make sure code is doing what I want it to***
  # Returns a graph sampled from the Mixed Membership Stochastic Block Model with
  # - Nv nodes
  # - K classes
  # - class membership weights alpha (1 x K vector)
  # - interclass interaction matrix B (K x K matrix)
  
  A <- matrix(0, nrow=Nv, ncol=Nv)
  
  pi <- rdirichlet(Nv, alpha)
  
  for(i in 1:(Nv - 1)){
    senderMemberships <- rmnom(Nv - i, 1, pi[i,])
    receiverMemberships <- rmnom(Nv - i, rep(1, Nv -i), pi[(i+1):Nv,])
    t <- rbinom(Nv - i, 1, as.vector(t(senderMemberships %*% B %*% t(receiverMemberships))))
    A[i, (i+1):Nv] <- t
    A[(i+1):Nv,i] <- t
  }
  
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  return(g)
  
}

distance <- function(x, y){
  # computes ||x - y||_2 for two vectors x, y in R^n 
  n <- length(x)
  t <- rep(0, n)
  for(i in 1:n){
    t[i] <- (x[i] - y[i])^2
  }
  return(sqrt(sum(t)))
}


explogit <- function(x){
  return(exp(x) / (1 + exp(x)))
}

LatentPositionClusterModel <- function(Nv, means, variances, classprobs, alpha){
  # Returns a graph sampled from the Latent Position Distance Cluster Model, with
  # - constants classprobs representing the prob a a node belonging to each of G clusters
  # - latent position of node i in cluster j sampled from MVN_d(means[[j]], variances[j]I)
  # - logodds(edge from i to j) = alpha - d(i,j)
  # NOTE: means must be a matrix with rows as co-ordinates
  
  A <- matrix(0, nrow=Nv, ncol=Nv) #adjacency matrix
  d <- length(variances)
  
  Z <- matrix(0, nrow=Nv, ncol=d) #matrix of latent positions (row i = pos of node i)
  classes <- rcat(Nv, classprobs, 1:length(variances)) # cluster memberships of each node,
  #  determining which MVN distribution to draw latent position from
  for(i in 1:Nv){
    Z[i,] <- mvrnorm(1, means[[classes[i]]], variances[classes[i]] * diag(d))
  }
  
  for(i in 1:(Nv - 1)){
    distances <- apply(as.matrix(Z[(i+1):Nv,]), 1, distance, y=Z[i,])
    logodds <- alpha - distances
    p <- explogit(logodds)
    t <- rbinom(Nv - i, 1, p)
    A[i,(i+1):Nv] <- t
    A[(i+1):Nv] <- t
  }
  
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  V(g)$class <- classes
  return(g)
  
}
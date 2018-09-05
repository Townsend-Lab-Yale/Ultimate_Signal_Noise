TreeCollapse.signal.noise<-function(freqs,subratevector,ratevector,treePath){
  require(ape) #ape is good for tree structures
  
  if (!file.exists(treePath)) {
    print("Input Tree does not exist")
    return(NULL)
  }
  
  if(sum(freqs)!=1){
    print("Base Frequencies sum to 1")
    return(NULL)
  }
  if(length(subratevector)!=6){
    print("subratevector must have length 6")
    print("Error in values a-e")
    return(NULL)
  }
  
  
  
  #unload command args
  a <- subratevector[1]
  b <- subratevector[2]
  c <- subratevector[3]
  d <- subratevector[4]
  e <- subratevector[5]
  f <- subratevector[6]
  piT <- freqs[1]
  piC <- freqs[2]
  piA <- freqs[3]
  piG <- freqs[4]
  
  frequ<-c(piT,piC,piA,piG)
  
  ##Shared internal functions
  calcMu <- function(a, b, c, d, e, f, piT, piA, piG, piC) {
    Mu_ <- (1 / 2) / (a * piT * piC + b * piT * piA + c * piT * piG +
                        d * piC * piA + e * piC * piG + f * piA * piG)
    return(Mu_)
  }
  #calc Q substitution rate matrix
  calcQ <- function(a, b, c, d, e, f, piT, piA, piG, piC) {
    #Dim 4X4
    Q <- matrix(nrow = 4, ncol = 4)
    Q[1, 1] <- ((-a) * piC) - (b * piA) - (c * piG)
    Q[1, 2] <- a * piC
    Q[1, 3] <- b * piA
    Q[1, 4] <- c * piG
    Q[2, 1] <- a * piT
    Q[2, 2] <- ((-a) * piT) - (d * piA) - (e * piG)
    Q[2, 3] <- d * piA
    Q[2, 4] <- e * piG
    Q[3, 1] <- b * piT
    Q[3, 2] <- d * piC
    Q[3, 3] <- ((-b) * piT) - (d * piC) - (f * piG)
    Q[3, 4] <- f * piG
    Q[4, 1] <- c * piT
    Q[4, 2] <- e * piC
    Q[4, 3] <- f * piA
    Q[4, 4] <- ((-c) * piT) - (e * piC) - (f * piA)
    return(Q)
  }
  #calcTrueQ is an internal function, just multiplies Mu_*Q. Returns final Q matrix
  calcTrueQ <- function(Q, Mu_) {
    return(Mu_ * Q)
  }
  #function to get the eigenvalues, vectors, and inverse of the vectors
  doEigenMaths <- function(Q) {
    #Obtain the eigenvalues and vectors
    evects <- eigen(Q)
    evalues <- evects$values
    evectors <- evects$vectors
    
    #Reorder in ascending order, swap rows.
    evalues <- evalues[c(4, 3, 2, 1)] #same as mathematica
    evectors <- evectors[c(4, 3, 2, 1), c(4, 3, 2, 1)]
    tev <- evectors
    #Get inverse
    itev <- solve(tev)
    return(list(evalues, evectors, itev))
  }
#### END SHARED INTERNAL FUNCTIONS
  
####Start Collapse Internal Functions
  
  #gets probability transition matrix
  P <- function(lamda, T_, eigenStuff) {
    evalues <- eigenStuff[[1]]
    tev <- eigenStuff[[2]]
    itev <- eigenStuff[[3]]
    return(tev %*% (diag(exp(
      evalues * lamda * T_
    )) %*% itev))
  }
  #get the joint probaility for the first conditional probability
  JointProb3 <-
    function(freqs,
             root,
             char1,
             char2,
             T1,
             T2,
             lamda,
             eigenStuff) {
      return(freqs[root] * P(lamda, T1, eigenStuff)[root, char1] * P(lamda, T2, eigenStuff)[root, char2])
    }
  #get the conditional probability for the first conditional entropy
  CondProb3 <-
    function(freqs,
             root,
             char1,
             char2,
             T1,
             T2,
             lamda,
             eigenStuff) {
      output <- JointProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff)
      outhold <- 0
      for (i in 1:length(freqs)) {
        outhold <-
          outhold + JointProb3(freqs, i, char1, char2, T1, T2, lamda, eigenStuff)
      }
      return(output / outhold)
    }
  #Calculate the conditional entropy for 2 branches at average rate lamda
  CondEntropy3 <- function(freqs, T1, T2, lamda, eigenStuff) {
    output <- 0
    len_freq <- 1:length(freqs)
    for (root in len_freq) {
      for (char1 in len_freq) {
        for (char2 in len_freq) {
          JP3 <- JointProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff)
          logCP3 <-
            log(CondProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff))
          output <- output + (JP3 * logCP3)
        }
      }
    }
    return(-1 * output)
  }
  
  #Compute the joint probability for the second conditional probability
  JointProb2 <-
    function(freqs,
             rootprime,
             charprime,
             Tprime,
             lamda,
             eigenStuff) {
      return(freqs[rootprime] * P(lamda, Tprime, eigenStuff)[rootprime, charprime])
    }
  #compute the conditional prob for the second conditional entropy
  CondProb2 <-
    function(freqs,
             rootprime,
             charprime,
             Tprime,
             lamda,
             eigenStuff) {
      hold1 <- JointProb2(freqs, rootprime, charprime, Tprime, lamda, eigenStuff)
      tmp <- 0
      for (i in 1:length(freqs)) {
        tmp <- tmp + JointProb2(freqs, i, charprime, Tprime, lamda, eigenStuff)
      }
      return(hold1 / tmp)
    }
  #get the second part of the conditional entropy
  CondEntropy2 <- function(freqs, Tprime, lamda, eigenStuff) {
    output <- 0
    len_freq <- 1:length(freqs)
    for (root in len_freq) {
      for (charprime in len_freq) {
        cprob <- log(CondProb2(freqs, root, charprime, Tprime, lamda, eigenStuff))
        output <-
          output + JointProb2(freqs, root, charprime, Tprime, lamda, eigenStuff) *
          cprob
      }
    }
    return(-1 * output)
  }
  #get mutual information/entropy for collapse by finding root
  InfoEquivalent <- function(freqs, T1, T2, lamda, eigenStuff) {
    cond_entropy3 <- CondEntropy3(freqs, T1, T2, lamda, eigenStuff)
    #internal function to minimize (pick Tprime such that
    # cond_entropy2-cond_entropy3 is 0)
    to_find_root <- function(Tprime) {
      cond_entropy2 <- CondEntropy2(freqs, Tprime, lamda, eigenStuff)
      return(cond_entropy2 - cond_entropy3)
    }
    #output<-uniroot(to_find_root,min(c(T1,T2)),maxiter = 500)
    #output<-uniroot(Vectorize(to_find_root),c(0.000000001,min(T1,T2)),maxiter = 10000)$root
    output <-
      uniroot((to_find_root), c(0.000000001, min(T1, T2)), maxiter = 10000)$root
    
    return(output)
  }
  
  #there isn't a good way to collapse the tree branches natively in ape
  #this function helps get around that!
  bind.tip <- function(tree,
                       tip.label,
                       edge.length = NULL,
                       where = NULL) {
    if (is.null(where))
      where <- length(tree$tip) + 1
    tip <- list(
      edge = matrix(c(2, 1), 1, 2),
      tip.label = tip.label,
      edge.length = edge.length,
      Nnode = 1
    )
    class(tip) <- "phylo"
    obj <- bind.tree(tree, tip, where = where)
    return(obj)
  }
  #using the current tree, collapse the outmost nested pair of branches
  getInfoEquiv <- function(freqs, tree, lamda, eigenStuff) {
    #for now, use average+0.1+nodeabove
    #this needs to be the infocollapse
    this <-
      InfoEquivalent(freqs,
                     tree$edge.length[1],
                     tree$edge.length[2],
                     lamda,
                     eigenStuff)
    #this<-mean(tree$edge.length[1]+tree$edge.length[2])+0.1
    nodeabove = tree$edge[1, 1] #do checks
    findCon = which(tree$edge[, 2] == nodeabove)
    addLen <- tree$edge.length[findCon]
    this <- this + addLen
    return(c(this, tree$edge[findCon, 1]))
  }
  #//////END INTERNAL FUNCTION DEFINITIONS////////
  
  ####Get Shared EigenInformation
  
  eigenStuff <-
    doEigenMaths(calcTrueQ(
      calcQ(a, b, c, d, e, f, piT, piA, piG, piC),
      calcMu(a, b, c, d, e, f, piT, piA, piG, piC)
    ))
  tree<-read.tree(treePath)
  #Internal function to evaluate lamda
  evalLambda <- function(lamda,eigenStuff,tree) {
    evalues <- eigenStuff[[1]]
    tev <- eigenStuff[[2]]
    itev <- eigenStuff[[3]]
    x <- tree
    #if the tree is rooted, unroot it.
    if (is.rooted(x)) {
      x <- unroot(x)
    }
    #sort our tree edges to get innermost pair of branches
    x = reorder(x, "postorder")
    iter = 0 #to keep track of writing out files (unique names)
    while (length(x$edge.length) > 5) {
      #until we have a quartet
      x = reorder(x, "postorder") #sort
      #get branch 1 and 2
      first <- x$edge[1, 2]
      second = x$edge[2, 2]
      #return branch location info and values for collapse
      pair <- getInfoEquiv(freqs, x, lamda, eigenStuff)
      #using the value in pair, collapse the tree one level
      x2 <-
        bind.tip(
          x,
          tip.label = paste0(first, "_", second),
          edge.length = pair[1],
          where = pair[2]
        )
      x3a <- drop.tip(x2, tip = second)
      x3b <- drop.tip(x3a, tip = first)
      x <- x3b
      iter <- iter + 1
    }
    internodeIndex <-
      which(x$edge[, 2] == 6) #definition of internode (node 5 to node 6)
    quartet_branches <- x$edge.length[-internodeIndex]
    internode <- c(quartet_branches, x$edge.length[internodeIndex])
    p <- list()
    p <- array(, dim = c(5, 4, 4))
    for (v in 1:length(internode)) {
      p[v, , ] <- (tev %*% (diag(exp(
        evalues * lamda * internode[v]
      )) %*% itev))
    }
    correct <- 0
    wrong1 <- 0
    wrong2 <- 0
    
    for (original_character in 1:4) {
      for (internode_character in 1:4) {
        for (leaf_character_1 in 1:4) {
          for (leaf_character_2 in 1:4) {
            if (leaf_character_1 != leaf_character_2) {
              correct <- correct + (frequ[original_character] *
                                      p[5, original_character, internode_character] *
                                      p[1, original_character, leaf_character_1] *
                                      p[2, original_character, leaf_character_1] *
                                      p[3, internode_character, leaf_character_2] *
                                      p[4, internode_character, leaf_character_2])
              wrong1 <-
                wrong1 + (frequ[original_character] *
                            p[5, original_character, internode_character] *
                            p[1, original_character, leaf_character_1] *
                            p[2, original_character, leaf_character_2] *
                            p[3, internode_character, leaf_character_1] *
                            p[4, internode_character, leaf_character_2])
              wrong2 <-
                wrong2 + (frequ[original_character] *
                            p[5, original_character, internode_character] *
                            p[1, original_character, leaf_character_1] *
                            p[2, original_character, leaf_character_2] *
                            p[3, internode_character, leaf_character_2] *
                            p[4, internode_character, leaf_character_1])
            }
          }
        }
      }
    }
    all <- c(correct, wrong1, wrong2)
    return(all)
  }
  #Initialize blanks
  eYsum <- 0
  eX1sum <- 0
  eX2sum <- 0
  eY2sum <- 0
  eX12sum <- 0
  eX22sum <- 0
  eX1Ysum <- 0
  eX2Ysum <- 0
  eX1X2sum <- 0
  
  for (lmbda in ratevector) {
    all <- evalLambda(lmbda,eigenStuff,tree)
    y <- all[1]
    x1 <- all[2]
    x2 <- all[3]
    eYsum <- eYsum + y
    eX1sum <- eX1sum + x1
    eX2sum <- eX2sum + x2
    
    eY2sum <- eY2sum + (y ^ 2)
    eX12sum <- eX12sum + (x1 ^ 2)
    eX22sum <- eX22sum + (x2 ^ 2)
    
    eX1Ysum <- eX1Ysum + (x1 * y)
    eX2Ysum <- eX2Ysum + (x2 * y)
    eX1X2sum <- eX1X2sum + (x1 * x2)
  }
  
  Mu_1 <- eYsum - eX1sum
  Mu_2 <- eYsum - eX2sum
  
  
  Sigma_1 <-
    sqrt(eX1sum + eYsum - eX12sum - eY2sum + 2 * eX1Ysum)
  Sigma_2 <-
    sqrt(eX2sum + eYsum - eX22sum - eY2sum + 2 * eX2Ysum)
  Rho_ <-
    (-eX1X2sum + eX1Ysum + eX2Ysum + eYsum - eY2sum) / (Sigma_1 * Sigma_2)
  
  #Internal function for integration
  FofT <- function(t) {
    F1ofT = ((1 / Sigma_1) * dnorm((t - Mu_1) / Sigma_1) * pnorm(Rho_ * (t - Mu_1) /
                                                                   (Sigma_1 * sqrt(1 - Rho_ * Rho_)) - (t - Mu_2) / (Sigma_2 * sqrt(1 - Rho_ *
                                                                                                                                      Rho_))))
    F2ofT = ((1 / Sigma_2) * dnorm((t - Mu_2) / Sigma_2) *
               pnorm(Rho_ * (t - Mu_2) / (Sigma_2 * sqrt(1 - Rho_ * Rho_)) - (t - Mu_1) /
                       (Sigma_1 * sqrt(1 - Rho_ * Rho_))))
    return(F1ofT + F2ofT)
  }
  
  princtree <- integrate(FofT,-Inf,-.5)
  prpolytomy = integrate(FofT,-.5, .5)
  prcortree  = integrate(FofT, .5, Inf)
  
  print(paste0("Probablility Correct: ", prcortree$value))
  print(paste0("Probability Incorrect: ", princtree$value))
  print(paste0("Probability Polytomy: ", prpolytomy$value))
  return(c(prcortree$value,princtree$value,prpolytomy$value))
}

#Independent test call
#TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=c(rep(0.05,500)),treePath="test.tre")

#Dependent tet calls
#TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=new,treePath="test.tre")
# set.seed(66)
# for(i in 1:5){
#   print("i")
#   print(i)
#   write.tree(rtree(n=5*i,rooted = FALSE),file = "rtree.tre")
#   for(j in 1:5){
#     new<-runif(n=50*j, min=0, max=.7)
#     print("j")
#     print(j)
#     print(sum(TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=new,treePath="rtree.tre")))
#     
#   }
# }

library(igraph)
library(RSQLite)
library(jsonlite)
library(R6)

VotingSimulation <- 
  R6Class ("VotingSimulation",
           public = list(
             #setup
             configuration="Small World,Degree=4,Rewire=0.01",
             DISTANCENEIGHBORS=2,  # average degree =  2*distance for small world graphs
             PROBREWIRE=0.01,  #0 for BA graph
             graphType="SW_min2", 
             N_VECTOR = c(1000),   #size of graph
             ITERATIONS=100, 
             epsilon = 0.001,     
             adjMatrix = NULL,
             invDiagDegMatrix = NULL,
             graph=NULL,
             degreeVector = NULL,
             MAXSTEP=100,
             seeds = c( 0.225, 0.25, 0.275, 0.3, 0.325),
             tbs = c(0.6, 0.8, 1.0),
             SEEDPERCENT_VECTOR=NULL,
             TBPERCENT_VECTOR =NULL,
             PERCENTDECAY_VECTOR = c(0.0,0.1,0.15,0.2),
             K=3,  #num candidates
             
             #simulation state
             n=1,
             step = 1,
             votingComplete = FALSE,
             percentDecay = 0,
             winner = 0,
             currentFreq = NULL,
             currentPrefVector = NULL,
             prevPrefVector = NULL,
             
             #simulation data structures
             thresholdMatrix = NULL,
             percentNeighborsMatrix = NULL,
             PreferenceVectorT1 = NULL,
             PreferenceVector = NULL,
             PreferenceMatrixT1 = NULL,
             PreferenceMatrix = NULL,
             seedIndexVector = NULL,
             seedCountVector = NULL,
             
             #database
             db = NULL,
             
             
             initialize = function() {
               self$db <- dbConnect(SQLite(), dbname="voting.sqlite")
             },
               
             
             newGraph = function() {
               self$graph<- sample_smallworld(dim=1, size=self$n, nei=self$DISTANCENEIGHBORS, p=self$PROBREWIRE)
               #self$graph<- sample_pa( n=self$n, m=3,directed=FALSE)
               
               #ensure a minimum degree 2
               degreeVector <- degree(self$graph)
               for (i in 1:self$n)  {
                 if (degreeVector[i]==0) {
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                 }
                 else if (degreeVector[i]==1) {
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                 }
                
               }
               self$degreeVector = degree(self$graph)  #recompute
               self$adjMatrix= t(as.matrix(as_adjacency_matrix(self$graph)))  #transpose
               self$invDiagDegMatrix=diag(sapply(self$degreeVector, function(i) (1/i)),self$n,self$n)
             },
             
             
             newPreferenceVectorT1 = function() {
               self$PreferenceVectorT1 = rep((self$K+1), self$n)  #initialize to all undecided value K+1
               offset=0
               for (candidate in 1:self$K) {
                 seedsThisCandidate <- self$seedCountVector[candidate]
                 for (i in 1:seedsThisCandidate) {
                   #seedIndexVector is random ordering of vertices
                   r <- self$seedIndexVector[offset+i]  #get next voter row using random ordering of vertices
                   self$PreferenceVectorT1[r]=candidate    
                 }
                 offset <- offset + seedsThisCandidate
               }
               
             },
             
             newPreferenceMatrixT1 = function() {
               self$PreferenceMatrixT1 = matrix( 0, nrow = self$n, ncol = self$K)  #all 0, no Preference
               for (r in (1:self$n)) {
                 if (self$PreferenceVectorT1[r]<= self$K) {   #skip if undecided k+1
                   self$PreferenceMatrixT1[r,self$PreferenceVectorT1[r]] = 1   #seed
                 }
               }
              
             },
             
             
             newThresholdMatrix = function() {
               self$thresholdMatrix = matrix( -1, nrow = self$n, ncol = self$K )  #init all to -1
               offset=0
               for (candidate in 1:self$K) {
                 nSeeds = self$seedCountVector[candidate]
                 nTrueBelievers = ceiling(nSeeds * self$TBPERCENT_VECTOR[candidate])
                 if (nSeeds > 0)  {     #ugh, loop iterates once when  nSeeds 0!
                   for (i in 1:nSeeds) {
                     r = self$seedIndexVector[offset+i]
                     if (i<=nTrueBelievers) {
                       #true believer, threshold 0 for seed candidate,  threshold above 100% for other candidates
                       for (c in 1:self$K) {
                         self$thresholdMatrix[r,c] = ifelse(c==candidate,0,1.1)
                       }
                     }
                     else {
                       #adherent,  seed candidate threshold=min random/degree, other candidates unique random/degree
                       #UPDATE - pick k random numbers between 2..degree.  possible duplicate values if degree-1<k
                       if (self$degreeVector[r] < self$K) {
                         randomDegreeVector = rep(2,self$K)
                       }
                       else {
                         randomDegreeVector = sample(2:self$degreeVector[r], self$K, replace= (self$degreeVector[r]-1<self$K ))
                       }
                       
                       
                       index = which.min(randomDegreeVector)
                       #seed candidate should get min threshold, so swap with min
                       tmp=randomDegreeVector[index]
                       randomDegreeVector[index] = randomDegreeVector[candidate]
                       randomDegreeVector[candidate] = tmp
                       
                       for (c in 1:self$K) {
                         self$thresholdMatrix[r,c] = 1.0*randomDegreeVector[c]/self$degreeVector[r]  # avoid integer division
                       }
                     }
                   }
                 }
                 offset = offset + nSeeds
               }
               #for non-seed voters,  threshold = random(2..degree)/degree
              for (r in 1:self$n) {
                 if (self$degreeVector[r] < self$K) {
                   randomDegreeVector = rep(2,self$K)
                 }
                 else {
                   randomDegreeVector = sample(2:self$degreeVector[r], self$K, replace= TRUE )
                 }
                 if (self$thresholdMatrix[r,1] == -1)  {
                   for (c in 1:self$K) {
                     self$ thresholdMatrix[r,c]  = 1.0*randomDegreeVector[c]/self$degreeVector[r]
                   }
                 }
               }
               
             },
             
             
             PreferenceMatrixTn = function() {
               
               self$prevPrefVector = self$currentPrefVector   #HISTORY
                
               self$PreferenceMatrix[]=0    #reset values of existing matrix 
               self$PreferenceVector[]=(self$K+1)  #init to all undecided
               
               #assign based on neighbors Preference and threshold
               for (r in (1:self$n))   {
                 min<-1.1  #initialize to max threshold more than 100%
                 for (c in (1:self$K))  {
                   
                   #add epsilon in case rounding issues for proportions
                   if (self$percentNeighborsMatrix[r,c]+self$epsilon>=self$thresholdMatrix[r,c])  {
                   #threshold met, test for new min threshold
                     if (self$thresholdMatrix[r,c]<min)  {
                       min<-self$thresholdMatrix[r,c]
                     }
                   }
                 }
                   
                 if (min<1.1)  {
                   #could be several candidates, select random
                   indices <- c()
                   for (c in (1:self$K)) {
                     if (self$thresholdMatrix[r,c] == min & self$percentNeighborsMatrix[r,c]+self$epsilon >= min) {
                       indices <- c(indices,c)
                     }
                   }
                   #sample - if only one item in vector it generates random number from 1..item so need to test length
                   if (length(indices) <= 1) {
                     Preference<-indices[1]
                   } else {
                     Preference<-sample(indices,1)
                   }
                   
                   self$PreferenceMatrix[r,Preference] = 1
                   self$PreferenceVector[r]=Preference  
                 }
               }
              
               self$currentPrefVector = self$PreferenceVector    #HISTORY
               
             },
             
             
             freqTable = function(cVector) {
               kplus1=self$K + 1
               n=self$n
               tmp=table(factor(cVector,levels=c(1:kplus1)))
               #round to 3 digits, proportion of n
               sapply(tmp,function(i) round(1.0*i/n, 3))
             },
             
             #check distance between current and previous preferences
             checkResult = function() { 
               self$currentFreq<- self$freqTable(self$currentPrefVector)
               prevFreq<-self$freqTable(self$prevPrefVector)
               distance<-euc.dist(prevFreq,self$currentFreq)  
               self$votingComplete = (distance<0.01) 
               if (self$votingComplete) {
                 self$winner=as.vector(which.max(self$currentFreq))
               }
               
             },
             
             #weight edges for decay
             randomWeightVector = function ()  {
               weightVector <- rep(1, self$n)  #all voters have default weight 1
               m=ceiling(self$percentDecay * self$n)
               randomSample <- sample(1:self$n, m, replace=FALSE )  #select m voters
               for (i in 1:m) {
                 r <- randomSample[i]  #voter row
                 weightVector[r]<-runif(1, 0.5, 0.9)    #set weight to random value between 0.5 and 0.9
               }
               weightVector
             },
             
             updatePercentNeighborsMatrix = function() {
               if (self$percentDecay > 0.0)  {
                 weightVector<- self$randomWeightVector()
                 weightDiagMatrix<- diag(weightVector) 
                 weightedPreferenceMatrix<- weightDiagMatrix  %*% self$PreferenceMatrix
                 self$percentNeighborsMatrix = self$invDiagDegMatrix %*% self$adjMatrix %*% weightedPreferenceMatrix
               }
               else {
                 self$percentNeighborsMatrix = self$invDiagDegMatrix %*% self$adjMatrix %*% self$PreferenceMatrix
               }
             },
             
             
             saveResult = function()  {
               result=ifelse(self$step ==self$MAXSTEP,"maxsteps","converge")
               ranking<-order(- self$currentFreq)
               winnerPercent<- self$currentFreq[ranking[1]]
               winnerLead <- winnerPercent - self$currentFreq[ranking[2]]
               seedJSON <-toJSON(self$SEEDPERCENT_VECTOR)
               tbJSON <- toJSON(self$TBPERCENT_VECTOR)
               degreeDistributionJSON <- toJSON(degree_distribution(self$graph))
               insertStmt <- sprintf("INSERT INTO session (configuration,iterations,k,type,DISTANCENEIGHBORS,probRewire,n,step,result,finalPercentages,seed,tb,decay,winner,winnerPercent,winnerLead,meanDistance,clusterCoeff,meanDegree, degreeDistribution) 
                                     VALUES ('%s',%d,%d,'%s',%d,%f,%d,%d,'%s','%s','%s','%s',%f,%d,%f,%f,%f,%f, %f, '%s')",
                                     self$configuration, self$ITERATIONS, self$K,self$graphType,self$DISTANCENEIGHBORS,self$PROBREWIRE, self$n,self$step,result,toString(self$currentFreq,sep=','), seedJSON, tbJSON, self$percentDecay,self$winner, winnerPercent, winnerLead,
                                     mean_distance(self$graph),transitivity(self$graph), mean(degree(self$graph)),degreeDistributionJSON)
               dbSendQuery(conn = self$db,insertStmt)
             },
             
             run = function()  {
               for (n in self$N_VECTOR) {
                 self$n=n   #used to be 100, 500, 1000.  Now just 1000
                 for (i in 1:self$ITERATIONS)  {
                    self$newGraph()
                    for (s in self$seeds)  {
                      self$SEEDPERCENT_VECTOR <- c(s, 0.1, 0.1)
                      
                      #dependent on n and seed % vector.  compute outside tb loop so all tb use same seed assignment
                      self$seedCountVector=sapply(self$SEEDPERCENT_VECTOR, function(i) (i*self$n))
                      
                      self$seedIndexVector = sample(1:self$n, sum(self$seedCountVector), replace=FALSE )
                      # compute so same vector used for each tb and decay
                      self$newPreferenceVectorT1()  
                      self$newPreferenceMatrixT1()
                      
                      for (t in self$tbs)  {
                        self$TBPERCENT_VECTOR <- c(t, 1.0, 1.0)
                        #cat("seed " , self$SEEDPERCENT_VECTOR, " tb ",self$TBPERCENT_VECTOR, n,i,"\n")
                    
                        self$newThresholdMatrix()  #based on seed and tb
                        for (percentDecay in self$PERCENTDECAY_VECTOR)  {
                          #cat("decay",percentDecay,"\n")
                          
                          self$percentDecay=percentDecay
                          self$step=1
                          self$votingComplete=FALSE
                          self$PreferenceVector = self$PreferenceVectorT1 #reset so each decay starts with same initial Preferences
                          self$PreferenceMatrix = self$PreferenceMatrixT1 #reset 
                          
                          self$currentPrefVector = self$PreferenceVectorT1   #HISTORY
                          self$prevPrefVector = NULL                         #HISTORY
                          
                          self$updatePercentNeighborsMatrix()
                     
                          #compute next iteration of simulation
                          while (!self$votingComplete && self$step<self$MAXSTEP) {
                            self$step=self$step+1
                            self$PreferenceMatrixTn()
                            self$updatePercentNeighborsMatrix()
                            self$checkResult()
                          }  
                          self$saveResult()  #store result in database
                        }
                      }
                 
                    }
                  }
                 }
              
             }
           )     
           
  )

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


main <- function() {
  mySimulation <- VotingSimulation$new()
  mySimulation$run()
}

#run the simulation
main()

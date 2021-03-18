# install.packages("pracma", repos="http://R-Forge.R-project.org")
# install.packages("ramify")

rm(list = ls())
library(dplyr)
library(tidyr)
library(purrr)
library(Ryacas)
library(stringr)
library(tictoc)
library(pracma)
library(ramify)
# source("create_sym_dict.R")


# ---- symbolic dictionary functions  ----
# calc_household_prob_symbolic as in "create_sym_dict.R" 
# create_symbolic_dictionary as in "create_sym_dict.R" with permutations adjustment

calc_household_prob_symbolic <- function(dict,idx,N_names,I_names) {
  
  
  ag_num <- length(N_names)
  
  stopifnot(length(N_names)==ag_num)
  stopifnot(length(I_names)==ag_num)
  stopifnot(idx<=nrow(dict))
  
  
  Q <- paste0("Q",1:ag_num) %>% ysym()
  Qtext <- paste0(",{",  paste0("Q",1:ag_num,collapse = ","),"})")
  
  
  if(is.na(dict[idx,]$prob)) {
    
    N <- as.numeric(dict[idx,N_names])
    I <- as.numeric(dict[idx,I_names])
    # print(N)
    # print(I)
    n <- sum(N)
    k <- sum(I)
    if(k<n) {
      ni <- N-I
      bin_fact <- prod(choose(N,I))
      prob <-  bin_fact*prod( Q^( ni*(k+1) )) 
      id <- paste0(c(I,I),collapse=',')
      idx2 <- which(dict$id2==id)
      if(length(idx2) > 0) {
        prob <- prob*dict[idx2,]$prob %>%
          yac_str()
        prob_str <- paste0("Expand(",prob,Qtext)
        prob <- ysym(prob_str)
      }
    } else {
      id <- paste0(N,collapse=',')
      idx2 <- which(dict$id1==id)
      idx2 <- setdiff(idx2,idx)
      prob <- ysym(dict[idx2[1],]$prob)
      if(length(idx2)>1) {
        for(i in 2:length(idx2)) {
          prob <- prob + ysym(dict[idx2[i],]$prob)
          prob_str <- paste0("Expand(",prob,Qtext)
          prob <- ysym(prob_str)
        }
      }
      prob <- 1-prob 
    }
    dict[idx,]$prob <- yac_str(prob)
  }
  return (dict)
}


create_symbolic_dictionary <- function(max_members,past_dict = NULL) {
  ysym("MaxEvalDepth(100000)")
  
  # Creating an empty dictionary with all the option for number of 
  # people and number of infected
  ag_names <- 1:length(max_members)
  names(max_members) <- ag_names
  N_ag <- lapply(max_members,function(m) 0:m)
  I_ag <- N_ag
  N_names <- paste0('N.',ag_names)
  I_names <- paste0('I.',ag_names)
  names(N_ag) <- N_names
  names(I_ag) <- I_names
  
  dict <- expand.grid(c(N_ag,I_ag))
  dict <- dict[-1,]
  dict <- dict[-which(apply(dict[,N_names]-dict[,I_names],1,function(r) any(r<0))),]
  rownames(dict) <- 1:nrow(dict)
  dict$id1 <- apply(dict[,N_names],1,function(r) paste0(r,collapse=','))
  dict$id2 <- apply(dict[,c(N_names,I_names)],1,function(r) paste0(r,collapse=','))
  
  # if there is a past dictionary add it, if not add an empty column for prob
  if(!is.null(past_dict)){
    dict <- left_join(dict, select(past_dict,prob,idx),by = c("id2"="idx"))
  } else {
    dict$prob <- NA
  }
  
  saving <- 0
  for(i in 1:nrow(dict)) {
    if(is.na(dict[i,]$prob)){
      dict <- calc_household_prob_symbolic(dict,i,N_names,I_names)
      perm_prob <- prob_permutations(dict[i,]$prob, dict[i,]$id2)
      # perm_prob = ['index', 'permutation', 'id2', 'prob']
      for(j in 1:nrow(perm_prob)){
        perm_index <- match(perm_prob[j,]$id2, dict$id2, nomatch = NA)
        if(!is.na(perm_index) & is.na(dict[perm_index,]$prob)){
          dict[perm_index,]$prob <- perm_prob[j,]$prob
          saving <- saving + 1
        }
      }
      saveRDS(dict,"tempDict1.RDS")
    }
  }
  disp('permutations saved',saving, 'calls for calc_household_prob_symbolic')
  return (as_tibble(dict))
}


# ---- functions ----

make_function_robust <- function(body){
  args <- paste(param_names, collapse = ', ')
  eval(parse(text = paste('f <- function(', args, ') { return(' , body , ')}', sep='')))
}


permutations <- function(n, original_call=TRUE){
  # create all possible permutations over [n] as list of lists
  if( (floor(n)!=n) | (n<1) ) stop('n must be a natural number')
  perms_list <- vector(mode="list", length=factorial(n))
  if (n==1) {
    perms_list[1] <- list(1)
    return(perms_list)
  }
  prev_list <- permutations(n-1, original_call = FALSE)
  ind <- 1
  for(perm in prev_list){
    perms_list[[ind]] <- append(n, perm)
    ind <- ind+1
    for(loc in 2:n){
      perms_list[[ind]] <- append(perm, n, after=loc-1)
      ind <- ind+1
    }
  }
  if (isFALSE(original_call)){return(perms_list)}
  perms_df <- data.frame(index=1:factorial(n))
  perms_df$permutation <- apply(perms_df['index'], 1, function(r) as.list(perms_list[[r]]))
  return(perms_df)
}


id2_permutations <- function(id2){
  # return all distinct permutations except identity
  # id2 <- "2,1,0,0,0,0"
  id_2 <- strsplit(id2, split=',')[[1]]
  perm_dict <- permutations_dict[-nrow(permutations_dict),c('index', 'permutation')]
  perm_dict$id2 <- apply(perm_dict['permutation'], 1, function(r) paste0(id_2[c(unlist(r),unlist(r)+risk_groups)], collapse=','))
  perm_dict <- perm_dict[!duplicated.data.frame(perm_dict['id2']),]
  row.names(perm_dict) <- 1:nrow(perm_dict)
  return (perm_dict)
}


prob_permutations <- function(sym_prob, id2){
  # sym_prob <- "Q1*Q2"
  # id2 <- "1,1,0,0,0,0"
  perms <- id2_permutations(id2)
  outer <- param_names
  names(outer) <- if (isFALSE(B)) paste0('Q.', 1:risk_groups) else as.character(c(paste0('Q.', 1:risk_groups), paste0('B.', 1:risk_groups)))
  inner <- names(outer)
  perms$prob <- NA
  for (j in 1:nrow(perms)){
    names(inner) <- if(isFALSE(B)) param_names[c(unlist(perms[j,"permutation"]))] else param_names[c(unlist(perms[j,"permutation"]))]
    perms[j,"prob"] <- str_replace_all(str_replace_all(sym_prob,inner),outer)
  }
  return (perms)
}


create_imputation_possibilities <- function(id3_template){
  # Create a table for a family where not all house members were tested.
  # Each row in the table is a possible imputation.
  # id2_imputation is the id2 for the sub-family that did not tested
  # id2_imputed is the id2 for the entire family given the imputation
  
  # n=number of family members for each risk group
  # c=number of family members from each group that were tested
  # d=number of positive tests among the c that tested
  n <- id3_template[1:risk_groups]
  d <- id3_template[(risk_groups+1) : (2*risk_groups)]
  c <- id3_template[(2*risk_groups+1) : (3*risk_groups)]
  
  # sub family that tested (index case excluded)
  # testedId2 <- paste(c(id3_template[(2*risk_groups+1) : (3*risk_groups)],id3_template[(risk_groups+1) : (2*risk_groups)]), collapse = ',')
  # testedProb <- dict$prob[dict$id2 == testedId2]
  
  # sub family that did not tested
  maxImpute <- n-c
  
  # create a table with possible imputations
  # the k-th possible imputation - imputationTable[[k]] is
  # ( imputation  , numerator , denominator, imputed family Id2 , imputed family prob )
  # where numerator/denominator = imputation probability
  
  
  ag_names <- 1:risk_groups
  names(maxImpute) <- ag_names
  I_ag <- lapply(maxImpute,function(m) 0:m)
  I_names <- paste0('I.',ag_names)
  names(I_ag) <- I_names
  
  imputationTable <- expand.grid(c(I_ag))
  imputationTable$id2_imputation <- apply(imputationTable[,I_names],1,function(r) paste0(c(maxImpute, r),collapse=','))
  imputationTable$id2_imputed <- apply(imputationTable[,I_names],1,function(r) paste0(c(n, d+r),collapse=','))
  imputationTable$observation_index <- NA
  
  # adding probability for each imputation - normalizing the conditional distribution
  imputationTable$imputation_prob_numerator <- apply(imputationTable["id2_imputed"],1,function(r) dict$prob[dict$id2 == r])
  normalization_factor_sym <- paste(imputationTable["imputation_prob_numerator"][[1]], collapse = '+')
  imputationTable$imputation_prob <- apply(imputationTable["imputation_prob_numerator"],1,function(r) paste('(', r, ')/(', normalization_factor_sym, ')'))
  
  # create the numeric functions and evaluate at point Q
  # imputationTable$imputation_prob_f <- apply(imputationTable["imputation_prob"],1,function(r) make_function_robust(r, risk_groups))
  # imputationTable$imputation_prob_f[[1]](0.4,0.4,0.4)
  
  return (imputationTable[c("observation_index", "id2_imputation", "id2_imputed", "imputation_prob")])
}


create_imputations_table <- function(observationsTable){
  # observationsTable <- observations
  
  # any observation should be a list of size 3*r
  # if the observations presented as strings (obs="3,2,2,1,1,1,2,1,1") we need to call
  # create_imputation_possibilities(strsplit(obs, split=',')[[1]])
  
  # Create one table whose sub-tables are create_imputation_possibilities(observation)
  observationsNum <- nrow(observationsTable)
  imputationsTable <- create_imputation_possibilities(observationsTable[1,])
  imputationsTable[is.na(imputationsTable)] <- 1
  for (i in 2:observationsNum){
    imputationsTable <- rbind(imputationsTable, create_imputation_possibilities(observationsTable[i,]))
    imputationsTable[is.na(imputationsTable)] <- i
  }
  return(imputationsTable)
}


expectation_phase <- function(idsTable, Q, B=NULL){
  # Calculate E[idsTable|Q,B] where idsTable represent a distribution (not normalized, the sum is n).
  # col(idsTable) = (original observation index, id2_imputation, id2_imputed, symbolic probability for that imputation)
  # dict$prob[dict$id2 == id2_imputation] is the probability of the value id2_imputed.
  # the returned value is a data frame ["imputation_prob", "id2_imputed"] where each line represent a possible imputation
  #    for some observation with it probability given Q,B ()
  distributionTable <- data.frame(idsTable["imputation_f_prob"], idsTable["id2_imputed"])
  distributionTable$imputation_prob <- apply(distributionTable["imputation_f_prob"],1,function(r) do.call(what = r[[1]], args = as.list(Q)))
  return (distributionTable[c("imputation_prob", "id2_imputed")])
}


maximization_phase <- function(idsProbs, precision=0.001, barsNum=30, lower=NULL, upper=NULL, B=FALSE){
  # lower <- l_bound
  # upper <- u_bound
  tableToValue <- function(...){
    params <- c(...)
    params_list <- vector(mode = "list", length = length(param_names))
    for (i in seq_along(param_names)){
      params_list[[i]] <- params[i]
    }
    # evaluate log_likelihood for each possible family
    evaluated_functions <- apply(imputationsTable['imputed_f_prob'], 1, function(r) log(do.call(what=r[[1]], args=params_list)))
    # return weighted sum
    return(idsProbs$imputation_prob %*% evaluated_functions)
  }
  # creating rangeTable - grid to maximize over
  ## param_names <- if (isFALSE(B)) paste0('Q.', 1:risk_groups) else as.character(c(paste0('Q.', 1:risk_groups), paste0('B.', 1:risk_groups)))
  valid_lower <- !(is.null(lower) | (isTRUE(B) & length(lower) != 2*risk_groups) | (isFALSE(B) & length(lower) != risk_groups))
  valid_upper <- !(is.null(upper) | (isTRUE(B) & length(upper) != 2*risk_groups) | (isFALSE(B) & length(upper) != risk_groups))
  barsNames <- vector(mode = "character", length = length(param_names))
  for(i in 1:length(param_names)){
    low <- if (valid_lower) lower[i] else 0.001
    high <- if (valid_upper) upper[i] else 0.999
    # reduce the precision of the solution to avoid huge grid
    # actual_bars <- min(((high-low)/precision), barsNum)
    line_i <- linspace(low, high, min(((high-low)/precision), barsNum))
    assign(paste('X', i, sep=''), line_i, pos=sys.frame())
    barsNames[i] <- paste('X', i, sep='')
  }
  bars <- sapply(barsNames, get, simplify = FALSE, USE.NAMES = TRUE)
  names(bars) <- param_names
  rangeTable <- expand.grid(c(bars))
  # up till here we get a table whose rows are the points to maximize over.
  rangeTable$value <- apply(rangeTable[param_names],1,function(r) tableToValue(r))
  maximizerIndex <- argmax(as.matrix(rangeTable$value), rows = FALSE)
  return(rangeTable[maximizerIndex,])
}


# ---- main -----

# -- preparing for E-M

# loading partial observations
observations <- readRDS('observations.RDS')

# extract number of risk groups and family size bound according to the data
risk_groups <- length(observations[1,])/3
max_members <- unname(sapply(observations[1:risk_groups],max))

permutations_dict <- permutations(risk_groups)

# define model's parameters
B=FALSE
param_names <- if (isFALSE(B)) paste0('Q', 1:risk_groups) else as.character(c(paste0('Q', 1:risk_groups), paste0('B', 1:risk_groups)))

# create symbolic dictionary (according to the chosen model)
# max_members <- c(3,3,3) took 7 minutes and save 780 calls to calc_household_prob_symbylic
tic()
dict <- create_symbolic_dictionary(max_members)
toc()
# dict <- readRDS('dict.RDS')

# create all possible imputations, for each family, in terms of id2
imputationsTable <- create_imputations_table(observations)

# add the numeric probability for each imputation
# those functions don't appear in the family dictionary, they have been normalized over each (original) family
imputationsTable$imputation_f_prob <- apply(imputationsTable["imputation_prob"],1,function(r) (make_function_robust(r)))

# add the probabilities for each imputed family from dict (both symbolic and numeric)
imputationsTable$imputed_prob <- dict[match(imputationsTable$id2_imputed,dict$id2),]$prob
imputationsTable$imputed_f_prob <- apply(imputationsTable["imputed_prob"],1,function(r) (make_function_robust(r)))

# imputationsTable=[observation_index, id2_imputation, id2_imputed, imputation_prob(sym+num), imputed_prob(sym+num)]
# ready for E-M! imputationsTable holds all the data needed. 


# -- running E-M

# first call for E-M
# tuning parameters:
#    initial point to calculate expectation
#    bounds define the domain to maximize over
#    precision define the grin within bounds
#    bars number is an upper bound for the grid complexity
#    tolerance and max iteration are the stopping criteria
init_point <- c(0.3, 0.5, 0.7)
l_bound <- c(0.2, 0.4, 0.4)
u_bound <- c(0.6, 0.9, 0.7)
prec <- 0.001
barsNum <- 10
tol <- 0.0001
maxIter <- 10

# maximization over expectation
maxim <- maximization_phase(expectation_phase(imputationsTable, init_point), prec, barsNum, l_bound, u_bound)
maximizer <- maxim[param_names]
maximal_value <- maxim$value

iter <- 0
prev_maximal_value <- maximal_value + 2*tol

# continue
while((iter<10) & (abs(maximal_value-prev_maximal_value)>tol)){
  iter <- iter+1
  prev_maximal_value <- maximal_value
  prev_maximizer <- maximizer
  # l_bound <- c(0.2, 0.4, 0.4)
  # u_bound <- c(0.6, 0.9, 0.7)
  # prec <- 0.001
  # barsNum <- 6
  maxim <- maximization_phase(expectation_phase(imputationsTable, prev_maximizer), prec, barsNum, l_bound, u_bound)
  maximizer <- maxim[param_names]
  maximal_value <- maxim$value
}
disp(iter)
print(maximizer)
disp(maximal_value)




# ---- "working" main ----

# -- load the data - partial observations and extract number of risk groups
# observations <- data.frame(c(2,1),c(1,1),c(3,3),c(1,1),c(0,0),c(1,2),c(1,1),c(0,0),c(2,2))
# names(observations) <- c("n1", "n2", "n3","d1", "d2", "d3", "c1", "c2", "c3" )
# saveRDS(observations, 'observations.RDS')
observations <- readRDS('observations.RDS')
risk_groups <- length(observations[1,])/3


# -- define family size bound according to the data
# max_members <- c(3,3,3)
max_members <- unname(sapply(observations[1:risk_groups],max))


# -- create symbolic dictionary (according to the chosen model)
# tic()
# dict <- create_symbolic_dictionary(max_members)
# toc()
# saveRDS(dict, 'dict.RDS')
dict <- readRDS('dict.RDS')

# add probabilities functions from the symbolic probabilities (using "make_function" with proper variables)
### not nessesery at this point
# dict$f_prob <-apply(dict[,"prob"],1,function(r) make_function_robust(r))

N_names <- names(dict)[1:(risk_groups)]
I_names <- names(dict)[(1+risk_groups):(2*risk_groups)]


# -- create all possible imputations, for each family, in terms of id2
imputationsTable <- create_imputations_table(observations)


# -- add the numeric probability for each imputation
# those functions don't appear in the family dictionary, they have been normalized over each (original) family
imputationsTable$imputation_f_prob <- apply(imputationsTable["imputation_prob"],1,function(r) (make_function_robust(r)))


# -- add the probabilities for each imputed family from dict (both symbolic and numeric)
imputationsTable$imputed_prob <- dict[match(imputationsTable$id2_imputed,dict$id2),]$prob
imputationsTable$imputed_f_prob <- apply(imputationsTable["imputed_prob"],1,function(r) (make_function_robust(r)))

# -- imputationsTable=[observation_index, id2_imputation, id2_imputed, imputation_prob(sym+num), imputed_prob(sym+num)]
# at this point imputationsTable holds all the data needed for E-M 

# -- guess initial point
Q_init <- c(0.3, 0.5, 0.7)


# -- expectation phase
expectationTable <- expectation_phase(imputationsTable, Q_init)



# Need to choose step size, tolerance and precision
Q_min <- c(0.2, 0.4, 0.4)
Q_max <- c(0.6, 0.9, 0.7)
tolerance <- 0.0001
precision <- 0.001
barsNum <- 6
B=FALSE

param_names <- if (isFALSE(B)) paste0('Q.', 1:risk_groups) else as.character(c(paste0('Q.', 1:risk_groups), paste0('B.', 1:risk_groups)))

maxim <- maximization_phase(expectationTable, tolerance, precision, barsNum, lower=Q_min, upper=Q_max, B=FALSE)
maximizer <- maxim[param_names]
maximal_value <- maxim$value






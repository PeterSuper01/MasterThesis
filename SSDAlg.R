library(pracma)
library(parallel)

J_function <- function(design, t = 3) {
  sum_products <- 0
  m <- ncol(design)
  comb <- combn(m, t)
  for (i in 1:choose(m, t)) {
    sum_products <- sum_products + sum(apply(design[, comb[, i]], 1, prod))^2
  }
  return(sum_products / choose(m, t))
}

col_exchange <- function(M, c1, c2) {
  M_temp <- M
  M_temp[, c1] <- M[, c2]
  M_temp[, c2] <- M[, c1]
  return(M_temp)
}

row_exchange <- function(M, r1, r2) {
  M_temp <- M
  M_temp[r1, ] <- M[r2, ]
  M_temp[r2, ] <- M[r1, ]
  return(M_temp)
}

CreateSSD <- function(N, n_factor, lambda, balance = T, n_try = 100) {
  # N <- 12
  # n_factor <- 20
  # lambda <- 0 # 0~1
  
  design <- cbind(rep(c(-1, 1), each=N/2), rep(c(-1, 1, 1, -1), each=N/4))
  for (column in 1:(n_factor-2)) {
    # for (column in 1:(n-2)) {
    last_added <- Inf
    for (try in 1:n_try) {
      # randomly generate a balance length N vector
      if (balance) {
        col_added <- rep(c(-1, 1), each=(N/2))[randperm(N, N)]
      } else {
        col_added <- rbinom(N, 1, 0.5)*2-1
        if (abs(sum(col_added)) == N)
          col_added[1] <- -col_added[1]
      }
      delta <- Reduce("+", apply(design, 2, function(x) (x %*% t(x)+1)/2, simplify = F))
      delta_added <- (col_added %*% t(col_added) + 1) / 2
      n <- 2+column
      
      delta_new <- delta + delta_added
      delta_for_K <- delta_new[upper.tri(delta_new)]
      K1 <- sum(delta_for_K)
      K2 <- sum(delta_for_K^2)
      K3 <- sum(delta_for_K^3)
      # lambda * E(s^2) + (1-lambda) * E(J_3^2)
      target <- (1-lambda)*8/N/N/(N-2) * K3 + 
        (4*(N-2)*lambda-12*n*(1-lambda))/N/N/(N-2) * K2 + 
        ((6*n^2-6*n+4)*(1-lambda)-4*n*(N-2)*lambda)/N/N/(N-2) * K1
      
      while (T) {
        
        # exchange every pair of different-sign row in col_added
        pos_place <- which(col_added==1)
        neg_place <- which(col_added==-1)
        target_exchage <- matrix(0, length(pos_place), length(neg_place))
        for (i in 1:length(pos_place)) {
          for (j in 1:length(neg_place)) {
            delta_added_temp <- col_exchange(delta_added, pos_place[i], neg_place[j])
            delta_added_temp <- row_exchange(delta_added_temp, neg_place[j], pos_place[i])
            # trans_matrix <- diag(N)
            # trans_matrix[c(pos_place[i], neg_place[j]), ] <- 
            #   trans_matrix[c(neg_place[j], pos_place[i]), ]
            # delta_added_temp <- trans_matrix %*% delta_added %*% t(trans_matrix)
            delta_new_temp <- delta + delta_added_temp
            delta_for_K <- delta_new_temp[upper.tri(delta_new_temp)]
            K1 <- sum(delta_for_K)
            K2 <- sum(delta_for_K^2)
            K3 <- sum(delta_for_K^3)
            target_exchage[i, j] <- (1-lambda)*8/N/N/(N-2) * K3 + 
              (4*(N-2)*lambda-12*n*(1-lambda))/N/N/(N-2) * K2 + 
              ((6*n^2-6*n+4)*(1-lambda)-4*n*(N-2)*lambda)/N/N/(N-2) * K1
          }
        }
        
        # if no improvement, then break
        last_value <- min(target_exchage)
        if (min(target_exchage) >= target) {
          break
        } else {
          # exchange the pair of rows that has the minimum target value
          exchange_places <- which(target_exchage == min(target_exchage), arr.ind=TRUE)[1, ]
          col_added[c(pos_place[exchange_places[1]], neg_place[exchange_places[2]])] <-
            col_added[c(neg_place[exchange_places[2]], pos_place[exchange_places[1]])]
          # update target value
          delta_added <- col_exchange(delta_added, pos_place[exchange_places[1]], neg_place[exchange_places[2]])
          delta_added <- row_exchange(delta_added, neg_place[exchange_places[2]], pos_place[exchange_places[1]])
          delta_new <- delta + delta_added
          delta_for_K <- delta_new[upper.tri(delta_new)]
          K1 <- sum(delta_for_K)
          K2 <- sum(delta_for_K^2)
          K3 <- sum(delta_for_K^3)
          # lambda * E(s^2) + (1-lambda) * E(J_3^2)
          target <- (1-lambda)*8/N/N/(N-2) * K3 + 
            (4*(N-2)*lambda-12*n*(1-lambda))/N/N/(N-2) * K2 + 
            ((6*n^2-6*n+4)*(1-lambda)-4*n*(N-2)*lambda)/N/N/(N-2) * K1
        }
        last_value <- min(target_exchage)
      }
      if (last_value < last_added) {
        last_added <- last_value
        best_col_added <- col_added
      }
    }
    
    # append design matrix by column
    design <- cbind(design, best_col_added)
  }
  colnames(design) <- 1:ncol(design)
  return(design)
}

CreateSSD1 <- function(N, n_factor, lambda, balance = T, n_try = 100) {
  design <- CreateSSD(N, n_factor+1, lambda, balance, n_try)
  target <- numeric(n_factor+1)
  for (i in 1:(n_factor+1)) {
    sub_design <- design[, -i]
    target[i] <- lambda*J_function(sub_design, 2) + (1-lambda)*J_function(sub_design, 3)
  }
  design <- design[, -which.min(target)]
  return(design)
}

# try 20 times CreateSSD1 to get the min target value
CreateSSD_min <- function(N, n_factor, lambda, balance = F, n_try = 1000, m_try = 20) {
  target_best <- Inf
  for (i in 1:m_try) {
    design <- CreateSSD1(N, n_factor, lambda, balance, n_try)
    target <- lambda*J_function(design, 2) + (1-lambda)*J_function(design, 3)
    if (target < target_best) {
      target_best <- target
      design_best <- design
    }
  }
  return(design_best)
}

lambda <- 0.2
# parallel version
cluster <- makeCluster(5)
clusterExport(cluster, c("col_exchange", "row_exchange", "CreateSSD", "CreateSSD1", "CreateSSD_min", "J_function", "lambda"))

result <- parLapply(cluster, 1:5, function(i) {
  set.seed(i)
  return(CreateSSD_min(12, 16, lambda, F, 1000, 20))
})
# result <- clusterCall(cluster, CreateSSD_min, 12, 16, lambda, F, 100, 20)
# clusterCall(cluster, function(x=2){x+1}, 5)
stopCluster(cluster)

target <- numeric(5)
for (i in 1:5) {
  target[i] <- lambda*J_function(result[[i]], 2)+(1-lambda)*J_function(result[[i]], 3)
}
best_design <- result[[which.min(target)]]
best_design

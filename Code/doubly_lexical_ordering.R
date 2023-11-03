## Implementation of doubly lexical matrix ordering algorithm as proposed by
## Hoffman et al., 1985, https://doi.org/10.1137/0606070 as described in
## Spinrad, 2000, BOOK, Efficient Graph Representations, pg. 121
## - has O(i^2,j)
## Konrad Herbst, 2022-12-07

binarize_matrix <- function(M){
  return(matrix(as.numeric(!is.na(M)), ncol = ncol(M)))
}

row_vals <- function(M, R, C){
  x <- lapply(R, function(r){
    sapply(C, function(c){
      return(sum(M[r,c]))
    })
  })
  names(x) <- R
  return(x)
}

lexsort <- function(row_vals){
  rows <- as.numeric(names(row_vals))
  lex <- as.numeric(sapply(row_vals, function(x) paste0(rev(x), collapse = "")))
  return(rows[order(lex, decreasing = TRUE)])
}

N <- function(r, M){
  return(  which(M[r, ] == 1)  )
}

dlordering <- function(M){
  ## the time complexity of this algorithm is O(i^2, j) so it's favorable
  ## to use the smaller dimension as matrix rows
  transposed <- nrow(M)>ncol(M)
  if(transposed) M <- t(M)
  R <- 1:nrow(M)
  C <- list(1:ncol(M))
  Rret <- list()
  while(length(R)>0L){
    rowvals <- row_vals(M, R, C)
    partrow <- lexsort(rowvals)[[1]]
    C <- unlist(lapply(C, function(Cc){
      ret <- list()
      C_1 <- intersect(Cc, N(partrow, M))
      C_2 <- setdiff(Cc, C_1)
      if((length(C_2)>0L)){
        ret <- list(C_2)
      }
      if(length(C_1)>0L){
        ret <- c(ret, list(C_1))
      }
      return(ret)
    }), recursive = FALSE)
    
    Rret <- c(partrow, Rret)
    R <- setdiff(R, partrow)
  }
  if(!transposed){
    ret <- list("R" = unlist(Rret), "C" = unlist(C))
  } else {
    ret <- list("R" = unlist(C), "C" = unlist(Rret))
  }
  return(ret)
}


M <- matrix(c(
  1,0,1,1,0,
  0,1,1,0,0,
  1,1,1,0,0,
  0,1,0,1,1,
  1,0,0,1,0
), nrow = 5, byrow = TRUE)

dlo <- dlordering(M)
M
M[dlo[[1]], dlo[[2]]]

M2 <- matrix(c(
  0,1,0,1,1,0,
  1,0,1,1,1,0,
  0,1,1,1,1,0,
  1,1,0,0,0,1,
  0,1,0,1,0,1,
  1,0,0,1,1,1
), nrow = 6, byrow = TRUE)
dlo <- dlordering(M2)
M2
M2[dlo[[1]], dlo[[2]]]
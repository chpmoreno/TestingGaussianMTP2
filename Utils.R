# Function to determine if a Matrix is an M-Matrix
is_m_matrix <- function(K) {
      n_values <- prod(dim(K))
      if(sum(K <= 0) == (n_values - dim(K)[1])) {
            is_it = 1
      } else {
            is_it = 0
      }
      return(is_it)
}
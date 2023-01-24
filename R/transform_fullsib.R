#' Compute expected data conditional on unbalanced data
#'
#' @param y_data An object inheriting "`fullsibdata`" containing unbalanced data
#' @param params A list of objects containing a field called `means` used to fill
#' in the unobserved data of `y_data`
#'
#' @return A balanced `fullsibdata` object.
#' 
#' @export
balance_data <- function(y_data, params) {

  balanced_tables <- list()
  J <- y_data$n_ind

  for(i in 1:length(y_data$tables)) {

    # Fill the bottom of the matrix with repeated means
    conditional_means <- matrix(
      rep(params[[i]]$mean, times = J - nrow(y_data$tables[[i]])),
      ncol  = y_data$n_traits,
      byrow = TRUE
    )
    
    balanced_tables[[i]] <- rbind(
      y_data$tables[[i]],
      conditional_means
    )
  }

  return(fullsibdata(balanced_tables))
  
}

#' Compute sum-of-squares matrices
#'
#' @export
ss_oneway <- function(y_data) {
  centered_tables <- lapply(
    y_data$tables,
    \(X) scale(X, center = TRUE, scale = FALSE)
  )

  ind_mean_addends <- lapply(
    centered_tables,
    \(X) t(X) %*% X
  )

  ind_ss <- Reduce("+", ind_mean_addends)

  sires_means <- t(sapply(
    y_data$tables,
    colMeans
  ))

  sire_means_centered <- scale(sires_means,
                               center = TRUE,
                               scale = FALSE)
  
  sire_ss <- y_data$n_ind * t(sire_means_centered) %*% sire_means_centered

  list(
    A_E = ind_ss,
    A_A = sire_ss
  )
}

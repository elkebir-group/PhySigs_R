#' Normalize mutational signature feature matrix
#'
#' This function takes as input feature matrix where each row is a clone and each
#' column is a trinucleotide context (e.g., "A(C > A)G"). It then adjusts the counts
#' based on the trinucleotide context fraction. The backend of this function relies
#' getTriContextFraction in deconstuctSigs. Possible normalization options include:
#' on "default", "genome", and "exome". Please see deconstructSigs documentation
#' for more information on the normalization procedure.
#'
#' @param feat_mat data.frame
#' @param norm_method character
#'
#' @return data.frame
#' @export
#'
#' @examples
#' normalizeFeatureMatrix(P, "genome")
normalizeFeatureMatrix <- function(feat_mat, norm_method) {

  # For each node in the tree
  for (row in row.names(feat_mat)) {

    # Update row with normalized trinucleotide counts
    # Note that we multiply by the sum of the row to keep counts
    feat_mat[row, ] <- deconstructSigs::getTriContextFraction(
      mut.counts.ref = feat_mat[row, ],
      trimer.counts.method = norm_method
    ) * sum(feat_mat[row, ])
  }

  # Return updated data frame
  return(feat_mat)
}

#' Get matrix factorization reconstruction error
#'
#' Given a normalized feature matrix, estimated exposures, and signature matrix, this function returns
#' the Frobineous norm between the reconstructed feature matrix and the normalized feature matrix. In
#' the feature matrix, each row is a clone and each column is a trinucleotide context. In the exposure
#' matrix, each row corresponds to a mutational signature and each column corresponds to a cluster of
#' clones. Finally, in the signature vector, the signatures of interest are given (e.g., "Signature.1").
#'
#' @param feat_mat data.frame
#' @param exp_mat data.frame
#' @param sigs_filter character vector
#'
#' @return numeric
#' @export
#'
#' @examples getError(P_Norm, E_list[[1]], S)
getError <- function(feat_mat, exp_mat, sigs_filter) {

  # Expand node clusters into a column for each node
  expanded_exp_mat <- expandClusters(exp_mat)

  # Ensure the node exposures are also in the feature matrix
  expanded_exp_mat <- expanded_exp_mat[row.names(feat_mat)]

  # Change exposure matrix to be in terms of counts not percentage
  for (node in row.names(feat_mat)) {
    expanded_exp_mat[node] <- expanded_exp_mat[node] * sum(feat_mat[node, ])
  }

  # Compute a difference matrix where each entry computes the differnce between the original feature counts
  # and the counts reconstructed after multiplying the exposure to mutational signatures
  diff_mat <- feat_mat - (t(as.matrix(expanded_exp_mat)) %*% as.matrix(deconstructSigs::signatures.cosmic[sigs_filter, ]))

  # Take the sum of squares for error
  error <- sum(abs(diff_mat)^2)

  # Return reconstruction error
  return(error)
}

#' Compute the Bayesian Information Criterion
#'
#' Given a normalized feature matrix, estimated exposures, and signature matrix, this function returns
#' the Bayesian information criterion based on the error between the reconstructed feature matrix and
#' the normalized feature matrix. In the feature matrix, each row is a clone and each column is a
#' trinucleotide context. In the exposure matrix, each row corresponds to a mutational signature and
#' each column corresponds to a cluster of clones. Finally, in the signature vector, the signatures of
#' interest are given (e.g., "Signature.1").
#'
#' @param feat_mat data.frame
#' @param exp_mat data.frame
#' @param sigs_filter character vector
#'
#' @return numeric
#' @export
#'
#' @examples getBIC(P_Norm, E_list[[1]], S)
getBIC <- function(feat_mat, exp_mat, sigs_filter) {

  # Get the reconstruction error between the feature matrix
  # and the exposures multiplied by the signature matrix
  error <- getError(feat_mat, exp_mat, sigs_filter)

  # Get the number of elements of feature matrix
  # (i.e., # nodes * 96 features)
  n <- prod(dim(feat_mat))

  # Get the number of elements of exposure matrix
  # (i.e. # active signatures * # exposure shifts)
  k <- prod(dim(exp_mat))

  # Return BIC
  return(n * log(error / n) + k * log(n))
}

#' Compute best exposure for k exposure shifts
#'
#' This function takes as input a phylogeny, the signature vector, and the normalized feature matrix.
#' It also takes as input k, the desired number of exposure shifts to be found in the input phylogeny.
#' The phylogeny is represented with a graphNel object, and the normalized feature matrix has rows
#' corresponding to tree nodes and columns corresponding to trinucleotide contexts (e.g., "A(C > A)G").
#' It returns a list of exposure matrices whose entries correspond to the best partition of the
#' phylogeny for each possible number of clusters. Each exposure matrix has a row for each requested
#' signature and a column for each cluster of clones induced by the partition of the phylogeny.
#'
#' @param tree graphNel
#' @param feat_mat data.frame
#' @param k numeric
#' @param sigs_filter character vector
#'
#' @return data.frame
#' @export
#'
#' @examples treeExposures(T_tree, P_Norm, 1, S)
treeExposures <- function(tree, feat_mat, k, sigs_filter) {

  # Get nodes and edges of tumor tree
  V <- graph::nodes(tree)
  E <- graph::edgeMatrix(tree)

  # Reformat edges into list of tuples
  E_list <- list()
  nrEdges <- dim(E)[[2]]
  for (i in 1:nrEdges) {
    E_list[[i]] <- cbind(E[1, i], E[2, i])
  }

  # Generate all combination of k edges to create all possible
  # expsoure shifts on tree
  C <- utils::combn(E_list, k)

  # Initialize best clustering
  best_CC <- NULL
  best_error <- Inf
  best_sample_exp <- NULL

  # For each clustering
  for (idx in 1:dim(C)[2]) {

    # Copy over input tree and remove selected edges for this clustering
    # The remaining connected components are then the node clusters
    cpy_tree <- tree
    if (k > 0) {
      for (i in 1:k) {
        edge <- C[i, idx][[1]]
        cpy_tree <- graph::removeEdge(V[edge[1]], V[edge[2]], cpy_tree)
      }
    }

    # Initialize exposure data frame with columns labeled by active signatures
    sample_exp <- data.frame(matrix(0L, nrow = length(sigs_filter), ncol = 0))
    row.names(sample_exp) <- sigs_filter

    # Identify node clusters from connected components
    CC <- graph::connComp(cpy_tree)
    error <- 0

    # For each connected component
    for (CCC in CC) {

      # We reduce the TE problem to the SE problem
      feat_mat_CCC <- sum(feat_mat[CCC[1], ]) * feat_mat[CCC[1], ]
      if (length(CCC) >= 2) {
        for (i in 2:length(CCC)) {
          feat_mat_CCC <- feat_mat_CCC + (sum(feat_mat[CCC[i], ]) * feat_mat[CCC[i], ])
        }
      }
      row.names(feat_mat_CCC) <- paste(CCC, collapse = ";")
      feat_mat_CCC <- as.data.frame(feat_mat_CCC)

      # Get exposure for this connected component
      sample_exp_CCC <- deconstructSigs::whichSignatures(
        tumor.ref = feat_mat_CCC,
        signatures.ref = deconstructSigs::signatures.cosmic,
        associated = sigs_filter,
        contexts.needed = TRUE,
        signature.cutoff = 0.0001,
        tri.counts.method = "default"
      )

      # Add unknown signature if found by deconstructSigs
      sample_exp_CCC$weights$Signature.unknown <- sample_exp_CCC[["unknown"]]
      active <- sample_exp_CCC$weights[sigs_filter]
      sample_exp <- cbind(sample_exp, t(active))
    }

    # Compute reconstruction error between feature matrix and exposures times signatures
    error <- getError(feat_mat, sample_exp, sigs_filter)

    # If the error is lower than previously observed error, keep this clustering
    if (error < best_error) {
      best_error <- error
      best_CC <- CC
      best_sample_exp <- sample_exp
    }
  }

  # Print the minimum reconstruction error for k exposure shifts
  print(paste("k:", k, "; error:", best_error))

  # Return best exposures
  return(best_sample_exp)
}

#' Compute all tree exposures for different number of clusters
#'
#' This function takes as input a phylogeny, the signature vector, and the normalized feature matrix.
#' The phylogeny is represented with a graphNel object, and the normalized feature matrix has rows
#' corresponding to tree nodes and columns corresponding to trinucleotide contexts (e.g., "A(C > A)G").
#' It returns a list of exposure matrices whose entries correspond to the best partition of the
#' phylogeny for each possible number of clusters. Each exposure matrix has a row for each requested
#' signature and a column for each cluster of clones induced by the partition of the phylogeny.
#' Start and stop are optional parameters to specify the number of exposure shifts to try.
#'
#' @param tree graphNel
#' @param feat_mat data.frame
#' @param sigs_filter character vector
#' @param start numeric
#' @param stop numeric
#'
#' @return list of data.frame objects
#' @export
#'
#' @examples allTreeExposures(T_tree, P_Norm, S, 0, 1)
allTreeExposures <- function(tree, feat_mat, sigs_filter, start = NULL, stop = NULL) {

  # Initialize list of exposures
  exp_list <- list()

  # Set range of k to loop over
  if (is.null(start)){
    start = 0
  }

  if (is.null(stop)){
    stop = length(graph::nodes(tree)) - 1
  }

  # Loop over all possible number of exposure shifts
  for (k in start:stop) {

    # Find best exposures for fixed number of clusters k
    exp_list[[as.character(k)]] <- treeExposures(tree, feat_mat, k, sigs_filter)
  }

  # Return list of exposures found for each k
  return(exp_list)
}

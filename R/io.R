#' Format exposure output data frame for given number of exposure shifts
#'
#' Given a sample ID (e.g., patient tumor ID), tree ID for this sample (e.g., 2), list of exposure
#' data frames for all numbers of exposure shifts, and the number of exposure shifts of interest,
#' this function returned a dataframe where each row corresponds to a node in the sample tree. The
#' columns are as follows: Sample ID, Tree ID, Node ID, Signature 1 exposure ... Signature X exposure.
#' This formatting is used in the input to the corresponding visualization tool for PhySigs.
#'
#' @param sample character
#' @param tree_idx numeric
#' @param exp_list list of data.frame objects
#' @param k numeric
#'
#' @return data.frame
#' @export
#'
#' @examples outputExposures("Tumor1", 0, E_list, 1)
outputExposures <- function(sample, tree_idx, exp_list, k){

  # Initialize data frame
  out_df <- data.frame("Sample"=character())

  # For cluster in clustering with k exposure shifts
  for (c in names(exp_list[[k]])){

    # Obtain exposures to mutational signatures for this node cluster
    current_exp <- exp_list[[k]][c]

    # For each node in cluster
    for (n in unlist(strsplit(c, "[;]"))){

      # Add row to dataframe showing the exposure for this node in this patient tree
      x <- data.frame("Sample" = sample, "Tree_ID" = tree_idx, "k" = k, "Node" = n)
      x <- cbind(x,t(current_exp))
      out_df <- fastmerge(out_df, x)
    }
  }
  # Clear row names
  rownames(out_df) <- c()

  # Output data frame with the following columns
  # Sample; Tree ID (1-based index); Node (Clone; 1-based index); Signature 1 exposure ... Signature X exposure
  return(out_df)
}

#' Format tree output data frame
#'
#' Given a sample ID (e.g., patient tumor ID), tree ID for this sample (e.g., 2), and a tree object,
#' This function returns a dataframe with one row corresponding to the input tree. The columns are
#' as follows: Sample ID; Tree ID; Number of clones/nodes; Nodes; Edges (source -> target).
#' This formatting is used in the input to the corresponding visualization tool for PhySigs.
#'
#' @param sample character
#' @param tree_idx numeric
#' @param tree graphNel
#'
#' @return data.frame
#' @export
#'
#' @examples outputTrees("Tumor1", 0, T_tree)
outputTrees <- function(sample, tree_idx, tree){

  # Get nodes and edges of tree
  V <- graph::nodes(tree)
  E <- graph::edgeMatrix(tree)

  # Make list of noes separated by semicolon
  nrNodes <- graph::numNodes(tree)
  for (i in 1:nrNodes) {
    if (i==1){
      v <- c(V[i])
    }
    else{
      v <- paste(v, V[i], sep = ";")
    }
  }

  # Make list of edges separated by semicolon
  # Each edge is formatted source arrow target
  nrEdges <- dim(E)[[2]]
  for (i in 1:nrEdges) {
    source <- V[[E[1,i]]]
    target <- V[[E[2,i]]]
    if (i==1){
      e <- paste(source, target, sep="->")
    }
    else{
      e <- paste(e, paste(source, target, sep="->"), sep = ";")
    }
  }

  # Output data frame with the following columns
  # Sample ID; Tree ID (1-based index); Number of clones/nodes (integer); Nodes; Edges (source -> target)
  return(data.frame("Sample" = sample, "Tree_ID" = tree_idx, "Num_Nodes" = nrNodes, "Nodes" = v, "Edges" = e))
}

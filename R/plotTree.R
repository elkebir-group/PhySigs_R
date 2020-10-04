#' Plot tree with exposure pie chart at each node
#'
#' Given a tree represented by a graphNel object, this function prints the tree with each
#' node labeled by a pie chart showing the exposure to mutational signatures at that node.
#' Patient, title, and tree_idx fields just label the figure. The exposure matrix should have
#' a row for each active mutational signature and a column for each clone cluster.
#'
#' @param patient character
#' @param title character
#' @param tree graphNel
#' @param exp_mat data.frame
#' @param tree_idx numeric
#'
#' @export
#'
#' @examples plotTree("Tumor1", "Example Plot", T_tree, E_list[[1]], tree_idx=1)
plotTree <- function(patient, title, tree, exp_mat, tree_idx = 0) {

  if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
    stop("Package \"Rgraphviz\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Expand matrix with clone cluster columns into matrix with clone columns
  expanded_exp_mat <- expandClusters(exp_mat)

  # Subset to clones in present tree (in case some have been trimmed out)
  expanded_exp_mat <- expanded_exp_mat[, graph::nodes(tree)]

  # Initialize edge attributes for tree graphic
  eAttrs <- list()
  eAttrs$color <- list()

  # Get vertex and edge sets from tree
  V <- graph::nodes(tree)
  E <- graph::edgeMatrix(tree)
  nrEdges <- dim(E)[[2]]

  # Create attributes for each edge
  for (i in 1:nrEdges) {
    source <- V[[E[1, i]]]
    target <- V[[E[2, i]]]

    ok <- FALSE

    # Find if source and target nodes are in the same cluster
    for (C in names(exp_mat)) {
      CC <- strsplit(C, "[;]")[[1]]

      if ((source %in% CC) && (target %in% CC)) {
        ok <- TRUE
      }
    }

    # Change formatting of edge depending on if there is an exposure shift
    if (ok) {
      eAttrs$color[[paste(source, "~", target, sep = "")]] <- "black"
    }
    else {
      eAttrs$color[[paste(source, "~", target, sep = "")]] <- "black"
    }
  }

  # Set color palette
  fill <- RColorBrewer::brewer.pal(8, "Set1")
  names(fill) <- row.names(exp_mat)

  # Format pie chart nodes
  g1layout <- Rgraphviz::agopen(tree, name = "foo")
  makeNodeDrawFunction <- function(x, fill, patient) {
    force(x)
    function(node, ur, attrs, radConv) {

      # Remove labels
      names(x) <- vector(mode = "character", length = length(names(x)))

      # Get node locations
      nc <- Rgraphviz::getNodeCenter(node)

      # Get consistent color scheme
      pal <- fill[which(x > 0)]

      # Make plots
      Rgraphviz::pieGlyph(x[which(x > 0)], xpos = Rgraphviz::getX(nc), ypos = Rgraphviz::getY(nc), radius = Rgraphviz::getNodeRW(node), col = pal)
    }
  }
  drawFunc <- apply(expanded_exp_mat, 2, makeNodeDrawFunction, fill = fill, patient = patient)

  # Make plot title
  title_id <- paste(patient, sep = "")
  if (tree_idx > 0) {
    title_id <- paste(patient, letters[tree_idx], sep = "")
  }

  # Truncate signature names
  sig_ids <- strsplit(as.character(row.names(exp_mat)), "[.]")
  sig_ids <- as.character(sapply(sig_ids, "[", 2))

  # Make plot
  Rgraphviz::plot(tree,
    drawNode = drawFunc, edgeAttrs = eAttrs,
    attrs = list(
      node = list(height = 2, width = 2, fontsize = 2),
      edge = list(fontsize = 5)
    ), mai = c(0.15, 0.15, 0.15, 1),
    main = paste("ID:", title_id, "--", title, sep = " ")
  )

  # Make plot legend
  graphics::legend("topright", inset = c(-0.09, .03), title = "Signatures", sig_ids, fill = fill, xpd = TRUE)
}

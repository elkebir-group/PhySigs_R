#' Concatenate two dataframes
#'
#' Create a new dataframe from two existing dataframes by stacking them. If, a column
#' exists in one dataframe but not the other, NA is added. This code was adapted
#' from https://stackoverflow.com/questions/8169323/r-concatenate-two-dataframes.
#'
#' @param d1 data.frame
#' @param d2 data.frame
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fastmerge(data.frame("X" = 1:2, "Y" = c("a", "b")), data.frame("X" = 1:2, "Z" = c(5, 7)))
fastmerge <- function(d1, d2) {

  # Get column names
  d1.names <- names(d1)
  d2.names <- names(d2)

  # Columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)

  # Columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)

  # Add blank columns to d2
  if (length(d2.add) > 0 & nrow(d2) > 0) {
    for (i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }

  # Add blank columns to d1
  if (length(d1.add) > 0 & nrow(d1) > 0) {
    for (i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }

  # Return concatenated data frame
  return(rbind(d1, d2))
}

#' Expand exposure dataframe width by spliting semicolons
#'
#' Given a dataframe with column names containing semicolons,
#' splits column name at semicolon. The resulting string split generates new column names.
#' Each new column name is populated with a duplicate of its original column.
#' This is used to expand columns labeled by clusters of nodes into a column for each node.
#'
#' @param d1 data.frame
#'
#' @return data.frame
#' @export
#'
#' @examples
#' x <- data.frame(c(.2, .3, .5), c(.1, .7, .2))
#' colnames(x) <- c("a;b;c", "1;2;3")
#' expandClusters(x)
expandClusters <- function(d1) {

  # Create empty data frame to return with same number of rows and row names
  d2 <- data.frame(matrix(0L, nrow = length(row.names(d1)), ncol = 0))
  row.names(d2) <- row.names(d1)

  # Split the name of each column at the semicolons
  for (col in names(d1)) {
    s <- unlist(strsplit(col, ";"))

    # Make a copy of the column for every node in the column name
    for (node in s) {
      d2[as.character(node)] <- d1[col]
    }
  }

  # Return expanded data frame
  return(d2)
}

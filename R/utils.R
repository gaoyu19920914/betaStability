# function: get the n-th row of a data frame with error handling
get_nth_row <- function(df, n) {
    # Check if input is a data frame
    if (!is.data.frame(df)) {
        stop("Error: 'df' must be a data frame.")
    }

    # Check if n is numeric and integer-like
    if (!is.numeric(n) || length(n) != 1 || n != as.integer(n)) {
        stop("Error: 'n' must be a single integer.")
    }

    # Convert to integer in case it's numeric but whole number
    n <- as.integer(n)

    # Check bounds
    if (n < 1 || n > nrow(df)) {
        stop(
            "Error: Row number ",
            n,
            " is out of bounds. Data frame has ",
            nrow(df),
            " rows."
        )
    }

    # Return the n-th row as a data frame
    return(df[n, , drop = FALSE])
}

# function: calculate the absolute difference between the n-th rows
get_deltaenv_rows <- function(x, y, env) {
    row_x <- get_nth_row(env, x)
    row_y <- get_nth_row(env, y)
    if (is.null(row_x) || is.null(row_y)) {
        return(NULL)
    }
    return(abs(row_x - row_y))
}

# function: convert the predicted distance vector to a distance matrix
pred2matrix <- function(pred, site_ids) {
    n_sites <- length(site_ids)
    if (n_sites * (n_sites - 1) / 2 != length(pred)) {
        stop("The dimensions of variables are not matched!")
    }

    # Create an empty distance matrix
    dist_matrix <- matrix(0, nrow = n_sites, ncol = n_sites)
    rownames(dist_matrix) <- site_ids
    colnames(dist_matrix) <- site_ids
    n <- 1
    for (i in seq_len(n_sites - 1)) {
        for (j in (i + 1):n_sites) {
            if (n > length(pred)) stop("length of pred out of bourder")
            dist_matrix[i, j] <- pred[n]
            dist_matrix[j, i] <- pred[n]
            n <- n + 1
        }
    }
    return(dist_matrix)
}

# function: clean dataframe (remove same values and convert strings to integers)
df_prepare <- function(df) {
    # Remove columns where all elements are the same
    df <- df[, vapply(df, function(x) length(unique(x)) > 1, logical(1))]

    # Convert columns to integer if they contain different types of strings
    for (col in names(df)) {
        if (is.character(df[[col]]) && length(unique(df[[col]])) > 1) {
            df[[col]] <- as.integer(as.factor(df[[col]]))
        }
    }
    return(df)
}

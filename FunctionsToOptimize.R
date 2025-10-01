
treeBorders = function(tree){
  # --- Setup ---
  splits <- tree$split
  dimension <- tree$dimension
  n_splits <- nrow(splits)

  # Handle the case of no splits
  if (n_splits == 0) {
    Bounds <- data.frame(matrix(tree$border, nrow = 1))
    names(Bounds) <- c(str_c("L", 1:dimension), str_c("U", 1:dimension))
    return(Bounds)
  }

  # --- 1. Pre-allocate a matrix ---
  # The maximum possible number of final boundaries is 1 + the number of splits.
  max_rows <- 1 + n_splits
  # Using a matrix is much faster for iterative assignment than a data frame.
  bounds_mat <- matrix(NA_real_, nrow = max_rows, ncol = 2 * dimension)

  # Initialize with the starting border
  bounds_mat[1, ] <- as.vector(tree$border)

  # Use a pointer to track the next available empty row
  next_row_idx <- 2

  colnames(splits) = c('splits','dim','split_part')

  # --- 2. Loop and Fill ---
  # Iterate through the splits and fill the pre-allocated matrix

  for (i in 1:n_splits) {
    # Extract split info once per loop
    cur_split_row <- splits[i, ]
    cur_part <- cur_split_row$split_part
    cur_dim <- cur_split_row$dim
    new_bound <- cur_split_row$splits

    # Check for valid partition index (prevents errors)
    if (cur_part >= next_row_idx) next

    # Define column indices for clarity and speed
    lower_col <- cur_dim
    upper_col <- cur_dim + dimension

    # Get current boundaries from our matrix
    cur_L <- bounds_mat[cur_part, lower_col]
    cur_U <- bounds_mat[cur_part, upper_col]

    # Check if the split is valid
    if (new_bound > cur_L && new_bound < cur_U) {
      # Create the new sibling row by copying the parent's data
      sibling_row <- bounds_mat[cur_part, ]

      # Modify the parent row (becomes the 'lower' partition)
      bounds_mat[cur_part, upper_col] <- new_bound

      # Modify the sibling row (becomes the 'upper' partition)
      sibling_row[lower_col] <- new_bound

      # Place the new sibling row in the next available spot (NO rbind!)
      bounds_mat[next_row_idx, ] <- sibling_row

      # Increment the row pointer
      next_row_idx <- next_row_idx + 1
    }
  }

  # --- 3. Finalize ---
  # Trim the matrix to the number of rows we actually used
  final_bounds <- bounds_mat[1:(next_row_idx - 1), , drop = FALSE]

  # Convert to a data frame and assign names once at the end
  final_df <- as.data.frame(final_bounds)
  names(final_df) <- c(str_c("L", 1:dimension), str_c("U", 1:dimension))

  return(final_df)

}

createNewRegion = function(tree, random = T, split_part, split_dim, split_perc){

  newTree = tree

  cur_Borders = treeBorders(tree)

  if(random){
    split_part = sample(1:nrow(cur_Borders),size = 1,
                        prob = apply(apply(cur_Borders, MARGIN = 1, FUN = diff, lag = tree$dimension), MARGIN = 2, FUN = prod)/prod(diff(tree$border, lag = tree$dimension)))
    split_dim = sample(1:tree$dimension,size = 1, prob = diff(unlist(cur_Borders[split_part,]), lag = tree$dimension)/sum(diff(unlist(cur_Borders[split_part,]), lag = tree$dimension)))
    split_perc = rbeta(1,5,5)
  }
  if(split_part > nrow(cur_Borders)){
    print("Bruh there aren't that many partitions")
  } else if(split_dim > tree$dimension){
    print("Bruh the tree doesn't have that many possible predictors")
  } else if(split_perc >= 1 | split_perc <= 0){
    print("Not a valid percentage: Value must be in (0,1)")
  } else{



    cur_boundary = as.numeric(cur_Borders[split_part,c(split_dim,split_dim + tree$dimension)])
    new_boundary = split_perc*(cur_boundary[2] - cur_boundary[1]) + cur_boundary[1]

    new_splits = c(tree$split$splits, new_boundary)
    new_dim = c(tree$split$dim, split_dim)
    new_split_part = c(tree$split$split_part, split_part)

    newTree$split = data.frame(splits = new_splits, dim = new_dim, split_part = new_split_part)

    newTree
  }

}

generateTree <- function(n_splits, border) {
  # --- 1. Initialization ---
  dimension <- length(border) / 2

  # Handle the edge case of no splits
  if (n_splits == 0) {
    return(create_tree(
      boundaries = integer(0),
      split_dim_v = integer(0),
      split_label_v = integer(0),
      dim = dimension,
      border = border
    ))
  }

  # Pre-allocate a data frame for storing the split history. This is crucial for speed.
  splits_df <- data.frame(
    splits = numeric(n_splits),
    dim = integer(n_splits),
    split_part = integer(n_splits)
  )

  # Pre-allocate a matrix to store the boundaries of the current leaf nodes.
  # This avoids costly resizing (like rbind) inside the loop.
  max_regions <- n_splits + 1
  current_regions <- matrix(0.0, nrow = max_regions, ncol = 2 * dimension)
  current_regions[1, ] <- as.vector(border)

  # n_regions tracks the number of active leaf regions we have.
  n_regions <- 1

  # Define column indices once for faster access
  lower_cols <- 1:dimension
  upper_cols <- (dimension + 1):(2 * dimension)

  # --- 2. Main Loop ---
  # This loop builds the tree structure efficiently.

  for (i in 1:n_splits) {

    # Get the matrix of regions that are currently active leaves
    active_regions_mat <- current_regions[1:n_regions, , drop = FALSE]

    # Calculate hyper-volumes to use as sampling weights
    side_lengths <- matrix(active_regions_mat[, upper_cols] - active_regions_mat[, lower_cols], ncol = dimension)
    volumes <- apply(side_lengths, 1, prod)

    # Sample the index of the region we are going to split
    region_to_split_idx <- sample.int(n = n_regions, size = 1, prob = volumes)

    # Get the side lengths of only the chosen region
    chosen_region_sides <- side_lengths[region_to_split_idx, ]

    # Sample the dimension to split along, weighted by the side lengths
    dim_to_split <- sample.int(n = dimension, size = 1, prob = chosen_region_sides)

    # Generate a random split point using a Beta distribution
    split_perc <- rbeta(1, 5, 5)

    # Get the bounds of the chosen region along the chosen dimension
    cur_boundary <- c(
      current_regions[region_to_split_idx, lower_cols[dim_to_split]],
      current_regions[region_to_split_idx, upper_cols[dim_to_split]]
    )

    new_boundary_val <- split_perc * (cur_boundary[2] - cur_boundary[1]) + cur_boundary[1]

    # Record the split information in our pre-allocated data frame
    splits_df[i, ] <- c(new_boundary_val, dim_to_split, region_to_split_idx)

    # --- Update the Regions Matrix (No Recalculation!) ---

    # The new sibling region starts as a copy of its parent
    sibling_region <- current_regions[region_to_split_idx, ]

    # Modify the parent region (it becomes the "lower" part of the split)
    current_regions[region_to_split_idx, upper_cols[dim_to_split]] <- new_boundary_val

    # Modify the sibling region (it becomes the "upper" part of the split)
    sibling_region[lower_cols[dim_to_split]] <- new_boundary_val

    # Increment the region counter and place the new sibling in the next available row
    n_regions <- n_regions + 1
    current_regions[n_regions, ] <- sibling_region
  }

  # --- 3. Finalization ---
  # Use the collected splits to create the final tree object in your desired format
  final_tree <- create_tree(
    boundaries = splits_df$splits,
    split_dim_v = splits_df$dim,
    split_label_v = splits_df$split_part,
    dim = dimension,
    border = border
  )

  return(final_tree)
}










n_splits = 99
border = c(0,0,1,1)


fast_tree = generate_tree_fast(999, c(0,0,1,1))

plotTreeGrid(fast_tree)

#

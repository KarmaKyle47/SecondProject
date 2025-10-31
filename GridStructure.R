library(stringr)
library(ggplot2)

create_tree = function(boundaries, split_dim_v, split_label_v, dim, border){
  list(dimension = dim,
       split = data.frame(splits = boundaries, dim = split_dim_v, split_part = split_label_v),
       border = border)
}

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

split_domain = function(grid, tree, plot = T){
  grid = data.frame(grid)
  label = rep(1, nrow(grid))

  j=1
  if(nrow(tree$split) > 0){
    for(i in 1:(nrow(tree$split))){
      cur_tree = tree
      cur_tree$split = cur_tree$split[0:(j-1),]
      cur_Borders = treeBorders(cur_tree)
      cur_part = tree$split$split_part[i]
      cur_dim = tree$split$dim[i]
      split_boundary = tree$split$splits[i]

      if(cur_part <= j & cur_dim <= tree$dimension){
        if(split_boundary > cur_Borders[cur_part, cur_dim] & split_boundary < cur_Borders[cur_part, cur_dim + tree$dimension]){
          j=j+1
          cur_grid = grid[label == cur_part,cur_dim]

          label[label == cur_part][cur_grid > split_boundary] = j
        }
      }
    }
  }



  if(plot){
    grid.label = data.frame(cbind(grid,label))
    grid.label$label = as.factor(grid.label$label)
    Borders = treeBorders(tree)

    plot = ggplot(grid.label, aes(x = X1, y = X2, color = label)) + geom_point() +
      geom_segment(data = Borders, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
      geom_segment(data = Borders, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
      geom_segment(data = Borders, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
      geom_segment(data = Borders, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') #+
      #xlim(0,1) + ylim(0,1)

    list(labels = label,plot = plot)
  } else{
    label
  }
}

plotTreeGrid = function(tree){

  if(tree$dimension != 2){
    warning("Grid can only be plotted in 2D")
  } else{

    Borders = treeBorders(tree)
    ShadeData = data.frame(xmin = Borders$L1,
                           xmax = Borders$U1,
                           ymin = Borders$L2,
                           ymax = Borders$U2,
                           Region = str_c("Region", 1:nrow(Borders)))

    plot = ggplot() +
      geom_segment(data = Borders, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
      geom_segment(data = Borders, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
      geom_segment(data = Borders, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
      geom_segment(data = Borders, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
      geom_rect(data = ShadeData,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Region),
                inherit.aes = FALSE,
                alpha = 0.4) +
      xlab("Longitude") + ylab("Latitude") + theme(legend.position = "none")

    plot
  }


}

createNewRegion = function(tree, random = T, split_part, split_dim, split_perc){

  newTree = tree

  if(is.list(tree$split) & !is.data.frame(tree$split)){
    tree$split = data.frame(do.call(rbind,tree$split))
  }


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



    if(!is.data.frame(newTree$split)){
      newTree$split[[length(newTree$split)+1]] = c(new_boundary, split_dim, split_part)

    } else{
      new_splits = c(tree$split$splits, new_boundary)
      new_dim = c(tree$split$dim, split_dim)
      new_split_part = c(tree$split$split_part, split_part)

      newTree$split = data.frame(splits = new_splits, dim = new_dim, split_part = new_split_part)
    }
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
    split_perc <- rbeta(1,10,10)
    #rbeta(5,5)

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

generate_grid_tree <- function(grid_size, border) {
  # --- 1. Validation and Initialization ---
  if (length(grid_size) != 1 || !is.numeric(grid_size) || grid_size <= 0) {
    stop("'grid_size' must be a single positive number.")
  }

  dimension <- length(border) / 2
  if (dimension %% 1 != 0) {
    stop("Length of 'border' must be an even number.")
  }

  # --- Calculate grid counts from grid_size ---
  extents <- border[(dimension + 1):(2 * dimension)] - border[1:dimension]
  if (any(extents <= 0)) {
    stop("Border extents must be positive (e.g., xmax > xmin).")
  }

  # Calculate the number of cells in each dimension
  num_cells <- extents / grid_size

  # Check if the extents are clean multiples of grid_size
  if (any(abs(num_cells - round(num_cells)) > 1e-9)) {
    warning("Border extents may not be integer multiples of grid_size. \nRounding to the nearest number of cells.")
  }
  grid_sizes <- as.integer(round(num_cells))

  if (any(grid_sizes < 1)) {
    stop("Calculated grid cells are less than 1 in at least one dimension. \n'grid_size' may be larger than the border extent.")
  }

  # --- The rest of the logic is the same as the previous function ---
  total_splits <- prod(grid_sizes) - 1
  max_regions <- total_splits + 1

  if (total_splits == 0) {
    return(create_tree(
      boundaries = numeric(0),
      split_dim_v = integer(0),
      split_label_v = integer(0),
      dim = dimension,
      border = border
    ))
  }

  splits_df <- data.frame(
    splits = numeric(total_splits),
    dim = integer(total_splits),
    split_part = integer(total_splits)
  )

  current_regions <- matrix(0.0, nrow = max_regions, ncol = 2 * dimension)
  current_regions[1, ] <- as.vector(border)
  n_regions <- 1
  split_idx <- 0

  lower_cols <- 1:dimension
  upper_cols <- (dimension + 1):(2 * dimension)

  for (d in 1:dimension) {
    n_subdivisions <- grid_sizes[d]
    if (n_subdivisions <= 1) {
      next
    }
    regions_to_process_indices <- 1:n_regions
    for (r_idx in regions_to_process_indices) {
      parent_region_bounds <- current_regions[r_idx, ]
      dim_min <- parent_region_bounds[d]
      dim_max <- parent_region_bounds[d + dimension]
      cell_width <- (dim_max - dim_min) / n_subdivisions
      region_to_split_label <- r_idx
      for (k in 1:(n_subdivisions - 1)) {
        split_idx <- split_idx + 1
        split_val <- dim_min + k * cell_width
        splits_df$splits[split_idx] <- split_val
        splits_df$dim[split_idx] <- d
        splits_df$split_part[split_idx] <- region_to_split_label
        sibling_region <- current_regions[region_to_split_label, ]
        current_regions[region_to_split_label, upper_cols[d]] <- split_val
        sibling_region[lower_cols[d]] <- split_val
        n_regions <- n_regions + 1
        current_regions[n_regions, ] <- sibling_region
        region_to_split_label <- n_regions
      }
    }
  }

  final_tree <- create_tree(
    boundaries = splits_df$splits,
    split_dim_v = splits_df$dim,
    split_label_v = splits_df$split_part,
    dim = dimension,
    border = border
  )

  final_tree$boundaries = treeBorders(final_tree)
  final_tree$boundaries = final_tree$boundaries[order(final_tree$boundaries[,2], final_tree$boundaries[,1]),]

  return(final_tree)
}

grid_tree = generate_grid_tree(0.01, c(0,0,1,1))

orderBoundaries(tree)[[2]]

tree = generateTree(999, border = c(0,0,10,10))
plotTreeGrid(tree)
plot(randomTree[[2]])

t2 = Sys.time()
t2-t1
plotTreeGrid(generateTree(100, border = c(0,0,10,10)))

Areas_EqualProb = apply(apply(treeBorders(randomTree), MARGIN = 1, FUN = diff, lag = randomTree$dimension), MARGIN = 2, FUN = prod)
Areas_WeightedProb = apply(apply(treeBorders(randomTree_wProb), MARGIN = 1, FUN = diff, lag = randomTree_wProb$dimension), MARGIN = 2, FUN = prod)
Areas_WeightedProb_wProb_Dim = apply(apply(treeBorders(randomTree_wProb_Dim), MARGIN = 1, FUN = diff, lag = randomTree_wProb_Dim$dimension), MARGIN = 2, FUN = prod)
Areas_WeightedProb_wProb_Dim_Beta = apply(apply(treeBorders(randomTree_wProb_Dim_Beta), MARGIN = 1, FUN = diff, lag = randomTree_wProb_Dim_Beta$dimension), MARGIN = 2, FUN = prod)

Areas = apply(apply(treeBorders(tree), MARGIN = 1, FUN = diff, lag = tree$dimension), MARGIN = 2, FUN = prod)

hist(log(Areas))
hist(log(Areas_WeightedProb))
boxplot(log(Areas_EqualProb), log(Areas_WeightedProb), log(Areas_WeightedProb_wProb_Dim), log(Areas_WeightedProb_wProb_Dim_Beta))
treeBorders(randomTree)

log(100/1000)

mean(sqrt(Areas))

sqrt(100/1000)


tree = create_tree(boundaries = c(0.25,0.75,0.5,0.25),
                   split_dim_v = c(1,2,1,2),
                   split_label_v = c(1,2,2,2) ,dim = 2,
                   border = c(0,0,2,5))

tree_modified = create_tree_modified(boundaries = c(0.25,0.75,0.5,0.25),
                                     split_dim_v = c(1,2,1,2),
                                     split_label_v = c(1,2,2,2) ,dim = 2,
                                     border = c(0,0,2,5))
newTree = createNewRegion(tree)


treeBorders(tree)
ArgoBart::treeBorders(tree)
plotTreeGrid(tree)

split_domain(matrix(c(0,0), nrow = 1), tree)


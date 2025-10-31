#Libraries
library(stringr)
library(ggplot2)
library(ggpubr)
library(MASS)
library(plotly)
library(MCMCpack)


#Creating Tree
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

generate_grid_tree = function(grid_size, border) {
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

  return(final_tree)
}

#Sampling Models

get_lowest_left_neighbor = function(treeBoundaries, index){

  all_left = which(treeBoundaries$U1 == treeBoundaries$L1[index] & (treeBoundaries$U2 >= treeBoundaries$L2[index] & treeBoundaries$L2 <= treeBoundaries$U2[index]) & treeBoundaries$order == 0)

  all_left[which.min(treeBoundaries$L2[all_left])]


}
get_lowest_down_neighbor = function(treeBoundaries, index){

  all_down = which(treeBoundaries$U2 == treeBoundaries$L2[index] & (treeBoundaries$U1 >= treeBoundaries$L1[index] & treeBoundaries$L1 <= treeBoundaries$U1[index]) & treeBoundaries$order == 0)

  all_down[which.min(treeBoundaries$L1[all_down])]


}
get_lowest_right_neighbor = function(treeBoundaries, index){

  all_right = which(treeBoundaries$L1 == treeBoundaries$U1[index] & (treeBoundaries$U2 >= treeBoundaries$L2[index] & treeBoundaries$L2 <= treeBoundaries$U2[index]) & treeBoundaries$order == 0)

  all_right[which.min(treeBoundaries$L2[all_right])]


}
get_lowest_up_neighbor = function(treeBoundaries, index){

  all_up = which(treeBoundaries$L2 == treeBoundaries$U2[index] & (treeBoundaries$U1 >= treeBoundaries$L1[index] & treeBoundaries$L1 <= treeBoundaries$U1[index]) & treeBoundaries$order == 0)

  all_up[which.min(treeBoundaries$L1[all_up])]


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

orderBoundaries_GeminiCleaned = function(tree, plot = T) {

  treeBoundaries <- treeBorders(tree)
  treeBoundaries$order <- 0

  # 1. Initialize by finding the bottom-left boundary and labeling it '1'
  which_BL <- which(treeBoundaries$L1 == tree$border[1] & treeBoundaries$L2 == tree$border[2])
  treeBoundaries$order[which_BL] <- 1

  # Initialize the 'heads' for the right- and up-seeking strategies
  latest_R <- which_BL
  latest_U <- which_BL

  # 2. Loop through and assign an order to each remaining boundary
  for (i in 2:nrow(treeBoundaries)) {

    found_R <- FALSE

    #--------------------------------------------------------------------#
    #  STRATEGY 1: Look for the next 'Right' candidate
    #--------------------------------------------------------------------#

    # Find an initial candidate. Priority is: right neighbor -> up neighbor -> then all the way left.
    which_next <- get_lowest_right_neighbor(treeBoundaries, index = latest_R)
    if (length(which_next) == 0) {
      which_next <- get_lowest_up_neighbor(treeBoundaries, index = latest_R)

      # INLINE LOGIC: Find the leftmost boundary in the row.
      # This loop runs as long as we can keep finding a neighbor to the left.
      while (TRUE) {
        if (length(which_next) == 0) break
        left_neighbor <- get_lowest_left_neighbor(treeBoundaries, index = which_next)
        if (length(left_neighbor) == 0) {
          break # Found the leftmost, exit.
        }
        which_next <- left_neighbor # Move left and continue.
      }
    }

    # Loop to test candidates until one is validated or we run out of options
    while (TRUE) {
      # First, guard against an empty index to prevent an error.
      if (length(which_next) == 0) break
      # Now, safely check if we've wrapped around to the leftmost wall.
      if (treeBoundaries$L1[which_next] == tree$border[1]) break

      # CORE VALIDATION LOGIC (RIGHT):
      # Check if the candidate's upper Y-coord is not higher than the highest
      # already-ordered boundary in the same column. This prevents 'jumping ahead'.
      max_U2_in_col <- max(treeBoundaries$U2[treeBoundaries$order != 0 & treeBoundaries$U1 == treeBoundaries$L1[which_next]])

      if (treeBoundaries$U2[which_next] <= max_U2_in_col) {
        # Valid candidate found. Label it, update the 'latest_R' head, and exit.
        treeBoundaries$order[which_next] <- i
        latest_R <- which_next
        found_R <- TRUE
        break
      } else {
        # Invalid candidate. Find the leftmost boundary in its row to try next.
        while (TRUE) {
          if (length(which_next) == 0) break
          left_neighbor <- get_lowest_left_neighbor(treeBoundaries, index = which_next)
          if (length(left_neighbor) == 0) {
            break # Found the leftmost, exit inner loop.
          }
          which_next <- left_neighbor # Move left and continue.
        }
      }
    }

    # If we successfully found a boundary with the 'Right' strategy, move to the next 'i'
    if (found_R) {
      next
    }

    #--------------------------------------------------------------------#
    #  STRATEGY 2: If Right failed, look for the next 'Up' candidate
    #--------------------------------------------------------------------#

    # Find a potential candidate. Priority: up neighbor -> right neighbor.
    which_next <- get_lowest_up_neighbor(treeBoundaries, index = latest_U)
    if (length(which_next) == 0) {
      which_next <- get_lowest_right_neighbor(treeBoundaries, index = latest_U)
    }

    # Loop to test candidates until one is validated or we run out
    while (TRUE) {
      if (length(which_next) == 0) break

      # CORE VALIDATION LOGIC (UP):
      # Check if the candidate's upper X-coord is not further right than the rightmost
      # already-ordered boundary in the same row.
      max_U1_in_row <- max(treeBoundaries$U1[treeBoundaries$order != 0 & treeBoundaries$U2 == treeBoundaries$L2[which_next]])

      if (treeBoundaries$U1[which_next] <= max_U1_in_row) {
        # Valid candidate found. Label it, update the 'latest_U' head, and exit.
        treeBoundaries$order[which_next] <- i
        latest_U <- which_next
        break
      } else {
        # Invalid candidate. Find the bottom-most boundary in its column to try next.
        while (TRUE) {
          if (length(which_next) == 0) break
          down_neighbor <- get_lowest_down_neighbor(treeBoundaries, index = which_next)
          if (length(down_neighbor) == 0) {
            break # Found the bottom-most, exit inner loop.
          }
          which_next <- down_neighbor # Move down and continue.
        }
      }
    }
  } # End of main for-loop

  # 3. Plotting (unchanged)

  if(plot){

    plot_noLabels <- plotTreeGrid(tree)
    plot_wLabels <- plot_noLabels + geom_text(
      data = treeBoundaries,
      aes(
        x = (L1 + U1) / 2,
        y = (L2 + U2) / 2,
        label = order
      ),
      color = "black",
      size = 2
    )

    list(treeBoundaries, plot_wLabels)

  } else{

    treeBoundaries

  }

}

find_neighbors = function(all_boundaries, index, sep = F, include_corners = F) {

  # Coordinates of the current boundary
  b <- all_boundaries[index, ]

  # Find boundaries that share an edge
  # Note: The small tolerance `1e-9` handles potential floating-point inaccuracies.
  if(include_corners){

    right_neighbors <- which(abs(all_boundaries$L1 - b$U1) <= 0 & all_boundaries$L2 <= b$U2 & all_boundaries$U2 >= b$L2)
    left_neighbors  <- which(abs(all_boundaries$U1 - b$L1) <= 0 & all_boundaries$L2 <= b$U2 & all_boundaries$U2 >= b$L2)
    up_neighbors    <- which(abs(all_boundaries$L2 - b$U2) <= 0 & all_boundaries$L1 <= b$U1 & all_boundaries$U1 >= b$L1)
    down_neighbors  <- which(abs(all_boundaries$U2 - b$L2) <= 0 & all_boundaries$L1 <= b$U1 & all_boundaries$U1 >= b$L1)

  } else{
    right_neighbors <- which(abs(all_boundaries$L1 - b$U1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
    left_neighbors  <- which(abs(all_boundaries$U1 - b$L1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
    up_neighbors    <- which(abs(all_boundaries$L2 - b$U2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)
    down_neighbors  <- which(abs(all_boundaries$U2 - b$L2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)
  }

  right_edges = apply(cbind(all_boundaries$L2[right_neighbors],
                            all_boundaries$U2[right_neighbors],
                            rep(b$L2, length(right_neighbors)),
                            rep(b$U2, length(right_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  left_edges = apply(cbind(all_boundaries$L2[left_neighbors],
                           all_boundaries$U2[left_neighbors],
                           rep(b$L2, length(left_neighbors)),
                           rep(b$U2, length(left_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  up_edges = apply(cbind(all_boundaries$L1[up_neighbors],
                         all_boundaries$U1[up_neighbors],
                         rep(b$L1, length(up_neighbors)),
                         rep(b$U1, length(up_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  down_edges = apply(cbind(all_boundaries$L1[down_neighbors],
                           all_boundaries$U1[down_neighbors],
                           rep(b$L1, length(down_neighbors)),
                           rep(b$U1, length(down_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})

  all_neighbors = c(right_neighbors, left_neighbors, up_neighbors, down_neighbors)
  all_edges = c(right_edges, left_edges, up_edges, down_edges)

  if(sep){

    return(list(Right = cbind(right_neighbors, right_edges),
                Left = cbind(left_neighbors, left_edges),
                Up = cbind(up_neighbors, up_edges),
                Down = cbind(down_neighbors, down_edges)))

  } else{
    return(cbind(all_neighbors, all_edges))
  }


}

sample_models_one_pass = function(tree, num_models, baseWeight = 0.1){

  treeBoundaries = data.frame(orderBoundaries_GeminiCleaned(tree)[[1]], model1 = 0, model2 = 0)

  baseTransMat = (baseWeight/num_models)*(matrix(rep(1, num_models^2), nrow = num_models) - diag(num_models))

  for(i in 1:nrow(treeBoundaries)){

    curTransMat = baseTransMat

    curIndex = which(treeBoundaries$order == i)
    neighbors_index_edge = find_neighbors(treeBoundaries, curIndex)

    definedNeighbors = which(treeBoundaries$order[neighbors_index_edge[,1]] < i)

    sum_of_defined_edges = sum(neighbors_index_edge[definedNeighbors, 2])

    if(sum_of_defined_edges > 0){

      for(j in definedNeighbors){

        curTransMat[treeBoundaries$model1[neighbors_index_edge[j,1]], treeBoundaries$model2[neighbors_index_edge[j,1]]] = curTransMat[treeBoundaries$model2[neighbors_index_edge[j,1]], treeBoundaries$model1[neighbors_index_edge[j,1]]] = curTransMat[treeBoundaries$model1[neighbors_index_edge[j,1]], treeBoundaries$model2[neighbors_index_edge[j,1]]] + neighbors_index_edge[j,2]/sum_of_defined_edges


      }

    }


    model1 = sample.int(num_models, size = 1, prob = rowSums(curTransMat))
    model2 = sample.int(num_models, size = 1, prob = curTransMat[model1,])

    treeBoundaries$model1[curIndex] = model1
    treeBoundaries$model2[curIndex] = model2


  }

  treeBoundaries

}

#Sampling Trajectories

evaluateCubicPatchValue = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(1, (y-y_0), (x-x_0), (x-x_0)*(y-y_0),
                (x-x_0)^2, (x-x_0)^2*(y-y_0), (x-x_0)^3, (x-x_0)^3*(y-y_0),
                (y-y_0)^2, (y-y_0)^3, (x-x_0)*(y-y_0)^2, (x-x_0)*(y-y_0)^3,
                (x-x_0)^2*(y-y_0)^2, (x-x_0)^2*(y-y_0)^3, (x-x_0)^3*(y-y_0)^2, (x-x_0)^3*(y-y_0)^3)

  sum(coef * polyTerms)

}
evaluateCubicPatchParX = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 0, 1, (y-y_0),
                2*(x-x_0), 2*(x-x_0)*(y-y_0), 3*(x-x_0)^2, 3*(x-x_0)^2*(y-y_0),
                0, 0, (y-y_0)^2, (y-y_0)^3,
                2*(x-x_0)*(y-y_0)^2, 2*(x-x_0)*(y-y_0)^3, 3*(x-x_0)^2*(y-y_0)^2, 3*(x-x_0)^2*(y-y_0)^3)

  sum(coef * polyTerms)
}
evaluateCubicPatchParY = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 1, 0, (x-x_0),
                0, (x-x_0)^2, 0, (x-x_0)^3,
                2*(y-y_0), 3*(y-y_0)^2, 2*(x-x_0)*(y-y_0), 3*(x-x_0)*(y-y_0)^2,
                2*(x-x_0)^2*(y-y_0), 3*(x-x_0)^2*(y-y_0)^2, 2*(x-x_0)^3*(y-y_0), 3*(x-x_0)^3*(y-y_0)^2)

  sum(coef * polyTerms)


}
evaluateCubicPatchParXY = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 0, 0, 1,
                0, 2*(x-x_0), 0, 3*(x-x_0)^2,
                0, 0, 2*(y-y_0), 3*(y-y_0)^2,
                4*(x-x_0)*(y-y_0), 6*(x-x_0)*(y-y_0)^2, 6*(x-x_0)^2*(y-y_0), 9*(x-x_0)^2*(y-y_0)^2)

  sum(coef * polyTerms)
}

calculatePatch_KnownDerivates = function(border, CornerValues, CornerParXs, CornerParYs, CornerParXYs){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  c00 = CornerValues[1]
  c01 = CornerParYs[1]
  c10 = CornerParXs[1]
  c11 = CornerParXYs[1]

  c20 = (-1*CornerParXs[2]*dx + 3*CornerValues[2] - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*CornerParXYs[2]*dx + 3*CornerParYs[2] - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (CornerParXs[2]*dx - 2*CornerValues[2] + 2*c00 + c10*dx)/(dx^3)
  c31 = (CornerParXYs[2]*dx - 2*CornerParYs[2] + 2*c01 + c11*dx)/(dx^3)

  c02 = (-1*CornerParYs[3]*dy + 3*CornerValues[3] - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*CornerParXYs[3]*dy + 3*CornerParXs[3] - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (CornerParYs[3]*dy - 2*CornerValues[3] + 2*c00 + c01*dy)/(dy^3)
  c13 = (CornerParXYs[3]*dy - 2*CornerParXs[3] + 2*c10 + c11*dy)/(dy^3)

  temp_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, 0, 0, 0, 0)

  Cz = evaluateCubicPatchValue(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cx = evaluateCubicPatchParX(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cy = evaluateCubicPatchParY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cxy = evaluateCubicPatchParXY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))

  c22 = ((CornerParXYs[4] - Cxy)*dx*dy - 3*(CornerParYs[4] - Cy)*dy - 3*(CornerParXs[4] - Cx)*dx + 9*(CornerValues[4] - Cz))/(dx^2*dy^2)
  c23 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 3*(CornerParYs[4] - Cy)*dy + 2*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx^2*dy^3)
  c32 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 2*(CornerParYs[4] - Cy)*dy + 3*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx^3*dy^2)
  c33 = ((CornerParXYs[4] - Cxy)*dx*dy - 2*(CornerParYs[4] - Cy)*dy - 2*(CornerParXs[4] - Cx)*dx + 4*(CornerValues[4] - Cz))/(dx^3*dy^3)


  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, c22, c23, c32, c33)

  known_coefs

}

calculateSurface_KnownCorners = function(boundaries, GridValues, GridParXs, GridParYs, GridParXYs){

  x_grid = sort(unique(c(boundaries$L1, boundaries$U1)))
  y_grid = sort(unique(c(boundaries$L2, boundaries$U2)))

  x_grid_len = length(x_grid)
  y_grid_len = length(y_grid)

  Coefs = matrix(nrow = nrow(boundaries), ncol = 16)

  for(i in 1:nrow(boundaries)){

    cur_border = as.numeric(boundaries[i,1:4])

    x_index = which(x_grid == cur_border[1])
    y_index = which(y_grid == cur_border[2])

    cur_corner_indices = c(0:1 + (y_index-1)*x_grid_len + x_index, 0:1 + (y_index)*x_grid_len + x_index)

    cur_CornerValues = GridValues[cur_corner_indices]
    cur_CornerParXs = GridParXs[cur_corner_indices]
    cur_CornerParYs = GridParYs[cur_corner_indices]
    cur_CornerParXYs = GridParXYs[cur_corner_indices]

    Coefs[i,] = calculatePatch_KnownDerivates(cur_border, cur_CornerValues, cur_CornerParXs, cur_CornerParYs, cur_CornerParXYs)

  }

  Coefs

}

tree = baseTree
sampledModels = base_models
model = 2
trans_prop


sampleTransitionSurface_ByGrid = function(tree, sampledModels, model, trans_prop, k, prior_mean = 1){

  #Sampling coefficients by sampling values and partial derivatives at the grid intersections

  baseBoundaries = tree$boundaries
  baseBoundaries$Type = "On"
  baseBoundaries$Type[sampledModels[,1] != model & sampledModels[,2] != model] = "Off"

  baseGridValues = c()
  baseGridParXs = c()
  baseGridParYs = c()
  baseGridParXYs = c()

  base_x_grid_points = sort(unique(c(baseBoundaries$L1, baseBoundaries$U1)))
  base_y_grid_points = sort(unique(c(baseBoundaries$L2, baseBoundaries$U2)))

  base_x_grid_len = length(base_x_grid_points)
  base_y_grid_len = length(base_y_grid_points)

  base_grid_points = expand.grid(base_x_grid_points, base_y_grid_points)

  for(i in 1:nrow(base_grid_points)){

    cur_BL_index = which(baseBoundaries$U1 == base_grid_points[i,1] & baseBoundaries$U2 == base_grid_points[i,2])
    cur_BR_index = which(baseBoundaries$L1 == base_grid_points[i,1] & baseBoundaries$U2 == base_grid_points[i,2])
    cur_TL_index = which(baseBoundaries$U1 == base_grid_points[i,1] & baseBoundaries$L2 == base_grid_points[i,2])
    cur_TR_index = which(baseBoundaries$L1 == base_grid_points[i,1] & baseBoundaries$L2 == base_grid_points[i,2])

    cur_neighbor_types = baseBoundaries$Type[c(cur_BL_index, cur_BR_index, cur_TL_index, cur_TR_index)]

    if(sum(cur_neighbor_types == "On") > 0){

      baseGridValues = c(baseGridValues, rnorm(1, mean = prior_mean, sd = k))
      baseGridParXs = c(baseGridParXs, rnorm(1, mean = 0, sd = k))
      baseGridParYs = c(baseGridParYs, rnorm(1, mean = 0, sd = k))
      baseGridParXYs = c(baseGridParXYs, rnorm(1, mean = 0, sd = k))

    } else{

      baseGridValues = c(baseGridValues, 0)
      baseGridParXs = c(baseGridParXs, 0)
      baseGridParYs = c(baseGridParYs, 0)
      baseGridParXYs = c(baseGridParXYs, 0)

    }

  }

  #Constructing all possible boundary regions

  if(missing(trans_prop)){
    trans_prop = rbeta(1, 10, 90)
  }

  x_grid_size = tree$boundaries[1,"U1"] - tree$boundaries[1,"L1"]
  y_grid_size = tree$boundaries[1,"U2"] - tree$boundaries[1,"L2"]

  trans_length = (x_grid_size + y_grid_size - sqrt((x_grid_size + y_grid_size)^2 - 4*trans_prop*x_grid_size*y_grid_size))/4

  newBoundaries = tree$boundaries
  newBoundaries$Type = "On"
  newBoundaries$Type[sampledModels[,1] != model & sampledModels[,2] != model] = "Off"
  newBoundaries$OGIndex = 1:nrow(newBoundaries)

  newSplits = tree$split

  x_cells_num = (tree$border[3] - tree$border[1])/x_grid_size
  y_cells_num = (tree$border[4] - tree$border[2])/x_grid_size

  for(i in 1:x_cells_num){

    curSplitLocLow = (i-1)*x_grid_size + trans_length
    curSplitLocHigh = i*x_grid_size - trans_length
    cellsToSplit = which(newBoundaries$L1 < curSplitLocLow & newBoundaries$U1 > curSplitLocHigh)
    cellTypes = newBoundaries[cellsToSplit,]$Type
    cellOGIndex = newBoundaries[cellsToSplit,]$OGIndex
    n_cells = length(cellsToSplit)

    curSplits = data.frame(splits = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), dim = 1, split_part = c(cellsToSplit, 1:n_cells + nrow(newBoundaries)))
    newSplits = rbind(newSplits, curSplits)

    newBoundaries[cellsToSplit,]$U1 = curSplitLocLow
    newBoundaries[cellsToSplit,]$Type = "T"

    curBoundaries = data.frame(L1 = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), L2 = rep(newBoundaries[cellsToSplit,]$L2,2), U1 = c(rep(curSplitLocHigh, n_cells), rep(curSplitLocHigh + trans_length, n_cells)), U2 = rep(newBoundaries[cellsToSplit,]$U2,2), Type = c(cellTypes, rep("T", n_cells)), OGIndex = rep(cellOGIndex, 2))

    newBoundaries = rbind(newBoundaries, curBoundaries)

  }


  for(i in 1:y_cells_num){

    curSplitLocLow = (i-1)*y_grid_size + trans_length
    curSplitLocHigh = i*y_grid_size - trans_length
    cellsToSplit = which(newBoundaries$L2 < curSplitLocLow & newBoundaries$U2 > curSplitLocHigh)
    cellTypes = newBoundaries[cellsToSplit,]$Type
    cellOGIndex = newBoundaries[cellsToSplit,]$OGIndex
    n_cells = length(cellsToSplit)

    curSplits = data.frame(splits = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), dim = 2, split_part = c(cellsToSplit, 1:n_cells + nrow(newBoundaries)))
    newSplits = rbind(newSplits, curSplits)

    newBoundaries[cellsToSplit,]$U2 = curSplitLocLow
    newBoundaries[cellsToSplit,]$Type = "T"

    curBoundaries = data.frame(L1 = rep(newBoundaries[cellsToSplit,]$L1,2), L2 = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), U1 = rep(newBoundaries[cellsToSplit,]$U1,2), U2 = c(rep(curSplitLocHigh, n_cells), rep(curSplitLocHigh + trans_length, n_cells)), Type = c(cellTypes, rep("T", n_cells)), OGIndex = rep(cellOGIndex, 2))

    newBoundaries = rbind(newBoundaries, curBoundaries)

  }

  #Defining true boundary regions

  transition_index = which(newBoundaries$Type == "T")

  neighbors_FirstPass = apply(matrix(transition_index), MARGIN = 1,
                              FUN = function(i){
                                curNeighbors = find_neighbors(all_boundaries = newBoundaries, index = i, sep = F, include_corners = T)
                                neighborTypes = newBoundaries$Type[curNeighbors]

                                if(sum(neighborTypes == "On") > 0){
                                  "On"
                                } else{
                                  "T"
                                }

                              })

  newBoundaries$Type[transition_index] = neighbors_FirstPass

  transition_index = which(newBoundaries$Type == "T")

  neighbors_SecondPass = apply(matrix(transition_index), MARGIN = 1,
                               FUN = function(i){
                                 curNeighbors = find_neighbors(all_boundaries = newBoundaries, index = i, sep = F, include_corners = T)
                                 neighborTypes = newBoundaries$Type[curNeighbors]

                                 if(sum(neighborTypes == "On") == 0){
                                   "Off"
                                 } else{
                                   "T"
                                 }

                               })
  newBoundaries$Type[transition_index] = neighbors_SecondPass

  #Updating Grid quantities

  GridValues = c()
  GridParXs = c()
  GridParYs = c()
  GridParXYs = c()

  x_grid_points = sort(unique(c(newBoundaries$L1, newBoundaries$U1)))
  y_grid_points = sort(unique(c(newBoundaries$L2, newBoundaries$U2)))

  all_grid_points = expand.grid(x_grid_points, y_grid_points)

  for(i in 1:nrow(all_grid_points)){

    cur_grid_point = as.numeric(all_grid_points[i,])

    cur_BL_index = which(newBoundaries$U1 == cur_grid_point[1] & newBoundaries$U2 == cur_grid_point[2])
    cur_BR_index = which(newBoundaries$L1 == cur_grid_point[1] & newBoundaries$U2 == cur_grid_point[2])
    cur_TL_index = which(newBoundaries$U1 == cur_grid_point[1] & newBoundaries$L2 == cur_grid_point[2])
    cur_TR_index = which(newBoundaries$L1 == cur_grid_point[1] & newBoundaries$L2 == cur_grid_point[2])

    cur_neighbor_types = newBoundaries$Type[c(cur_BL_index, cur_BR_index, cur_TL_index, cur_TR_index)]

    if(sum(cur_neighbor_types == "On") > 0){

      curOGIndex = which(baseBoundaries$L1 <= cur_grid_point[1] & baseBoundaries$U1 >= cur_grid_point[1] & baseBoundaries$L2 <= cur_grid_point[2] & baseBoundaries$U2 >= cur_grid_point[2])
      curOGIndex = curOGIndex[baseBoundaries$Type[curOGIndex] == "On"][1]
      cur_OGBorder = as.numeric(baseBoundaries[curOGIndex, 1:4])

      cur_OG_xIndex = which(base_x_grid_points == cur_OGBorder[1])
      cur_OG_yIndex = which(base_y_grid_points == cur_OGBorder[2])

      cur_OG_cornerIndices = c(0:1 + (cur_OG_yIndex-1)*base_x_grid_len + cur_OG_xIndex, 0:1 + (cur_OG_yIndex)*base_x_grid_len + cur_OG_xIndex)

      OGPatchCoefs = calculatePatch_KnownDerivates(border = as.numeric(baseBoundaries[curOGIndex, 1:4]),
                                                   CornerValues = baseGridValues[cur_OG_cornerIndices],
                                                   CornerParXs = baseGridParXs[cur_OG_cornerIndices],
                                                   CornerParYs = baseGridParYs[cur_OG_cornerIndices],
                                                   CornerParXYs = baseGridParXYs[cur_OG_cornerIndices])

      GridValues = c(GridValues, evaluateCubicPatchValue(coef = OGPatchCoefs, border = cur_OGBorder, curPos = cur_grid_point))
      GridParXs = c(GridParXs, evaluateCubicPatchParX(coef = OGPatchCoefs, border = cur_OGBorder, curPos = cur_grid_point))
      GridParYs = c(GridParYs, evaluateCubicPatchParY(coef = OGPatchCoefs, border = cur_OGBorder, curPos = cur_grid_point))
      GridParXYs = c(GridParXYs, evaluateCubicPatchParXY(coef = OGPatchCoefs, border = cur_OGBorder, curPos = cur_grid_point))

    } else{

      GridValues = c(GridValues, 0)
      GridParXs = c(GridParXs, 0)
      GridParYs = c(GridParYs, 0)
      GridParXYs = c(GridParXYs, 0)

    }

  }


  Coefs = calculateSurface_KnownCorners(newBoundaries[,1:4], GridValues, GridParXs, GridParYs, GridParXYs)

  TransitionTree = list(split = newSplits,
                        boundaries = newBoundaries,
                        dimension = tree$dimension,
                        border = tree$border,
                        coefs = Coefs)
  TransitionTree

}

sampleFullCubicSurface_AllModels = function(border_length = 1, cells_per_dim = 10, num_models = 3, k = 1, base_weight = 0.1, trans_prop = 8/9, prior_mean = 1){

  tree = generate_grid_tree(grid_size = border_length/cells_per_dim , c(0,0,border_length, border_length))

  sampledModels = sample_models_one_pass(tree = tree, num_models = num_models, baseWeight = base_weight)[,6:7]

  sampledCoefs = list()
  regionTypes = c()

  for(i in 1:num_models){

    cur_SampledTree = sampleTransitionSurface_ByGrid(tree, sampledModels, model = i, trans_prop, k, prior_mean)
    sampledCoefs[[i]] = cur_SampledTree$coefs
    regionTypes = cbind(regionTypes, cur_SampledTree$boundaries$Type)

  }

  finalTree = cur_SampledTree

  finalTree$boundaries = finalTree$boundaries[,1:4]
  finalTree$coefs = sampledCoefs
  finalTree$models = data.frame(regionTypes)
  names(finalTree$models) = str_c("Model",1:num_models)

  finalTree

}

#Sampling Particles

get_compSpace_pos = function(GMM, curPos_phy){

  n_models = length(GMM$Weights)

  x_cdf = sum(apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
    pnorm(curPos_phy[1], mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]))
  }) * GMM$Weights)

  cond_log_weights_uw = apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
    dnorm(curPos_phy[1], mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]), log = T) + log(GMM$Weights[model])
  })

  max_cond_log_weights_uw = max(cond_log_weights_uw)

  cond_log_weights = cond_log_weights_uw - (max_cond_log_weights_uw + log(sum(exp(cond_log_weights_uw - max_cond_log_weights_uw))))

  cond_weights = exp(cond_log_weights)

  y_given_x_cdf = sum(apply(matrix(1:n_models), MARGIN = 1, FUN = function(model){
    pnorm(curPos_phy[2], mean = GMM$Mean[model, 2] + GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*(curPos_phy[1] - GMM$Mean[model, 1]), sd = sqrt(GMM$Cov[[model]][2,2] - GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*GMM$Cov[[model]][1,2]))
  }) * cond_weights)

  c(x_cdf, y_given_x_cdf)

}

compSpaceData = function(GMM, phySpace_Data){

  compSpace = data.frame(t(apply(matrix(1:nrow(phySpace_Data)), MARGIN = 1, FUN = function(i){get_compSpace_pos(GMM, c(phySpace_Data$X[i], phySpace_Data$Y[i]))})), phySpace_Data$Model)
  names(compSpace) = c("X","Y","Model")

  compSpace

}

evaluateSampledTreeValue = function(tree, coefs, curPos){

  cur_index = which(curPos[1] >= tree$boundaries$L1 & curPos[1] <= tree$boundaries$U1 & curPos[2] >= tree$boundaries$L2 & curPos[2] <= tree$boundaries$U2)[1]

  cur_border = as.numeric(tree$boundaries[cur_index, 1:4])
  cur_coefs = as.numeric(coefs[cur_index,])

  evaluateCubicPatchValue(cur_coefs, cur_border, curPos)

}

TrajWeightedBaseVectorFields = function(t, curPos, baseVectorFields, compPatchTree, GMM){

  cur_CompPos = as.numeric(compSpaceData(GMM, data.frame(X = curPos[1],
                                                         Y = curPos[2],
                                                         Model = NA))[c(1,2)])

  cur_traj_value = as.numeric(lapply(X = compPatchTree$coefs, FUN = evaluateSampledTreeValue, tree = compPatchTree, curPos = cur_CompPos))

  cur_ModelVel = baseVectorFields(t, curPos)

  t(cur_ModelVel %*% matrix(cur_traj_value))

}

RungeKutta = function(startTime, startPos, baseVectorFields, compPatchTree, vel_sigma = 0.1, GMM, endTime, t_step = 0.01){

  n_dim = length(startPos)

  t_sim = seq(startTime, endTime, t_step)

  pos_sim = matrix(startPos, nrow = 1, byrow = T)

  for(i in 1:(length(t_sim)-1)){

    cur_t = t_sim[i]
    cur_pos = pos_sim[i,]

    k1 = TrajWeightedBaseVectorFields(cur_t, cur_pos, baseVectorFields, compPatchTree, GMM)
    k2 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k1, baseVectorFields, compPatchTree, GMM)
    k3 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k2, baseVectorFields, compPatchTree, GMM)
    k4 = TrajWeightedBaseVectorFields(cur_t + t_step, cur_pos + t_step*k3, baseVectorFields, compPatchTree, GMM)

    pos_sim = rbind(pos_sim, cur_pos + t_step*(k1+2*k2+2*k3+k4)/6 + rnorm(n_dim, mean = 0, sd = vel_sigma)*sqrt(t_step))

  }

  full_sim = data.frame(cbind(t_sim, pos_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))

  full_sim
}
EulerMaruyama = function(startTime, startPos, baseVectorFields, compPatchTree, vel_sigma = 0.1, GMM, endTime, t_step = 0.01){

  n_dim = length(startPos)

  t_sim = seq(startTime, endTime, t_step)

  pos_sim = matrix(startPos, nrow = 1, byrow = T)

  for(i in 1:(length(t_sim)-1)){

    cur_t = t_sim[i]
    cur_pos = pos_sim[i,]

    drift = TrajWeightedBaseVectorFields(cur_t, cur_pos, baseVectorFields, compPatchTree, GMM)
    diffusion = rnorm(n_dim, mean = 0, sd = vel_sigma)

    pos_sim = rbind(pos_sim, cur_pos + drift*t_step + diffusion*sqrt(t_step))

  }

  full_sim = data.frame(cbind(t_sim, pos_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))

  full_sim
}

samplePhySpaceParticles = function(n_particles, startTime, endTime, phySpaceBorder, phySpaceBorderBuffer = 0.1, baseVectorFields, compPatchTree, GMM, t_step = 0.01, vel_sigma = 0.1, pos_sigma = 0.1){

  x_phy_width = phySpaceBorder[3] - phySpaceBorder[1]
  y_phy_width = phySpaceBorder[4] - phySpaceBorder[2]

  startPos = data.frame(X = runif(n_particles, min = phySpaceBorder[1] + phySpaceBorderBuffer*x_phy_width, max = phySpaceBorder[3] - phySpaceBorderBuffer*x_phy_width),
                        Y = runif(n_particles, min = phySpaceBorder[2] + phySpaceBorderBuffer*y_phy_width, max = phySpaceBorder[4] - phySpaceBorderBuffer*y_phy_width))

  particleData_List = list()

  for(i in 1:n_particles){

    particleData_List[[i]] = cbind(RungeKutta(startTime = startTime, startPos = c(startPos$X[i], startPos$Y[i]), baseVectorFields = baseVectorFields, compPatchTree = compPatchTree, vel_sigma = vel_sigma, GMM = GMM, endTime = endTime, t_step = t_step), str_c("Particle",i))

    svMisc::progress(i, n_particles)
  }


  particleData_True = data.frame(do.call(rbind, particleData_List))

  names(particleData_True) = c('t', 'X1','X2', 'Particle')

  particleData_PosError = matrix(rnorm(2*nrow(particleData_True), mean = 0, sd = pos_sigma), ncol = 2)

  particleData_Obs = particleData_True
  particleData_Obs[,c(2:3)] = particleData_True[,c(2:3)] + particleData_PosError

  particleData_Obs

}

#Sampling GMM

sampleGMM = function(phy_space_border, lambda = 19, border_buffer_perc = 0.1, dim_sigma = 0.1){

  phy_x_width = phy_space_border[3] - phy_space_border[1]
  phy_y_width = phy_space_border[4] - phy_space_border[2]

  n_comp = rpois(n = 1,lambda = lambda) + 1

  mu_x = runif(n_comp, min = phy_space_border[1] + border_buffer_perc*phy_x_width, max = phy_space_border[3] - border_buffer_perc*phy_x_width)
  mu_y = runif(n_comp, min = phy_space_border[2] + border_buffer_perc*phy_y_width, max = phy_space_border[4] - border_buffer_perc*phy_y_width)

  mu = matrix(c(mu_x, mu_y), byrow = F, ncol = 2)

  Sigmas = replicate(n_comp, riwish(4, matrix(c((dim_sigma*phy_x_width)^2, 0,0,(dim_sigma*phy_y_width)^2), nrow = 2)), F)

  wts = as.vector(rdirichlet(n = 1, alpha = rep(1,n_comp)))


  GMM = list(Mean = mu, Cov = Sigmas, Weights = wts)

  GMM

}

#Sampling From Prior

baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

sampleFromFullPrior = function(phySpaceBorder = c(-1,-1,1,1), baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}, n_particles = 5, regions_per_dim_CompSpace = 10, traj_k_params = c(1,1),
                               base_weight_params = c(100,900), trans_proportion_params = c(80,10), traj_mean = 1,
                               error_k_params = c(1,1), n_GMM_mixtures_param = 19, GMM_buffer_perc_params = c(10,40), GMM_cov_params = c(10,40), startTime = 0, endTime = 5,
                               particle_t_step = 0.01, pos_error_params = c(1,1), vel_error_params = c(1,1)){

  #Generate Computational Space Trajectories

  n_models = ncol(baseVectorFields(t = 0, curPos = c(0.01,0.01))) #Fixed

  #regions_per_dim_CompSpace -> maybe add extra hyperprior

  traj_k = sqrt(rinvgamma(traj_k_params[1], shape = traj_k_params[2])) #std on the value and first derivatives of the trajectory surfaces at the corners of the grid
  error_k = sqrt(rinvgamma(error_k_params[1], shape = error_k_params[2])) #std on the value and first derivatives of the error surfaces at the corners of the grid (not sure about this one)
  base_weight_model_transition = rbeta(1, base_weight_params[1], base_weight_params[2]) #controls the stickiness of the models in the computational space
  trans_proportion = rbeta(1, trans_proportion_params[1], trans_proportion_params[2]) #the percent of each grid cell to be allocated to be a transition from On to Off

  cat("Sampling Vector Field Trajectory Surfaces...")

  modelTraj = sampleFullCubicSurface_AllModels(border_length = 1, cells_per_dim = regions_per_dim_CompSpace,
                                               num_models = n_models, k = traj_k, base_weight = base_weight_model_transition, #Trajectory Surfaces for the Vector Fields
                                               trans_prop = trans_proportion, prior_mean = traj_mean)

  #Generate Physical Space Transformation

  cat("Sampling Gaussian Mixture Model for Transition...")

  phySpaceBorder_Buffer = rbeta(1, GMM_buffer_perc_params[1], GMM_buffer_perc_params[2])/2 #the percent of the physical space border to treat as a buffer
  GMM_Marginal_Variance_Perc = rbeta(1, GMM_cov_params[1], GMM_cov_params[2])/2 #the std for the marginal variance of each mixture component (proportion of the physical space width)

  phySpaceGMM = sampleGMM(phy_space_border = phySpaceBorder, lambda = n_GMM_mixtures_param,
                          border_buffer_perc = phySpaceBorder_Buffer, dim_sigma = GMM_Marginal_Variance_Perc) #Sampled Gaussian Mixture Model to act as the non-uniform density in the physical space


  #Generate Particle Data

  cat("Sampling Particle Data...")

  PosErrorVar = rinvgamma(pos_error_params[1], pos_error_params[2]) #variance for the positional errors in the data (not sure about this either)
  VelErrorVar = rinvgamma(vel_error_params[1], vel_error_params[2]) #variance for the velocity stochastic errors in the data (not sure about this either)


  ParticleData = samplePhySpaceParticles(n_particles = n_particles, startTime = startTime, endTime = endTime,
                                         phySpaceBorder = phySpaceBorder, phySpaceBorderBuffer = phySpaceBorder_Buffer, #Samples true particle locations
                                         baseVectorFields = baseVectorFields,
                                         compPatchTree = modelTraj, GMM = phySpaceGMM, t_step = particle_t_step, vel_sigma = sqrt(VelErrorVar), pos_sigma = sqrt(PosErrorVar))
  ParticleData = na.omit(ParticleData)

  hyperparameters = c(traj_k = traj_k,
                      base_weight_model_transition = base_weight_model_transition,
                      transition_proportion = trans_proportion,
                      phy_space_border_buffer = phySpaceBorder_Buffer,
                      gmm_marginal_variance_proportion = GMM_Marginal_Variance_Perc,
                      pos_error_variance = PosErrorVar, pos_error_variance = VelErrorVar)

  list(ObsData = ParticleData, TrajectorySurfaces = modelTraj, TransformationGMM = phySpaceGMM, Hyperparameters = hyperparameters)

}





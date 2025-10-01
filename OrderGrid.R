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

orderBoundaries = function(tree){

  treeBoundaries = treeBorders(tree)
  treeBoundaries$order = rep(0,nrow(treeBoundaries))

  which_BL = which(treeBoundaries$L1 == tree$border[1] & treeBoundaries$L2 == tree$border[2])
  treeBoundaries$order[which_BL] = 1

  latest_R = which_BL
  latest_U = which_BL

  for(i in 2:nrow(treeBoundaries)){

    found_R = F
    found_U = F

    ## Look for next right candidate

    trying_R = T

    which_next = get_lowest_right_neighbor(treeBoundaries, index = latest_R)

    if(length(which_next) == 0){
      which_next = get_lowest_up_neighbor(treeBoundaries, index = latest_R)

      valid_next = length(get_lowest_left_neighbor(treeBoundaries, which_next)) == 0

      while(!valid_next){

        which_next = get_lowest_left_neighbor(treeBoundaries, index = which_next)

        valid_next = length(get_lowest_left_neighbor(treeBoundaries, which_next)) == 0

      }

    }

    if(treeBoundaries$L1[which_next] == tree$border[1]){
      trying_R = F
    }

    while(trying_R){

      if(treeBoundaries$U2[which_next] <= max(treeBoundaries$U2[treeBoundaries$order != 0 & treeBoundaries$U1 == treeBoundaries$L1[which_next]])){

        treeBoundaries$order[which_next] = i
        latest_R = which_next
        trying_R = F
        found_R = T
      } else{

        valid_next = F

        while(!valid_next){

          which_next = get_lowest_left_neighbor(treeBoundaries, index = which_next)

          valid_next = length(get_lowest_left_neighbor(treeBoundaries, which_next)) == 0

        }



      }

      if(treeBoundaries$L1[which_next] == tree$border[1]){

        trying_R = F

      }


    }

    if(found_R){
      next
    }

    #Trying Up


    trying_U = T

    which_next = get_lowest_up_neighbor(treeBoundaries, index = latest_U)

    if(length(which_next) == 0){
      which_next = get_lowest_right_neighbor(treeBoundaries, index = latest_U)
    }

    if(length(which_next) == 0){
      trying_U = F
    }


    while(trying_U){

      if(treeBoundaries$U1[which_next] <= max(treeBoundaries$U1[treeBoundaries$order != 0 & treeBoundaries$U2 == treeBoundaries$L2[which_next]])){

        treeBoundaries$order[which_next] = i
        latest_U = which_next
        trying_U = F
        found_U = T

      } else{

        valid_next = F

        while(!valid_next){

          which_next = get_lowest_down_neighbor(treeBoundaries, index = which_next)

          valid_next = length(get_lowest_down_neighbor(treeBoundaries, which_next)) == 0

        }

      }


    }

    if(!found_U){
      print(str_c(i, ' did not get a label'))
    }

  }

  plot_noLabels = plotTreeGrid(tree)


  plot_wLabels = plot_noLabels + geom_text(
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


}

orderBoundaries_GeminiCleaned <- function(tree) {

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
}


tree = generateTree(21, c(0,0,1,1))

t1 = Sys.time()
OB = orderBoundaries(tree)
t2 = Sys.time()
OB_G = orderBoundaries_GeminiCleaned(tree)
t3 = Sys.time()
OB_O = orderBoundaries_optimized(tree)
t4 = Sys.time()

OB[[2]]
OB_G[[2]]
OB_O[[2]]

t4-t3
t3-t2
t2-t1


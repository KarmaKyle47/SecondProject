// File: SurfaceTrajectories.stan

functions {

  /**
   * Helper function to map a 2D grid coordinate (i,j) to the
   * 1D row index in the finalBoundaries matrix.
   */
  int get_row_index(int i, int j, int comp_res) {
    int i_zero = i - 1;
    int j_zero = j - 1;
    int orig_row = i_zero / 3;
    int orig_col = j_zero / 3;
    int sub_row = i_zero % 3;
    int sub_col = j_zero % 3;
    int n_orig_index = orig_row * comp_res + orig_col + 1;
    int sub_cell_index = sub_row * 3 + sub_col + 1;
    return (n_orig_index - 1) * 9 + sub_cell_index;
  }

  /**
   * Your main function with the 3x3 split and 2-pass neighbor post-processing.
   */
  matrix updateCompGridTransitions(real trans_prop,
                                     matrix baseCompGridBoundaries,
                                     int comp_res,
                                     matrix models,
                                     int numModels) {

    // --- 1. Declarations ---
    int N = rows(baseCompGridBoundaries);
    int final_rows = N * 9;
    matrix[final_rows, 5 + numModels] finalBoundaries;

    int row_counter = 1;
    real x_grid_size = baseCompGridBoundaries[1, 3] - baseCompGridBoundaries[1, 1];
    real y_grid_size = baseCompGridBoundaries[1, 4] - baseCompGridBoundaries[1, 2];
    real trans_length = (x_grid_size + y_grid_size - sqrt((x_grid_size + y_grid_size)^2
                          - 4 * trans_prop * x_grid_size * y_grid_size)) / 4;

    vector[numModels] transition_models = rep_vector(0.5, numModels);
    vector[numModels] center_models;

    // --- 3. Main Loop ---
    for (n in 1:N) {
      real L1 = baseCompGridBoundaries[n, 1];
      real L2 = baseCompGridBoundaries[n, 2];
      real U1 = baseCompGridBoundaries[n, 3];
      real U2 = baseCompGridBoundaries[n, 4];
      real og_index = n;

      real x_split_1 = L1 + trans_length;
      real x_split_2 = U1 - trans_length;
      real y_split_1 = L2 + trans_length;
      real y_split_2 = U2 - trans_length;

      for (m in 1:numModels) {
        center_models[m] = (models[n, 1] == m || models[n, 2] == m);
      }

      // Bottom Row
      finalBoundaries[row_counter, 1:5] = [L1, L2, x_split_1, y_split_1, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_1, L2, x_split_2, y_split_1, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_2, L2, U1, y_split_1, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;

      // Middle Row
      finalBoundaries[row_counter, 1:5] = [L1, y_split_1, x_split_1, y_split_2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_1, y_split_1, x_split_2, y_split_2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = center_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_2, y_split_1, U1, y_split_2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;

      // Top Row
      finalBoundaries[row_counter, 1:5] = [L1, y_split_2, x_split_1, U2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_1, y_split_2, x_split_2, U2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
      finalBoundaries[row_counter, 1:5] = [x_split_2, y_split_2, U1, U2, og_index];
      finalBoundaries[row_counter, 6:(5+numModels)] = transition_models';
      row_counter += 1;
    }

    // --- 4. Post-Processing Neighbor Logic ---
    int grid_dim = comp_res * 3;
    for (m in 1:numModels) {
      vector[final_rows] current_model_col = finalBoundaries[, 5 + m];
      vector[final_rows] processed_model_col = current_model_col;

      // Pass 1: "T" -> "On" (1.0)
      for (i in 1:grid_dim) {
        for (j in 1:grid_dim) {
          int row_index = get_row_index(i, j, comp_res);
          if (current_model_col[row_index] == 0.5) {
            int has_on_neighbor = 0;
            for (ii in max(1, i - 1):min(grid_dim, i + 1)) {
              for (jj in max(1, j - 1):min(grid_dim, j + 1)) {
                if (ii == i && jj == j) continue;
                int neighbor_row_index = get_row_index(ii, jj, comp_res);
                if (current_model_col[neighbor_row_index] == 1.0) {
                  has_on_neighbor = 1;
                  break;
                }
              }
              if (has_on_neighbor == 1) break;
            }
            if (has_on_neighbor == 1) {
              processed_model_col[row_index] = 1.0;
            }
          }
        }
      }

      current_model_col = processed_model_col;

      // Pass 2: "T" -> "Off" (0.0)
      for (i in 1:grid_dim) {
        for (j in 1:grid_dim) {
          int row_index = get_row_index(i, j, comp_res);
          if (current_model_col[row_index] == 0.5) {
            int has_on_neighbor = 0;
            for (ii in max(1, i - 1):min(grid_dim, i + 1)) {
              for (jj in max(1, j - 1):min(grid_dim, j + 1)) {
                if (ii == i && jj == j) continue;
                int neighbor_row_index = get_row_index(ii, jj, comp_res);
                if (current_model_col[neighbor_row_index] == 1.0) {
                  has_on_neighbor = 1;
                  break;
                }
              }
              if (has_on_neighbor == 1) break;
            }
            if (has_on_neighbor == 0) {
              processed_model_col[row_index] = 0.0;
            }
          }
        }
      }

      finalBoundaries[, 5 + m] = processed_model_col;
    }

    return finalBoundaries;
  }

  real evaluateCubicPatchValue(vector[16] coefs, vector[4] border, vector[2] curPos) {
    // --- 1. Calculate relative coordinates ---
    real xr = curPos[1] - border[1]; // x - x_0
    real yr = curPos[2] - border[2]; // y - y_0

    // --- 2. Pre-calculate powers ---
    real xr2 = xr * xr;
    real xr3 = xr2 * xr;
    real yr2 = yr * yr;
    real yr3 = yr2 * yr;

    // --- 3. Build polynomial terms vector ---
    // (This must match the order of your 'coefs' vector)
    vector[16] polyTerms = [ 1.0, yr, xr, xr*yr,
                             xr2, xr2*yr, xr3, xr3*yr,
                             yr2, yr3, xr*yr2, xr*yr3,
                             xr2*yr2, xr2*yr3, xr3*yr2, xr3*yr3 ];

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParX(vector[16] coefs, vector[4] border, vector[2] curPos) {
    // --- 1. Calculate relative coordinates ---
    real xr = curPos[1] - border[1]; // x - x_0
    real yr = curPos[2] - border[2]; // y - y_0

    // --- 2. Pre-calculate powers ---
    real xr2 = xr * xr;
    real yr2 = yr * yr;
    real yr3 = yr2 * yr;

    // --- 3. Build polynomial terms vector ---
    // (This must match the order of your 'coefs' vector)
    vector[16] polyTerms = [ 0.0, 0.0, 1.0, yr,
                             2*xr, 2*xr*yr, 3*xr2, 3*xr2*yr,
                             0.0, 0.0, yr2, yr3,
                             2*xr*yr2, 2*xr*yr3, 3*xr2*yr2, 3*xr2*yr3 ];

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParY(vector[16] coefs, vector[4] border, vector[2] curPos) {
    // --- 1. Calculate relative coordinates ---
    real xr = curPos[1] - border[1]; // x - x_0
    real yr = curPos[2] - border[2]; // y - y_0

    // --- 2. Pre-calculate powers ---
    real xr2 = xr * xr;
    real xr3 = xr2 * xr;
    real yr2 = yr * yr;

    // --- 3. Build polynomial terms vector ---
    // (This must match the order of your 'coefs' vector)
    vector[16] polyTerms = [ 0.0, 1.0, 0.0, xr,
                             0.0, xr2, 0.0, xr3,
                             2*yr, 3*yr2, 2*xr*yr, 3*xr*yr2,
                             2*xr2*yr, 3*xr2*yr2, 2*xr3*yr, 3*xr3*yr2 ];

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParXY(vector[16] coefs, vector[4] border, vector[2] curPos) {
    // --- 1. Calculate relative coordinates ---
    real xr = curPos[1] - border[1]; // x - x_0
    real yr = curPos[2] - border[2]; // y - y_0

    // --- 2. Pre-calculate powers ---
    real xr2 = xr * xr;
    real yr2 = yr * yr;

    // --- 3. Build polynomial terms vector ---
    // (This must match the order of your 'coefs' vector)
    vector[16] polyTerms = [ 0.0, 0.0, 0.0, 1.0,
                             0.0, 2*xr, 0.0, 3*xr2,
                             0.0, 0.0, 2*yr, 3*yr2,
                             4*xr*yr, 6*xr*yr2, 6*xr2*yr, 9*xr2*yr2 ];

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateSampledSurfaceValue(matrix Boundaries, matrix coefs, vector[2] curPos) {

    // Use 0 as a sentinel value (since Stan indices are 1-based)
    int cur_index = 0;

    // --- 1. Find the correct patch index ---
    for (i in 1:rows(Boundaries)) {
      // Use logical AND '&&'
      if (curPos[1] >= Boundaries[i, 1] && curPos[1] <= Boundaries[i, 3] &&
          curPos[2] >= Boundaries[i, 2] && curPos[2] <= Boundaries[i, 4]) {

        cur_index = i; // Set the index
        break;         // Exit the loop
      }
    }

    // --- 2. Check if the point was found ---
    if (cur_index == 0) {
      // Point was not in any boundary. This is an error.
      // We reject the sample to stop computation.
      reject("Point (", curPos[1], ", ", curPos[2], ") was not found in any boundary patch.");
      return 0.0; // This line is unreachable but good practice
    }

    // --- 3. Extract data and evaluate ---
    // Transpose ' is correct to convert row_vector slice to vector
    vector[4] cur_border = Boundaries[cur_index, 1:4]';
    vector[16] cur_coefs = coefs[cur_index, 1:16]';

    return evaluateCubicPatchValue(cur_coefs, cur_border, curPos);
  }

  vector[16] calculatePatch_KnownDerivates(vector[4] border, vector[4] CornerValues, vector[4] CornerParXs, vector[4] CornerParYs, vector[4] CornerParXYs){

    real dx = border[3] - border[1];
    real dy = border[4] - border[2];

    real dx2 = dx * dx;
    real dx3 = dx2 * dx;
    real dy2 = dy * dy;
    real dy3 = dy2 * dy;

    real c00 = CornerValues[1];
    real c01 = CornerParYs[1];
    real c10 = CornerParXs[1];
    real c11 = CornerParXYs[1];

    real c20 = (-1*CornerParXs[2]*dx + 3*CornerValues[2] - 3*c00 - 2*c10*dx)/(dx2);
    real c21 = (-1*CornerParXYs[2]*dx + 3*CornerParYs[2] - 3*c01 - 2*c11*dx)/(dx2);
    real c30 = (CornerParXs[2]*dx - 2*CornerValues[2] + 2*c00 + c10*dx)/(dx3);
    real c31 = (CornerParXYs[2]*dx - 2*CornerParYs[2] + 2*c01 + c11*dx)/(dx3);

    real c02 = (-1*CornerParYs[3]*dy + 3*CornerValues[3] - 3*c00 - 2*c01*dy)/(dy2);
    real c12 = (-1*CornerParXYs[3]*dy + 3*CornerParXs[3] - 3*c10 - 2*c11*dy)/(dy2);
    real c03 = (CornerParYs[3]*dy - 2*CornerValues[3] + 2*c00 + c01*dy)/(dy3);
    real c13 = (CornerParXYs[3]*dy - 2*CornerParXs[3] + 2*c10 + c11*dy)/(dy3);

    vector[16] temp_coefs = [c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, 0.0, 0.0, 0.0, 0.0]';

    real Cz = evaluateCubicPatchValue(temp_coefs, border, border[3:4]);
    real Cx = evaluateCubicPatchParX(temp_coefs, border, border[3:4]);
    real Cy = evaluateCubicPatchParY(temp_coefs, border, border[3:4]);
    real Cxy = evaluateCubicPatchParXY(temp_coefs, border, border[3:4]);

    real c22 = ((CornerParXYs[4] - Cxy)*dx*dy - 3*(CornerParYs[4] - Cy)*dy - 3*(CornerParXs[4] - Cx)*dx + 9*(CornerValues[4] - Cz))/(dx2*dy2);
    real c23 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 3*(CornerParYs[4] - Cy)*dy + 2*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx2*dy3);
    real c32 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 2*(CornerParYs[4] - Cy)*dy + 3*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx3*dy2);
    real c33 = ((CornerParXYs[4] - Cxy)*dx*dy - 2*(CornerParYs[4] - Cy)*dy - 2*(CornerParXs[4] - Cx)*dx + 4*(CornerValues[4] - Cz))/(dx3*dy3);


    vector[16] known_coefs = [c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, c22, c23, c32, c33]';

    return known_coefs;

  }


  matrix calculateSurface_KnownCorners(matrix boundaries,
                                       vector GridValues,
                                       vector GridParXs,
                                       vector GridParYs,
                                       vector GridParXYs) {

    // --- 1. Declarations ---
    int n_cells = rows(boundaries); // This is N * 9

    // Calculate original grid dimensions
    int n_orig_cells = n_cells / 9; // Integer division
    int orig_grid_length = round(sqrt(n_orig_cells)); // N_orig_side

    // Calculate final (sub-divided) grid dimensions
    int grid_length = orig_grid_length * 3; // N_orig_side * 3

    // This is the number of CORNERS per side for the final grid
    int x_grid_len = grid_length + 1;

    // Pre-allocate the output matrix
    matrix[n_cells, 16] coefs;

    // --- 2. Main Loop ---
    // Loop through every cell
    for (i in 1:n_cells) {
      // --- A. Declarations for this loop ---
      int cur_corner_indices[4];
      vector[4] cur_border;
      vector[4] cur_CornerValues;
      vector[4] cur_CornerParXs;
      vector[4] cur_CornerParYs;
      vector[4] cur_CornerParXYs;

      // --- B. NEW LOGIC: Decode 'i' to find its grid position ---
      // (All 1-based indices)

      // 1. Find which original cell (1..N) and sub-cell (1..9) this is
      int n_orig_index = (i - 1) / 9 + 1;
      int sub_cell_index = (i - 1) % 9 + 1;

      // 2. Find (row, col) of the ORIGINAL cell
      int orig_cell_row = (n_orig_index - 1) / orig_grid_length + 1;
      int orig_cell_col = (n_orig_index - 1) % orig_grid_length + 1;

      // 3. Find (row, col) of the SUB-CELL (1..3)
      //    (Matches the 1-9 filling order: 1,2,3... 4,5,6... 7,8,9)
      int sub_cell_row = (sub_cell_index - 1) / 3 + 1;
      int sub_cell_col = (sub_cell_index - 1) % 3 + 1;

      // 4. Calculate the FINAL (row, col) of this cell in the 3x3 grid
      int final_cell_row = (orig_cell_row - 1) * 3 + sub_cell_row;
      int final_cell_col = (orig_cell_col - 1) * 3 + sub_cell_col;

      // 5. This (row, col) is the (y, x) index of the cell's
      //    bottom-left CORNER in the final (grid_length+1)x(grid_length+1) grid.
      int y_index = final_cell_row;
      int x_index = final_cell_col;

      // 6. Calculate the 4 flat indices for the GridXXX vectors
      cur_corner_indices[1] = (y_index - 1) * x_grid_len + x_index; // Bottom-Left (BL)
      cur_corner_indices[2] = cur_corner_indices[1] + 1;            // Bottom-Right (BR)
      cur_corner_indices[3] = y_index * x_grid_len + x_index;       // Top-Left (TL)
      cur_corner_indices[4] = cur_corner_indices[3] + 1;            // Top-Right (TR)

      // --- C. Get data for this cell ---
      // Get the cell's 4 boundaries (L1, L2, U1, U2)
      cur_border = boundaries[i, 1:4]'; // transpose row_vector -> vector

      // Get the 4 corner values/derivatives
      cur_CornerValues = GridValues[cur_corner_indices];
      cur_CornerParXs  = GridParXs[cur_corner_indices];
      cur_CornerParYs  = GridParYs[cur_corner_indices];
      cur_CornerParXYs = GridParXYs[cur_corner_indices];

      // --- D. Calculate coefficients ---
      // calculatePatch... returns a vector[16]. Transpose '
      // to assign it to the row_vector coefs[i, 1:16].
      coefs[i, 1:16] = calculatePatch_KnownDerivates(cur_border,
                                                    cur_CornerValues,
                                                    cur_CornerParXs,
                                                    cur_CornerParYs,
                                                    cur_CornerParXYs)';
    }

    return coefs;
  }

  matrix updateCornerQuantities(vector baseCornerValues,
                                vector baseCornerParXs,
                                vector baseCornerParYs,
                                vector baseCornerParXYs,
				                        matrix baseBoundaries,
                                matrix updatedBoundaries,
                                int model){

       // --- 1. Declarations ---
    int n_cells = rows(updatedBoundaries); // This is N * 9

    // Calculate original grid dimensions
    int n_orig_cells = n_cells / 9; // Integer division
    int orig_grid_length = round(sqrt(n_orig_cells)); // N_orig_side
    int orig_x_grid_len = orig_grid_length + 1;

    // Calculate final (sub-divided) grid dimensions
    int grid_length = orig_grid_length * 3; // N_orig_side * 3

    // This is the number of CORNERS per side for the final grid
    int x_grid_len = grid_length + 1;

    vector[n_cells] RegionTypes = updatedBoundaries[,5+model];

    int n_corners = (grid_length + 1) * (grid_length + 1);
    matrix[n_corners, 4] updatedCornerQuantities = rep_matrix(0.0, n_corners, 4);

    for(i in 1:n_cells){

      if(RegionTypes[i] == 1.0){

        // 1. Find which original cell (1..N) and sub-cell (1..9) this is
        int n_orig_index = (i - 1) / 9 + 1;
        int sub_cell_index = (i - 1) % 9 + 1;

        // 2. Find (row, col) of the ORIGINAL cell
        int orig_cell_row = (n_orig_index - 1) / orig_grid_length + 1;
        int orig_cell_col = (n_orig_index - 1) % orig_grid_length + 1;

	      int orig_corner_indices[4];

        // 6. Calculate the 4 flat indices for the GridXXX vectors
        orig_corner_indices[1] = (orig_cell_row - 1) * orig_x_grid_len + orig_cell_col; // Bottom-Left (BL)
        orig_corner_indices[2] = orig_corner_indices[1] + 1;            // Bottom-Right (BR)
        orig_corner_indices[3] = orig_cell_row * orig_x_grid_len + orig_cell_col;       // Top-Left (TL)
        orig_corner_indices[4] = orig_corner_indices[3] + 1;            // Top-Right (TR)

        // 3. Find (row, col) of the SUB-CELL (1..3)
        //    (Matches the 1-9 filling order: 1,2,3... 4,5,6... 7,8,9)
        int sub_cell_row = (sub_cell_index - 1) / 3 + 1;
        int sub_cell_col = (sub_cell_index - 1) % 3 + 1;

        // 4. Calculate the FINAL (row, col) of this cell in the 3x3 grid
        int final_cell_row = (orig_cell_row - 1) * 3 + sub_cell_row;
        int final_cell_col = (orig_cell_col - 1) * 3 + sub_cell_col;

        // 5. This (row, col) is the (y, x) index of the cell's
        //    bottom-left CORNER in the final (grid_length+1)x(grid_length+1) grid.
        int y_index = final_cell_row;
        int x_index = final_cell_col;

	      int cur_corner_indices[4];

        // 6. Calculate the 4 flat indices for the GridXXX vectors
        cur_corner_indices[1] = (y_index - 1) * x_grid_len + x_index; // Bottom-Left (BL)
        cur_corner_indices[2] = cur_corner_indices[1] + 1;            // Bottom-Right (BR)
        cur_corner_indices[3] = y_index * x_grid_len + x_index;       // Top-Left (TL)
        cur_corner_indices[4] = cur_corner_indices[3] + 1;            // Top-Right (TR)

        vector[4] cur_updatedBorder = updatedBoundaries[i,1:4]';
	      vector[4] cur_baseBorder = baseBoundaries[n_orig_index,1:4]';

      	vector[2] cur_BL = [cur_updatedBorder[1],cur_updatedBorder[2]]';
      	vector[2] cur_BR = [cur_updatedBorder[3],cur_updatedBorder[2]]';
	      vector[2] cur_TL = [cur_updatedBorder[1],cur_updatedBorder[4]]';
	      vector[2] cur_TR = [cur_updatedBorder[3],cur_updatedBorder[4]]';

        vector[16] cur_OrigCoefs = calculatePatch_KnownDerivates(cur_baseBorder, baseCornerValues[orig_corner_indices], 													 baseCornerParXs[orig_corner_indices],
										 baseCornerParYs[orig_corner_indices],
										 baseCornerParXYs[orig_corner_indices]);
        updatedCornerQuantities[cur_corner_indices[1],1:4] = [evaluateCubicPatchValue(cur_OrigCoefs, cur_baseBorder, cur_BL),
							      evaluateCubicPatchParX(cur_OrigCoefs, cur_baseBorder, cur_BL),
							      evaluateCubicPatchParY(cur_OrigCoefs, cur_baseBorder, cur_BL),
							      evaluateCubicPatchParXY(cur_OrigCoefs, cur_baseBorder, cur_BL)];

        updatedCornerQuantities[cur_corner_indices[2],1:4] = [evaluateCubicPatchValue(cur_OrigCoefs, cur_baseBorder, cur_BR),
							      evaluateCubicPatchParX(cur_OrigCoefs, cur_baseBorder, cur_BR),
							      evaluateCubicPatchParY(cur_OrigCoefs, cur_baseBorder, cur_BR),
							      evaluateCubicPatchParXY(cur_OrigCoefs, cur_baseBorder, cur_BR)];

        updatedCornerQuantities[cur_corner_indices[3],1:4] = [evaluateCubicPatchValue(cur_OrigCoefs, cur_baseBorder, cur_TL),
							      evaluateCubicPatchParX(cur_OrigCoefs, cur_baseBorder, cur_TL),
							      evaluateCubicPatchParY(cur_OrigCoefs, cur_baseBorder, cur_TL),
							      evaluateCubicPatchParXY(cur_OrigCoefs, cur_baseBorder, cur_TL)];

        updatedCornerQuantities[cur_corner_indices[4],1:4] = [evaluateCubicPatchValue(cur_OrigCoefs, cur_baseBorder, cur_TR),
							      evaluateCubicPatchParX(cur_OrigCoefs, cur_baseBorder, cur_TR),
							      evaluateCubicPatchParY(cur_OrigCoefs, cur_baseBorder, cur_TR),
							      evaluateCubicPatchParXY(cur_OrigCoefs, cur_baseBorder, cur_TR)];

      }


    }

    return updatedCornerQuantities;

  }

  vector[2] getCompSpacePos(matrix GMM_means, array[] cov_matrix[2,2] GMM_cov, simplex GMM_weights, vector[2] curPos_Phy){
    int n_mixtures = rows(GMM_means);
    real x_phy = curPos_Phy[1];
    real y_phy = curPos_Phy[2];

    real x_cdf = 0.0;
    vector[n_mixtures] cond_log_weights_uw;

    for(i in 1:n_mixtures){
      x_cdf += normal_cdf(x_phy, GMM_means[i,1], sqrt(GMM_cov[i][1,1]))*GMM_weights[i];
      cond_log_weights_uw[i] = normal_lpdf(x_phy | GMM_means[i,1], sqrt(GMM_cov[i][1,1])) + log(GMM_weights[i]);
    }

    vector[n_mixtures] cond_log_weights = cond_log_weights_uw - log_sum_exp(cond_log_weights_uw);
    vector[n_mixtures] cond_weights = exp(cond_log_weights);

    real y_given_x_cdf = 0.0;

    for(i in 1:n_mixtures){
      real cur_cond_mean = GMM_means[i,2] + GMM_cov[i][2,1]*(1/GMM_cov[i][1,1])*(x_phy - GMM_means[i,1]);
      real cur_cond_sd = sqrt(GMM_cov[i][2,2] - GMM_cov[i][2,1]*(1/GMM_cov[i][1,1])*GMM_cov[i][1,2]);

      y_given_x_cdf += normal_cdf(y_phy, cur_cond_mean, cur_cond_sd)*cond_weights[i];
    }

    return [x_cdf, y_given_x_cdf]';

  }

  real energySelf(vector[2] models){
    real cell_energy = 0.0;
    if(models[1] == models[2]){
      cell_energy += 1000000;
    }

    return cell_energy;
  }

  real energyPairs(vector[2] models1, vector[2] models2){

    real card_model1 = 2 - (models1[1] == models1[2]);
    real card_model2 = 2 - (models2[1] == models2[2]);

    real model_int = max(0.0, (models1[1] == models2[1] || models1[1] == models2[2]) +
                              (models1[2] == models2[1] || models1[2] == models2[2]) -
                              (models1[1] == models1[2]));

    real model_union = card_model1 + card_model2 - model_int;

    return 1 - (model_int)/(model_union);

  }

real calculateBaseGridEnergy(matrix baseCompGridBoundaries, matrix models) {

  int n_cells = rows(baseCompGridBoundaries);
  int n_grid_length = round(sqrt(n_cells));

  real total_energy = 0.0;

  // Single loop over all cells
  for (i in 1:n_cells) {

    // --- 1. Add Self Energy (same as before) ---
    total_energy += energySelf(models[i, 1:2]');

    // --- 2. Add "Right" Pair Energy ---
    // Check if 'i' is NOT in the far-right column
    // (The modulo operator % is 0 for the last cell in a row)
    if (i % n_grid_length != 0) {
      total_energy += energyPairs(models[i, 1:2]', models[i + 1, 1:2]');
    }

    // --- 3. Add "Up" Pair Energy ---
    // Check if 'i' is NOT in the top row
    if (i <= (n_cells - n_grid_length)) {
      total_energy += energyPairs(models[i, 1:2]', models[i + n_grid_length, 1:2]');
    }
  }

  return total_energy;
}



}

data {
  int<lower=0> N_data;      // Number of observations
  int<lower=1> comp_res;      // Computational Grid Resolution
  matrix[comp_res^2,4] baseCompGridBoundaries; //Base Computational Grid Boundaries
  matrix[N,2] Data;         // Particle Positions
  int<lower=1> GMM_num; // Number of Gaussian Mixtures in the Transformation
  matrix[GMM_num,2] GMM_means; // Means of all Gaussian Mixtures;
  array[GMM_num] cov_matrix[2,2] GMM_cov; // Covariance matrices of all Mixtures
  simplex[GMM_num] GMM_weights; // Weights for Mixtures
}

parameters {
  real alpha;          // Intercept
  real beta;           // Slope
  real<lower=0> sigma; // Residual standard deviation
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 5); // A weakly informative prior for scale parameters

  // Likelihood
  y ~ normal(alpha + beta * x, sigma);
}

generated quantities {
  vector[N] y_rep; // For posterior predictive checks
  for (n in 1:N) {
    y_rep[n] = normal_rng(alpha + beta * x[n], sigma);
  }
}

// File: SurfaceTrajectories.stan

functions {

  matrix generate_baseComp_boundaries(data int comp_res){

    int comp_res2 = comp_res*comp_res;

    vector[comp_res+1] full_dim_boundaries = linspaced_vector(comp_res + 1, 0.0, 1.0);

    vector[comp_res2] L1;
    vector[comp_res2] L2;
    vector[comp_res2] U1;
    vector[comp_res2] U2;

    for(i in 1:comp_res){

      int start_i = (i-1)*comp_res + 1;
      int end_i = (i)*comp_res;

      L1[start_i:end_i] = full_dim_boundaries[1:comp_res];
      L2[start_i:end_i] = rep_vector(full_dim_boundaries[i], comp_res);
      U1[start_i:end_i] = full_dim_boundaries[2:(comp_res + 1)];
      U2[start_i:end_i] = rep_vector(full_dim_boundaries[i+1], comp_res);

    }

    matrix[comp_res2, 4] Boundaries;

    Boundaries[,1] = L1;
    Boundaries[,2] = L2;
    Boundaries[,3] = U1;
    Boundaries[,4] = U2;

    return Boundaries;

  }

  matrix generate_TransComp_boundaries(int baseComp_res, real trans_prop){

    int full_res = 3*baseComp_res;
    int full_res2 = full_res*full_res;

    vector[baseComp_res+1] base_dim_boundaries = linspaced_vector(baseComp_res + 1, 0.0, 1.0);

    vector[full_res+1] full_dim_boundaries;

    real base_cell_size = 1.0/baseComp_res;
    real trans_length = base_cell_size*(1 - sqrt(1-trans_prop))/2;

    for(i in 1:baseComp_res){

      real cur_low = base_dim_boundaries[i];
      real cur_mid = base_dim_boundaries[i] + trans_length;
      real cur_high = base_dim_boundaries[i+1] - trans_length;

      int start_index = 3*(i-1);

      full_dim_boundaries[start_index + 1] = cur_low;
      full_dim_boundaries[start_index + 2] = cur_mid;
      full_dim_boundaries[start_index + 3] = cur_high;

    }

    full_dim_boundaries[full_res+1] = base_dim_boundaries[baseComp_res+1];

    vector[full_res2] L1;
    vector[full_res2] L2;
    vector[full_res2] U1;
    vector[full_res2] U2;

    for(i in 1:full_res){

      int start_i = (i-1)*full_res + 1;
      int end_i = (i)*full_res;

      L1[start_i:end_i] = full_dim_boundaries[1:full_res];
      L2[start_i:end_i] = rep_vector(full_dim_boundaries[i], full_res);
      U1[start_i:end_i] = full_dim_boundaries[2:(full_res + 1)];
      U2[start_i:end_i] = rep_vector(full_dim_boundaries[i+1], full_res);

    }

    matrix[full_res2, 4] Boundaries;

    Boundaries[,1] = L1;
    Boundaries[,2] = L2;
    Boundaries[,3] = U1;
    Boundaries[,4] = U2;

    return Boundaries;

  }

  real evaluateCubicPatchValue(vector coefs, vector border, vector curPos) {
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
                             xr2*yr2, xr2*yr3, xr3*yr2, xr3*yr3 ]';

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParX(vector coefs, vector border, vector curPos) {
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
                             2*xr*yr2, 2*xr*yr3, 3*xr2*yr2, 3*xr2*yr3 ]';

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParY(vector coefs, vector border, vector curPos) {
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
                             2*xr2*yr, 3*xr2*yr2, 2*xr3*yr, 3*xr3*yr2 ]';

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateCubicPatchParXY(vector coefs, vector border, vector curPos) {
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
                             4*xr*yr, 6*xr*yr2, 6*xr2*yr, 9*xr2*yr2 ]';

    // --- 4. Return the dot product ---
    return dot_product(coefs, polyTerms);
  }

  real evaluateSampledSurfaceValue(matrix Boundaries, data int comp_res, real[,] coefs, data vector curPos) {

    int cur_col_index = to_int(comp_res * curPos[1] + 1);
    int cur_row_index = to_int(comp_res * curPos[2] + 1);

    int cur_index = (cur_row_index-1)*comp_res + cur_col_index;


    vector[4] cur_border = Boundaries[cur_index, 1:4]';
    vector[16] cur_coefs = to_vector(coefs[cur_index, 1:16]);

    return evaluateCubicPatchValue(cur_coefs, cur_border, curPos);
  }

  vector calculatePatch_KnownDerivates(vector border, vector CornerValues, vector CornerParXs, vector CornerParYs, vector CornerParXYs){

    real dx = border[3] - border[1];
    real dy = border[4] - border[2];

    real inv_dx = 1.0 / dx;
    real inv_dx2 = inv_dx * inv_dx;
    real inv_dx3 = inv_dx2 * inv_dx;

    real inv_dy = 1.0 / dy;
    real inv_dy2 = inv_dy * inv_dy;
    real inv_dy3 = inv_dy2 * inv_dy;

    real c00 = CornerValues[1];
    real c01 = CornerParYs[1];
    real c10 = CornerParXs[1];
    real c11 = CornerParXYs[1];

    real c20 = (-1*CornerParXs[2]*dx + 3*CornerValues[2] - 3*c00 - 2*c10*dx)*(inv_dx2);
    real c21 = (-1*CornerParXYs[2]*dx + 3*CornerParYs[2] - 3*c01 - 2*c11*dx)*(inv_dx2);
    real c30 = (CornerParXs[2]*dx - 2*CornerValues[2] + 2*c00 + c10*dx)*(inv_dx3);
    real c31 = (CornerParXYs[2]*dx - 2*CornerParYs[2] + 2*c01 + c11*dx)*(inv_dx3);

    real c02 = (-1*CornerParYs[3]*dy + 3*CornerValues[3] - 3*c00 - 2*c01*dy)*(inv_dy2);
    real c12 = (-1*CornerParXYs[3]*dy + 3*CornerParXs[3] - 3*c10 - 2*c11*dy)*(inv_dy2);
    real c03 = (CornerParYs[3]*dy - 2*CornerValues[3] + 2*c00 + c01*dy)*(inv_dy3);
    real c13 = (CornerParXYs[3]*dy - 2*CornerParXs[3] + 2*c10 + c11*dy)*(inv_dy3);

    vector[16] temp_coefs = [c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, 0.0, 0.0, 0.0, 0.0]';

    real Cz = evaluateCubicPatchValue(temp_coefs, border, border[3:4]);
    real Cx = evaluateCubicPatchParX(temp_coefs, border, border[3:4]);
    real Cy = evaluateCubicPatchParY(temp_coefs, border, border[3:4]);
    real Cxy = evaluateCubicPatchParXY(temp_coefs, border, border[3:4]);

    real c22 = ((CornerParXYs[4] - Cxy)*dx*dy - 3*(CornerParYs[4] - Cy)*dy - 3*(CornerParXs[4] - Cx)*dx + 9*(CornerValues[4] - Cz))*(inv_dx2*inv_dy2);
    real c23 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 3*(CornerParYs[4] - Cy)*dy + 2*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))*(inv_dx2*inv_dy3);
    real c32 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 2*(CornerParYs[4] - Cy)*dy + 3*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))*(inv_dx3*inv_dy2);
    real c33 = ((CornerParXYs[4] - Cxy)*dx*dy - 2*(CornerParYs[4] - Cy)*dy - 2*(CornerParXs[4] - Cx)*dx + 4*(CornerValues[4] - Cz))*(inv_dx3*inv_dy3);


    vector[16] known_coefs = [c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, c22, c23, c32, c33]';

    return known_coefs;

  }

  matrix calculateSurface_KnownCorners(matrix boundaries,
                                       vector GridValues,
                                       vector GridParXs,
                                       vector GridParYs,
                                       vector GridParXYs) {

    // --- 1. Declarations ---
    int n_cells = rows(boundaries); // This is (N_orig_side * 3)^2

    // Calculate final (sub-divided) grid dimensions
    int grid_length = to_int(round(sqrt(n_cells))); // 3 * N_orig_side

    // This is the number of CORNERS per side for the final grid
    int x_grid_len = grid_length + 1; // 3 * N_orig_side + 1

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

      // --- B. NEW, SIMPLIFIED LOGIC ---
      // Decode 'i' (1-based) to find its 2D grid position.

      // 1. Calculate the FINAL (row, col) of this cell
      //    (Assumes row-major order)
      int final_cell_row = (i - 1) / grid_length + 1;
      int final_cell_col = (i - 1) % grid_length + 1;

      // 2. This (row, col) is the (y, x) index of the cell's
      //    bottom-left CORNER in the final (grid_length+1)x(grid_length+1) grid.
      int y_index = final_cell_row;
      int x_index = final_cell_col;

      // 3. Calculate the 4 flat indices for the GridXXX vectors
      cur_corner_indices[1] = (y_index - 1) * x_grid_len + x_index; // Bottom-Left (BL)
      cur_corner_indices[2] = cur_corner_indices[1] + 1;            // Bottom-Right (BR)
      cur_corner_indices[3] = y_index * x_grid_len + x_index;       // Top-Left (TL)
      cur_corner_indices[4] = cur_corner_indices[3] + 1;            // Top-Right (TR)

      // --- C. Get data for this cell (Unchanged) ---
      // Get the cell's 4 boundaries (L1, L2, U1, U2)
      cur_border = boundaries[i, 1:4]'; // transpose row_vector -> vector

      // Get the 4 corner values/derivatives
      cur_CornerValues = GridValues[cur_corner_indices];
      cur_CornerParXs  = GridParXs[cur_corner_indices];
      cur_CornerParYs  = GridParYs[cur_corner_indices];
      cur_CornerParXYs = GridParXYs[cur_corner_indices];

      // --- D. Calculate coefficients (Unchanged) ---
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

  matrix updateCornerQuantities(data int comp_res,
                                real trans_prop,
                                real[,,,] baseGridCornerQuantities, // [Corner, Qty, Cell, Model]
                                matrix baseGridBoundaries,
                                int model_num) {

    int baseGrid_res2 = comp_res * comp_res;
    int fine_grid_width = 3 * comp_res;
    int fine_point_width = fine_grid_width + 1;

    int TransGridRes = fine_grid_width * fine_grid_width; // (3*N)^2 cells
    int TransCornerRes = fine_point_width * fine_point_width; // (3*N + 1)^2 points

    // Call helper function (assuming this is correct)
    matrix[TransGridRes, 4] newBoundaries = generate_TransComp_boundaries(comp_res, trans_prop);

    // BUG FIX 1: Must use rep_vector to initialize zero-vectors
    vector[TransCornerRes] TransitionGridValue = rep_vector(0.0, TransCornerRes);
    vector[TransCornerRes] TransitionGridParX = rep_vector(0.0, TransCornerRes);
    vector[TransCornerRes] TransitionGridParY = rep_vector(0.0, TransCornerRes);
    vector[TransCornerRes] TransitionGridParXY = rep_vector(0.0, TransCornerRes);
    vector[TransCornerRes] TransitionGridCounts = rep_vector(0.0, TransCornerRes); // Use real for division

    for (i in 1:baseGrid_res2) {

      // This slicing looks correct based on your data structure
      vector[16] cur_coef = calculatePatch_KnownDerivates(
                            baseGridBoundaries[i, ]',
                            to_vector(baseGridCornerQuantities[, 1, i, model_num]),
                            to_vector(baseGridCornerQuantities[, 2, i, model_num]),
                            to_vector(baseGridCornerQuantities[, 3, i, model_num]),
                            to_vector(baseGridCornerQuantities[, 4, i, model_num]));

      // BUG FIX 2: Explicitly create R-style `vec + scalar` arrays
      int col_offset = 3 * ((i - 1) % comp_res);
      int row_offset = 3 * ((i - 1) / comp_res); // Stan does integer division

      int cur_col_indices[3] = {1 + col_offset, 2 + col_offset, 3 + col_offset};
      int cur_row_indices[3] = {1 + row_offset, 2 + row_offset, 3 + row_offset};

      // BUG FIX 3: Index arrays must be `int[]`, not `vector`.
      // And must be built explicitly, not with R's `vec + scalar` magic.
      int cur_grid_indices[9];
      cur_grid_indices[1] = (cur_row_indices[1] - 1) * fine_grid_width + cur_col_indices[1];
      cur_grid_indices[2] = (cur_row_indices[1] - 1) * fine_grid_width + cur_col_indices[2];
      cur_grid_indices[3] = (cur_row_indices[1] - 1) * fine_grid_width + cur_col_indices[3];
      cur_grid_indices[4] = (cur_row_indices[2] - 1) * fine_grid_width + cur_col_indices[1];
      cur_grid_indices[5] = (cur_row_indices[2] - 1) * fine_grid_width + cur_col_indices[2];
      cur_grid_indices[6] = (cur_row_indices[2] - 1) * fine_grid_width + cur_col_indices[3];
      cur_grid_indices[7] = (cur_row_indices[3] - 1) * fine_grid_width + cur_col_indices[1];
      cur_grid_indices[8] = (cur_row_indices[3] - 1) * fine_grid_width + cur_col_indices[2];
      cur_grid_indices[9] = (cur_row_indices[3] - 1) * fine_grid_width + cur_col_indices[3];

      matrix[16, 2] cur_pts;

      cur_pts[1,]  = newBoundaries[cur_grid_indices[1], 1:2];
      cur_pts[5,]  = newBoundaries[cur_grid_indices[4], 1:2];
      cur_pts[9,]  = newBoundaries[cur_grid_indices[7], 1:2];

      // Slices that were [3, 2] -> [U1, L2] (Bottom-Right)
      cur_pts[2,]  = [newBoundaries[cur_grid_indices[1], 3], newBoundaries[cur_grid_indices[1], 2]];
      cur_pts[3,]  = [newBoundaries[cur_grid_indices[2], 3], newBoundaries[cur_grid_indices[2], 2]];
      cur_pts[4,]  = [newBoundaries[cur_grid_indices[3], 3], newBoundaries[cur_grid_indices[3], 2]];
      cur_pts[6,]  = [newBoundaries[cur_grid_indices[4], 3], newBoundaries[cur_grid_indices[4], 2]];
      cur_pts[7,]  = [newBoundaries[cur_grid_indices[5], 3], newBoundaries[cur_grid_indices[5], 2]];
      cur_pts[8,]  = [newBoundaries[cur_grid_indices[6], 3], newBoundaries[cur_grid_indices[6], 2]];
      cur_pts[10,] = [newBoundaries[cur_grid_indices[7], 3], newBoundaries[cur_grid_indices[7], 2]];
      cur_pts[11,] = [newBoundaries[cur_grid_indices[8], 3], newBoundaries[cur_grid_indices[8], 2]];
      cur_pts[12,] = [newBoundaries[cur_grid_indices[9], 3], newBoundaries[cur_grid_indices[9], 2]];

      // Slices that were [1, 4] -> [L1, U2] (Top-Left)
      cur_pts[13,] = [newBoundaries[cur_grid_indices[7], 1], newBoundaries[cur_grid_indices[7], 4]];

      // Slices that were [3, 4] -> [U1, U2] (Top-Right)
      cur_pts[14,] = newBoundaries[cur_grid_indices[7], 3:4];
      cur_pts[15,] = newBoundaries[cur_grid_indices[8], 3:4];
      // (And fixing the nwwBoundaries typo)
      cur_pts[16,] = newBoundaries[cur_grid_indices[9], 3:4];

      // Note: I simplified the R code for cur_pts.
      // The R code had 16 unique-looking calls.
      // Your Stan code has many duplicates (e.g., [7],[3,4] is used 3 times).
      // I have preserved YOUR Stan logic, but you should double-check
      // it against the R code's `matrix(c(...))` call, which was harder to read.
      // The R code seemed to be picking the 16 *corners* of the 3x3 cell grid.

      // BUG FIX 2 (again): Explicit array creation
      int col_offset_pts = 3 * ((i - 1) % comp_res);
      int row_offset_pts = 3 * ((i - 1) / comp_res);

      int cur_col_indices_pts[4] = {1 + col_offset_pts, 2 + col_offset_pts, 3 + col_offset_pts, 4 + col_offset_pts};
      int cur_row_indices_pts[4] = {1 + row_offset_pts, 2 + row_offset_pts, 3 + row_offset_pts, 4 + row_offset_pts};

      // BUG FIX 3 (again): Must be int[]
      int cur_pts_indices[16];
      cur_pts_indices[1]  = (cur_row_indices_pts[1] - 1) * fine_point_width + cur_col_indices_pts[1];
      cur_pts_indices[2]  = (cur_row_indices_pts[1] - 1) * fine_point_width + cur_col_indices_pts[2];
      cur_pts_indices[3]  = (cur_row_indices_pts[1] - 1) * fine_point_width + cur_col_indices_pts[3];
      cur_pts_indices[4]  = (cur_row_indices_pts[1] - 1) * fine_point_width + cur_col_indices_pts[4];
      cur_pts_indices[5]  = (cur_row_indices_pts[2] - 1) * fine_point_width + cur_col_indices_pts[1];
      cur_pts_indices[6]  = (cur_row_indices_pts[2] - 1) * fine_point_width + cur_col_indices_pts[2];
      cur_pts_indices[7]  = (cur_row_indices_pts[2] - 1) * fine_point_width + cur_col_indices_pts[3];
      cur_pts_indices[8]  = (cur_row_indices_pts[2] - 1) * fine_point_width + cur_col_indices_pts[4];
      cur_pts_indices[9]  = (cur_row_indices_pts[3] - 1) * fine_point_width + cur_col_indices_pts[1];
      cur_pts_indices[10] = (cur_row_indices_pts[3] - 1) * fine_point_width + cur_col_indices_pts[2];
      cur_pts_indices[11] = (cur_row_indices_pts[3] - 1) * fine_point_width + cur_col_indices_pts[3];
      cur_pts_indices[12] = (cur_row_indices_pts[3] - 1) * fine_point_width + cur_col_indices_pts[4];
      cur_pts_indices[13] = (cur_row_indices_pts[4] - 1) * fine_point_width + cur_col_indices_pts[1];
      cur_pts_indices[14] = (cur_row_indices_pts[4] - 1) * fine_point_width + cur_col_indices_pts[2];
      cur_pts_indices[15] = (cur_row_indices_pts[4] - 1) * fine_point_width + cur_col_indices_pts[3];
      cur_pts_indices[16] = (cur_row_indices_pts[4] - 1) * fine_point_width + cur_col_indices_pts[4];

      // BUG FIX 1 (again): Initialization
      vector[16] curValues = rep_vector(0.0, 16);
      vector[16] curParXs = rep_vector(0.0, 16);
      vector[16] curParYs = rep_vector(0.0, 16);
      vector[16] curParXYs = rep_vector(0.0, 16);

      for (j in 1:16) {
        // Assuming these helper functions are correct
        curValues[j] = evaluateCubicPatchValue(cur_coef, baseGridBoundaries[i, ]', cur_pts[j, ]');
        curParXs[j]  = evaluateCubicPatchParX(cur_coef, baseGridBoundaries[i, ]', cur_pts[j, ]');
        curParYs[j]  = evaluateCubicPatchParY(cur_coef, baseGridBoundaries[i, ]', cur_pts[j, ]');
        curParXYs[j] = evaluateCubicPatchParXY(cur_coef, baseGridBoundaries[i, ]', cur_pts[j, ]');
      }

      // BUG FIX 9 (THE BIG ONE): Stan cannot "scatter-add"
      // `TransitionGridValue[cur_pts_indices] = ...` is invalid.
      // You MUST use an explicit loop.
      for (k in 1:16) {
        int idx = cur_pts_indices[k]; // Get the point's 1D index

        TransitionGridValue[idx] = TransitionGridValue[idx] + curValues[k];
        TransitionGridParX[idx]  = TransitionGridParX[idx] + curParXs[k];
        TransitionGridParY[idx]  = TransitionGridParY[idx] + curParYs[k];
        TransitionGridParXY[idx] = TransitionGridParXY[idx] + curParXYs[k];
        TransitionGridCounts[idx] = TransitionGridCounts[idx] + 1;
      }
    }

    // BUG FIX 10: Use `./` for element-wise division
    TransitionGridValue = TransitionGridValue ./ TransitionGridCounts;
    TransitionGridParX = TransitionGridParX ./ TransitionGridCounts;
    TransitionGridParY = TransitionGridParY ./ TransitionGridCounts;
    TransitionGridParXY = TransitionGridParXY ./ TransitionGridCounts;

    matrix[TransCornerRes, 4] curModelTransitionQuantities;

    curModelTransitionQuantities[, 1] = TransitionGridValue;
    curModelTransitionQuantities[, 2] = TransitionGridParX;
    curModelTransitionQuantities[, 3] = TransitionGridParY;
    curModelTransitionQuantities[, 4] = TransitionGridParXY;

    return curModelTransitionQuantities;
  }


  vector getCompSpacePos(data matrix GMM_means, data real[,,] GMM_cov, data vector GMM_weights, data vector curPos_Phy,
                         data vector GMM_sd_x, data vector GMM_sd_y_cond, data vector GMM_cond_slope){
    int n_mixtures = rows(GMM_means);
    real x_phy = curPos_Phy[1];
    real y_phy = curPos_Phy[2];

    real x_cdf = 0.0;
    vector[n_mixtures] cond_log_weights_uw;

    for(i in 1:n_mixtures){
      x_cdf += normal_cdf(x_phy, GMM_means[i,1], GMM_sd_x[i])*GMM_weights[i];
      cond_log_weights_uw[i] = normal_lpdf(x_phy | GMM_means[i,1], GMM_sd_x[i]) + log(GMM_weights[i]);
    }

    vector[n_mixtures] cond_log_weights = cond_log_weights_uw - log_sum_exp(cond_log_weights_uw);
    vector[n_mixtures] cond_weights = exp(cond_log_weights);

    real y_given_x_cdf = 0.0;

    for(i in 1:n_mixtures){
      real cur_cond_mean = GMM_means[i,2] + GMM_cond_slope[i]*(x_phy - GMM_means[i,1]);
      real cur_cond_sd = GMM_sd_y_cond[i];

      y_given_x_cdf += normal_cdf(y_phy, cur_cond_mean, cur_cond_sd)*cond_weights[i];
    }

    return [x_cdf, y_given_x_cdf]';

  }

  real energySelf(row_vector cell_probs, data real penalty, data int n_models){

    // row_vector[n_models] cell_probs;
    // cell_probs = inv_logit(cell_logits);

    real cell_energy = penalty*(sum(cell_probs) - 2)^2;
    return cell_energy;

  }

  real calculateBaseGridEnergy(matrix baseCompGridBoundaries, matrix ModelProbs, data real selfPenalty, data int n_models) {

    int n_cells = rows(baseCompGridBoundaries);
    int n_grid_length = to_int(round(sqrt(n_cells)));

    real total_energy = 0.0;

    // Single loop over all cells
    for (i in 1:n_cells) {

      // --- 1. Add Self Energy (same as before) ---
      total_energy += energySelf(ModelProbs[i, 1:n_models], selfPenalty, n_models);

      // --- 2. Add "Right" Pair Energy ---
      // Check if 'i' is NOT in the far-right column
      // (The modulo operator % is 0 for the last cell in a row)
      if (i % n_grid_length != 0) {
        total_energy += mean(square(ModelProbs[i, 1:n_models] - ModelProbs[i + 1, 1:n_models]));
      }

      // --- 3. Add "Up" Pair Energy ---
      // Check if 'i' is NOT in the top row
      if (i <= (n_cells - n_grid_length)) {
        total_energy += mean(square(ModelProbs[i, 1:n_models] - ModelProbs[i + n_grid_length, 1:n_models]));
      }
    }

    return total_energy;
  }

  real modelLogits_lpdf(matrix ModelProbs, matrix baseCompGridBoundaries, real logit_temp, data real selfPenalty, data int n_models){

    real energy;

    energy = calculateBaseGridEnergy(baseCompGridBoundaries, ModelProbs, selfPenalty, n_models);

    return -1*energy/logit_temp;


  }

  matrix baseVectorFields(data real t, data vector curPos){

    matrix[2,2] VF;

    real inv_norm = inv_sqrt(curPos[1] * curPos[1] + curPos[2] * curPos[2]);

    VF[1,1] = curPos[2]*inv_norm;
    VF[2,1] = -1*curPos[1]*inv_norm;
    VF[1,2] = curPos[1]*inv_norm;
    VF[2,2] = curPos[2]*inv_norm;

    return VF;

  }

  vector TrajWeightedBaseVectorFields(data real t, data vector phySpacePos, data vector compSpacePos,
                                      matrix Boundaries, data int trans_res, real[,,] coefs, data int N_models){


    vector[N_models] TrajValues;

    for(i in 1:N_models){

      TrajValues[i] = evaluateSampledSurfaceValue(Boundaries, trans_res, coefs[,,i], compSpacePos);

    }

    matrix[2,N_models] ModelVels = baseVectorFields(t, phySpacePos);

    vector[2] TrajWeightedVel = ModelVels * TrajValues;

    return TrajWeightedVel;

  }


}

data {
  int<lower=0> N_data;      // Number of observations
  int<lower=1> comp_res;      // Computational Grid Resolution
  matrix[N_data,5] Data;         // Particle Velocities with Positions for now (t, x, y, v_x, v_y)
  int<lower=1> GMM_num; // Number of Gaussian Mixtures in the Transformation
  matrix[GMM_num,2] GMM_means; // Means of all Gaussian Mixtures;
  real GMM_cov[2,2,GMM_num]; // Covariance matrices of all Mixtures
  simplex[GMM_num] GMM_weights; // Weights for Mixtures
  real trans_prop;
  real logit_temp;
}

transformed data {
  real<lower=0> traj_k_alpha;      // traj_k alpha for inv_gamma
  real<lower=0> traj_k_beta;       // traj_k beta for inv_gamma
  //real<lower=0> logit_temp_alpha;  // logit_temp alpha for inv_gamma
  //real<lower=0> logit_temp_beta;   // logit_temp beta for inv_gamma
  //real<lower=0> trans_prop_alpha;  // trans_prop alpha for beta
  //real<lower=0> trans_prop_beta;   // trans_prop beta for beta
  real<lower=0> sigma_vel_alpha;  // logit_temp alpha for inv_gamma
  real<lower=0> sigma_vel_beta;   // logit_temp beta for inv_gamma
  traj_k_alpha = 0.1;
  traj_k_beta = 0.1;
  //logit_temp_alpha = 1;
  //logit_temp_beta = 24;
  sigma_vel_alpha = 0.1;
  sigma_vel_beta = 0.1;
  //trans_prop_alpha = 1;
  //trans_prop_beta = 1;

  real<lower=0> selfPenalty; // self penalty for the modelLogits energy (sum = 2)
  selfPenalty = 10;

  real prior_traj_mean; // prior mean for the value of trajectories which are "on"
  prior_traj_mean = 1;

  real<lower=0> off_traj_sd; // sd for traj quantities which are "off"
  off_traj_sd = 0.01;

  int trans_res = 3*comp_res;
  int comp_res2 = comp_res*comp_res;
  int trans_res2 = trans_res*trans_res;
  int trans_corner_res = trans_res+1;
  int trans_corner_res2 = trans_corner_res*trans_corner_res;

  matrix[comp_res2,4] baseBoundaries = generate_baseComp_boundaries(comp_res);

  int N_models = 2;    // Number of models

  // --- Pre-calculate GMM values ---
  vector[GMM_num] GMM_sd_x;        // sqrt(cov[1,1])
  vector[GMM_num] GMM_sd_y_cond;   // Sqrt of conditional variance
  vector[GMM_num] GMM_cond_slope;  // cov[2,1] / cov[1,1]

  for(i in 1:GMM_num) {
    real inv_cov_11 = 1.0 / GMM_cov[1,1,i];
    real cond_slope = GMM_cov[2,1,i] * inv_cov_11;
    GMM_sd_x[i] = sqrt(GMM_cov[1,1,i]);
    GMM_cond_slope[i] = cond_slope;
    GMM_sd_y_cond[i] = sqrt(GMM_cov[2,2,i] - cond_slope * GMM_cov[1,2,i]);
  }

  matrix[N_data, 2] compSpacePos_data;
  for (i in 1:N_data) {
    compSpacePos_data[i, 1:2] = getCompSpacePos(GMM_means, GMM_cov, GMM_weights, Data[i, 2:3]',
                                                GMM_sd_x, GMM_sd_y_cond, GMM_cond_slope)';
  }

  matrix[trans_res2, 4] cur_transBoundaries = generate_TransComp_boundaries(comp_res, trans_prop);

}

parameters {
  real<lower=0> traj_k;      // SD of the trajectory corner quantities, if on
  //real<lower=0> logit_temp;    // Temperature for the boltzmann distribution on the model logits
  //real<lower=0,upper=1> trans_prop; // Proportion of base cell to act as transition
  matrix<lower=0,upper=1>[comp_res*comp_res, N_models] modelProbs; // The logit (equivelantly probability) of each model being on in each base region
  real baseCornerQuanities[4,4,comp_res*comp_res,N_models]; // Corner quanities for the base grid
  real<lower=0> sigma_vel; // Velocity Sigma
}

transformed parameters {

  real UpdatedCornerQuantities[trans_corner_res2, 4, N_models];

  for(i in 1:N_models){

    UpdatedCornerQuantities[1:trans_corner_res2, 1:4, i] = to_array_2d(updateCornerQuantities(comp_res, trans_prop, baseCornerQuanities, baseBoundaries, i));

  }

  real cur_Coefs[trans_res2, 16, N_models];

  for(i in 1:N_models){

    cur_Coefs[1:trans_res2, 1:16, i] = to_array_2d(calculateSurface_KnownCorners(cur_transBoundaries,
                                                                    to_vector(UpdatedCornerQuantities[,1,i]),
                                                                    to_vector(UpdatedCornerQuantities[,2,i]),
                                                                    to_vector(UpdatedCornerQuantities[,3,i]),
                                                                    to_vector(UpdatedCornerQuantities[,4,i])));
  }

  // matrix[comp_res2, N_models] log_prob_on;
  // matrix[comp_res2, N_models] log_prob_off;
  //
  // for (m in 1:N_models) {
  //   // modelLogits[, m] is a 'vector'
  //   log_prob_on[, m] = log_inv_logit(modelLogits[, m]);
  //   log_prob_off[, m] = log1m_inv_logit(modelLogits[, m]);
  // }

}

model {
  // Priors
  traj_k ~ inv_gamma(traj_k_alpha, traj_k_beta); // approx jeffrey's
  //logit_temp ~ inv_gamma(logit_temp_alpha, logit_temp_beta); //approx jeffrey's
  //trans_prop ~ beta(trans_prop_alpha, trans_prop_beta); //uniform
  sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta); //approx jeffrey's

  target += modelLogits_lpdf(modelProbs | baseBoundaries, logit_temp, selfPenalty, N_models);

  for(m in 1:N_models){

    for(c in 1:comp_res2){

      real cur_log_prob_on = log(modelProbs[c,m]);
      real cur_log_prob_off = log(1 - modelProbs[c,m]);

      target += log_sum_exp(cur_log_prob_on + normal_lpdf(baseCornerQuanities[,1,c,m]| prior_traj_mean, traj_k), // values
                            cur_log_prob_off + normal_lpdf(baseCornerQuanities[,1,c,m]| 0, off_traj_sd));

      target += log_sum_exp(cur_log_prob_on + normal_lpdf(baseCornerQuanities[,2,c,m]| 0, traj_k), //parX
                            cur_log_prob_off + normal_lpdf(baseCornerQuanities[,2,c,m]| 0, off_traj_sd));

      target += log_sum_exp(cur_log_prob_on + normal_lpdf(baseCornerQuanities[,3,c,m]| 0, traj_k), //parY
                            cur_log_prob_off + normal_lpdf(baseCornerQuanities[,3,c,m]| 0, off_traj_sd));

      target += log_sum_exp(cur_log_prob_on + normal_lpdf(baseCornerQuanities[,4,c,m]| 0, traj_k), //parXY
                            cur_log_prob_off + normal_lpdf(baseCornerQuanities[,4,c,m]| 0, off_traj_sd));

    }

  }

  // Likelihood

  for(i in 1:N_data){

    vector[2] cur_TrajWeightedVel = TrajWeightedBaseVectorFields(Data[i,1], Data[i,2:3]', compSpacePos_data[i, 1:2]',
                                                       cur_transBoundaries, trans_res, cur_Coefs, N_models);

    target += normal_lpdf(Data[i,4:5]| cur_TrajWeightedVel, sigma_vel);

  }

}

generated quantities {
    // Something here
}

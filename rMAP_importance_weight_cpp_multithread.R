library(Rcpp)
library(RcppArmadillo)

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Inline helper to calculate drift, spatial Jacobian (A), and parameter gradient (B)
inline void single_drift_jacobian_gradient(
    double x, double y, const arma::mat& beta_mat,
    const arma::vec& border, const arma::vec& omega,
    int M, int M_sq, double Lx, double Ly,
    double& dx, double& dy, arma::mat& A, arma::mat& B) {

    double log_T1 = 0.0, log_T2 = 0.0;
    double dT1_dx = 0.0, dT1_dy = 0.0;
    double dT2_dx = 0.0, dT2_dy = 0.0;

    double cx[128], cy[128], sx[128], sy[128];

    // O(M) Precalculation
    for(int i = 0; i <= M; ++i) {
        double scale = (i == 0) ? 1.0 : std::sqrt(2.0);
        double ang_x = omega[i] * (x - border[0]) / Lx;
        double ang_y = omega[i] * (y - border[1]) / Ly;

        cx[i] = scale * std::cos(ang_x);
        cy[i] = scale * std::cos(ang_y);
        sx[i] = -scale * std::sin(ang_x) * (omega[i] / Lx);
        sy[i] = -scale * std::sin(ang_y) * (omega[i] / Ly);
    }

    arma::vec Phi(M_sq, arma::fill::zeros);

    // O(M^2) computation for basis functions and their spatial derivatives
    for(int j = 0; j <= M; ++j) {
        for(int i = 0; i <= M; ++i) {
            int idx = j * (M + 1) + i;

            double phi_val = cx[i] * cy[j];
            double dphi_dx = sx[i] * cy[j];
            double dphi_dy = cx[i] * sy[j];

            Phi[idx] = phi_val;

            double b1 = beta_mat(idx, 0);
            double b2 = beta_mat(idx, 1);

            log_T1 += phi_val * b1;
            log_T2 += phi_val * b2;

            dT1_dx += dphi_dx * b1;
            dT1_dy += dphi_dy * b1;
            dT2_dx += dphi_dx * b2;
            dT2_dy += dphi_dy * b2;
        }
    }

    // Base Vector Fields and their Spatial Derivatives
    double c = std::sqrt(x*x + y*y); // + 1e-6;
    double c3 = c * c * c;

    double f1x = y / c;
    double f1y = -x / c;
    double f2x = x / c;
    double f2y = y / c;

    double df1x_dx = -(x * y) / c3;
    double df1x_dy =  (x * x) / c3;
    double df1y_dx = -(y * y) / c3;
    double df1y_dy =  (x * y) / c3;

    double df2x_dx =  (y * y) / c3;
    double df2x_dy = -(x * y) / c3;
    double df2y_dx = -(x * y) / c3;
    double df2y_dy =  (x * x) / c3;

    double eT1 = std::exp(log_T1);
    double eT2 = std::exp(log_T2);

    // Drift Output
    dx = eT1 * f1x + eT2 * f2x;
    dy = eT1 * f1y + eT2 * f2y;

    // Spatial Jacobian A (2x2)
    A(0,0) = eT1 * dT1_dx * f1x + eT1 * df1x_dx + eT2 * dT2_dx * f2x + eT2 * df2x_dx;
    A(0,1) = eT1 * dT1_dy * f1x + eT1 * df1x_dy + eT2 * dT2_dy * f2x + eT2 * df2x_dy;
    A(1,0) = eT1 * dT1_dx * f1y + eT1 * df1y_dx + eT2 * dT2_dx * f2y + eT2 * df2y_dx;
    A(1,1) = eT1 * dT1_dy * f1y + eT1 * df1y_dy + eT2 * dT2_dy * f2y + eT2 * df2y_dy;

    // Parameter Gradient B (2 x 2*M_sq)
    for(int k = 0; k < M_sq; ++k) {
        B(0, k)        = f1x * eT1 * Phi[k];
        B(1, k)        = f1y * eT1 * Phi[k];
        B(0, M_sq + k) = f2x * eT2 * Phi[k];
        B(1, M_sq + k) = f2y * eT2 * Phi[k];
    }
}


// [[Rcpp::export]]
double calculate_log_importance_cpp(
    const arma::vec& beta,
    int M_sq,
    const arma::mat& start_t_pos_mat,
    const arma::mat& end_t_pos_true_mat,
    const arma::vec& t_steps,
    int N_prop_steps,
    const arma::vec& border,
    double pos_sd,
    const arma::mat& prior_beta_sigma,
    int n_threads = 8) {

    // Reconstruct Beta Matrix (M_sq x 2)
    arma::vec beta_1 = beta.subvec(0, M_sq - 1);
    arma::vec beta_2 = beta.subvec(M_sq, 2 * M_sq - 1);
    arma::mat beta_mat = arma::join_rows(beta_1, beta_2);

    int N_data = start_t_pos_mat.n_rows;
    int N_params = 2 * M_sq;
    double Lx = border[2] - border[0];
    double Ly = border[3] - border[1];
    int M = std::round(std::sqrt(M_sq)) - 1;
    arma::vec omega = arma::regspace(0, M) * arma::datum::pi;

    // Output Memory
    arma::mat del_G = arma::zeros(2 * N_data, N_params);
    arma::vec end_t_pos_prop = arma::zeros(2 * N_data);

    // Flatten true end positions (interleaving x and y)
    arma::vec end_t_pos_true = arma::vectorise(end_t_pos_true_mat.cols(1, 2).t());

#ifdef _OPENMP
    if(n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
#endif

    // ---------------------------------------------------------
    // MULTITHREADED RK4 PROPAGATION (State & Sensitivity)
    // ---------------------------------------------------------
    #pragma omp parallel for schedule(static)
    for(int r = 0; r < N_data; ++r) {

        double cur_x = start_t_pos_mat(r, 1);
        double cur_y = start_t_pos_mat(r, 2);
        double dt = t_steps[r] / N_prop_steps;
        double dt_half = dt / 2.0;

        arma::mat cur_g = arma::zeros(2, N_params);

        // Thread-local temporary matrices
        arma::mat A1(2,2), A2(2,2), A3(2,2), A4(2,2);
        arma::mat B1(2, N_params), B2(2, N_params), B3(2, N_params), B4(2, N_params);

        for(int j = 0; j < N_prop_steps; ++j) {
            double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;

            // k1
            single_drift_jacobian_gradient(cur_x, cur_y, beta_mat, border, omega, M, M_sq, Lx, Ly, k1x, k1y, A1, B1);
            arma::mat k1_g = A1 * cur_g + B1;

            // k2
            arma::mat g2 = cur_g + k1_g * dt_half;
            single_drift_jacobian_gradient(cur_x + k1x*dt_half, cur_y + k1y*dt_half, beta_mat, border, omega, M, M_sq, Lx, Ly, k2x, k2y, A2, B2);
            arma::mat k2_g = A2 * g2 + B2;

            // k3
            arma::mat g3 = cur_g + k2_g * dt_half;
            single_drift_jacobian_gradient(cur_x + k2x*dt_half, cur_y + k2y*dt_half, beta_mat, border, omega, M, M_sq, Lx, Ly, k3x, k3y, A3, B3);
            arma::mat k3_g = A3 * g3 + B3;

            // k4
            arma::mat g4 = cur_g + k3_g * dt;
            single_drift_jacobian_gradient(cur_x + k3x*dt, cur_y + k3y*dt, beta_mat, border, omega, M, M_sq, Lx, Ly, k4x, k4y, A4, B4);
            arma::mat k4_g = A4 * g4 + B4;

            // Update State
            cur_x += (k1x + 2.0*k2x + 2.0*k3x + k4x) * dt / 6.0;
            cur_y += (k1y + 2.0*k2y + 2.0*k3y + k4y) * dt / 6.0;

            // Update Sensitivity
            cur_g += (k1_g + 2.0*k2_g + 2.0*k3_g + k4_g) * dt / 6.0;
        }

        // Safely write results back to shared memory
        end_t_pos_prop(2*r)     = cur_x;
        end_t_pos_prop(2*r + 1) = cur_y;
        del_G.rows(2*r, 2*r + 1) = cur_g;
    }

    // ---------------------------------------------------------
    // FAST LINEAR ALGEBRA ALONG rMAP WEIGHT FORMULA
    // ---------------------------------------------------------

    // 1. Construct Block Diagonal Prior Covariance C in C++
    arma::mat C = arma::zeros(N_params, N_params);
    C.submat(0, 0, M_sq - 1, M_sq - 1) = prior_beta_sigma;
    C.submat(M_sq, M_sq, N_params - 1, N_params - 1) = prior_beta_sigma;

    double var_inv = 1.0 / (pos_sd * pos_sd);

    // 2. Scaled Misfit Vector (K)
    arma::vec K = ((end_t_pos_prop - end_t_pos_true) + del_G * beta) * var_inv;

    // 3. Precision Matrix Inverse (H)
    arma::mat L_inv = var_inv * arma::eye(2 * N_data, 2 * N_data);
    arma::mat H_inv = L_inv + (del_G * C * del_G.t()) * (var_inv * var_inv);

    // inv_sympd is hyper-optimized for symmetric positive definite matrices
    arma::mat H = arma::inv_sympd(H_inv);

    // 4. Gauss-Newton Jacobian Determinant (|J|)
    arma::mat J_approx = arma::eye(N_params, N_params) + (C * del_G.t() * del_G) * var_inv;

    double log_det_J;
    double sign;
    arma::log_det(log_det_J, sign, J_approx);

    // 5. Final Log Importance Weight
    arma::mat K_t_H_K = K.t() * H * K;
    double log_importance = -0.5 * K_t_H_K(0,0) - 0.5 * log_det_J;

    return log_importance;
}
'

Rcpp::sourceCpp(code = cpp_code)

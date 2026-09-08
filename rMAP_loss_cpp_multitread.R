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

// Inline helper to calculate drift for a SINGLE particle.
// This completely avoids allocating the large N x num_bases Phi matrix.
inline void single_drift(double x, double y, const arma::mat& beta_mat,
                         const arma::vec& border, const arma::vec& omega,
                         int M, double Lx, double Ly, double& dx, double& dy) {

    double log_T1 = 0.0;
    double log_T2 = 0.0;

    // Use stack arrays to avoid expensive heap allocation inside the RK4 loop
    // 128 supports up to 16,384 basis functions. Increase if your M > 127.
    double cx[128];
    double cy[128];

    // O(M) Precalculation: Remove transcendental math from the nested loop
    for(int i = 0; i <= M; ++i) {
        double scale = (i == 0) ? 1.0 : std::sqrt(2.0);
        cx[i] = scale * std::cos(omega[i] * (x - border[0]) / Lx);
        cy[i] = scale * std::cos(omega[i] * (y - border[1]) / Ly);
    }

    // O(M^2) computation using only simple multiplication
    // Note: i is the inner loop to ensure consecutive memory access for col_idx
    for(int j = 0; j <= M; ++j) {
        for(int i = 0; i <= M; ++i) {
            double phi_val = cx[i] * cy[j];
            int col_idx = j * (M + 1) + i; // Now increments by exactly 1 per loop

            log_T1 += phi_val * beta_mat(col_idx, 0);
            log_T2 += phi_val * beta_mat(col_idx, 1);
        }
    }

    // Base Vector Field + Weighting
    double c = std::sqrt(x*x + y*y) + 1e-6;
    double f1x = y / c;
    double f1y = -x / c;
    double f2x = x / c;
    double f2y = y / c;

    double exp_L1 = std::exp(log_T1);
    double exp_L2 = std::exp(log_T2);

    dx = exp_L1 * f1x + exp_L2 * f2x;
    dy = exp_L1 * f1y + exp_L2 * f2y;
    }


// [[Rcpp::export]]
double rMAP_loss_cpp_multi(const arma::vec& beta,
                     int M_sq,
                     const arma::mat& aug_data_starts,
                     const arma::mat& aug_data_ends,
                     const arma::vec& t_steps,
                     const arma::mat& rand_vel_1,
                     const arma::mat& rand_vel_2,
                     const arma::vec& border,
                     double pos_sd,
                     const arma::mat& prior_beta_sigma,
                     const arma::vec& start_beta,
                     int n_threads = 8) {

    // Reconstruct Beta Matrices
    arma::vec beta_1 = beta.subvec(0, M_sq - 1);
    arma::vec beta_2 = beta.subvec(M_sq, 2 * M_sq - 1);
    arma::mat beta_mat = arma::join_rows(beta_1, beta_2);

    int N = aug_data_starts.n_rows;
    int N_prop_steps = rand_vel_1.n_rows;
    arma::vec sqrt_t_steps = arma::sqrt(t_steps);

    arma::mat curPos = aug_data_starts;

    double Lx = border[2] - border[0];
    double Ly = border[3] - border[1];
    int num_bases = beta_mat.n_rows;
    int M = std::round(std::sqrt(num_bases)) - 1;
    arma::vec omega = arma::regspace(0, M) * arma::datum::pi;

#ifdef _OPENMP
    if(n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
#endif

    // Loop Inversion: Parallelize over independent particles
    #pragma omp parallel for schedule(static)
    for(int r = 0; r < N; ++r) {
        double cur_t = aug_data_starts(r, 0);
        double cur_x = aug_data_starts(r, 1);
        double cur_y = aug_data_starts(r, 2);

        double dt = t_steps[r];
        double sdt = sqrt_t_steps[r];
        double dt_half = dt / 2.0;

        // Propagate this specific particle through all time steps
        for(int j = 0; j < N_prop_steps; ++j) {
            double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;

            // k1
            single_drift(cur_x, cur_y, beta_mat, border, omega, M, Lx, Ly, k1x, k1y);

            // k2
            single_drift(cur_x + k1x * dt_half, cur_y + k1y * dt_half,
                         beta_mat, border, omega, M, Lx, Ly, k2x, k2y);

            // k3
            single_drift(cur_x + k2x * dt_half, cur_y + k2y * dt_half,
                         beta_mat, border, omega, M, Lx, Ly, k3x, k3y);

            // k4
            single_drift(cur_x + k3x * dt, cur_y + k3y * dt,
                         beta_mat, border, omega, M, Lx, Ly, k4x, k4y);

            double rk4_dx = (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
            double rk4_dy = (k1y + 2.0*k2y + 2.0*k3y + k4y) / 6.0;

            cur_t += dt;
            // Accessing rand_vel matrices safely: rand_vel(j, r)
            cur_x += rk4_dx * dt + rand_vel_1(j, r) * sdt;
            cur_y += rk4_dy * dt + rand_vel_2(j, r) * sdt;
        }

        // Write final result back to memory
        curPos(r, 0) = cur_t;
        curPos(r, 1) = cur_x;
        curPos(r, 2) = cur_y;
    }

    // 4. Calculate Final Loss
    arma::vec diff_X = curPos.col(1) - aug_data_ends.col(1);
    arma::vec diff_Y = curPos.col(2) - aug_data_ends.col(2);
    double likelihood_loss = (arma::dot(diff_X, diff_X) + arma::dot(diff_Y, diff_Y)) / (2.0 * pos_sd * pos_sd);

    arma::vec d_beta1 = beta_1 - start_beta.subvec(0, M_sq - 1);
    arma::vec d_beta2 = beta_2 - start_beta.subvec(M_sq, 2 * M_sq - 1);

    double prior_loss = arma::dot(d_beta1, prior_beta_sigma * d_beta1) +
                        arma::dot(d_beta2, prior_beta_sigma * d_beta2);

    return likelihood_loss + prior_loss;
}
'

Rcpp::sourceCpp(code = cpp_code)

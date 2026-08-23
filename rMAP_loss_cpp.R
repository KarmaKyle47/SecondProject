library(Rcpp)
library(RcppArmadillo)

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper to calculate the combined weighted vector field
arma::mat calc_drift(const arma::mat& pos_t, const arma::mat& beta_mat, const arma::vec& border) {
    int N = pos_t.n_rows;
    arma::mat drift(N, 2);
    double Lx = border[2] - border[0];
    double Ly = border[3] - border[1];

    int num_bases = beta_mat.n_rows;
    int M = std::round(std::sqrt(num_bases)) - 1;
    arma::vec omega = arma::regspace(0, M) * arma::datum::pi;

    // 1. Calculate Phi matrix
    arma::mat Phi(N, num_bases, arma::fill::ones);

    for(int i = 0; i <= M; ++i) {
        for(int j = 0; j <= M; ++j) {
            double scale = 1.0;
            if(i != 0) scale *= std::sqrt(2.0);
            if(j != 0) scale *= std::sqrt(2.0);

            int col_idx = i * (M + 1) + j;
            for(int r = 0; r < N; ++r) {
                double cx = std::cos(omega[i] * (pos_t(r, 1) - border[0]) / Lx);
                double cy = std::cos(omega[j] * (pos_t(r, 2) - border[1]) / Ly);
                Phi(r, col_idx) = scale * cx * cy;
            }
        }
    }

    arma::mat log_Traj = Phi * beta_mat;

    // 2. Base Vector Field + Weighting
    for(int r = 0; r < N; ++r) {
        double x = pos_t(r, 1);
        double y = pos_t(r, 2);
        // Added 1e-6 epsilon to prevent division by zero near origin
        double c = std::sqrt(x*x + y*y) + 1e-6;

        double f1x = y / c;
        double f1y = -x / c;
        double f2x = x / c;
        double f2y = y / c;

        double exp_L1 = std::exp(log_Traj(r, 0));
        double exp_L2 = std::exp(log_Traj(r, 1));

        drift(r, 0) = exp_L1 * f1x + exp_L2 * f2x;
        drift(r, 1) = exp_L1 * f1y + exp_L2 * f2y;
    }

    return drift;
}

// [[Rcpp::export]]
double rMAP_loss_cpp(arma::vec beta,
                     int M_sq,
                     arma::mat aug_data_starts,
                     arma::mat aug_data_ends,
                     arma::vec t_steps,
                     arma::mat rand_vel_1,
                     arma::mat rand_vel_2,
                     arma::vec border,
                     double pos_sd,
                     arma::mat prior_beta_sigma,
                     arma::vec start_beta) {

    // Reconstruct Beta Matrices
    arma::vec beta_1 = beta.subvec(0, M_sq - 1);
    arma::vec beta_2 = beta.subvec(M_sq, 2 * M_sq - 1);
    arma::mat beta_mat = arma::join_rows(beta_1, beta_2);

    int N = aug_data_starts.n_rows;
    int N_prop_steps = rand_vel_1.n_rows;
    arma::vec sqrt_t_steps = arma::sqrt(t_steps);

    arma::mat curPos = aug_data_starts;

    // 3. RK4 Propagation Loop
    for(int j = 0; j < N_prop_steps; ++j) {
        arma::mat k1 = calc_drift(curPos, beta_mat, border);

        // In Armadillo, % does element-wise multiplication
        arma::mat pos_k2 = curPos;
        pos_k2.col(0) += t_steps / 2.0;
        pos_k2.col(1) += k1.col(0) % (t_steps / 2.0);
        pos_k2.col(2) += k1.col(1) % (t_steps / 2.0);
        arma::mat k2 = calc_drift(pos_k2, beta_mat, border);

        arma::mat pos_k3 = curPos;
        pos_k3.col(0) += t_steps / 2.0;
        pos_k3.col(1) += k2.col(0) % (t_steps / 2.0);
        pos_k3.col(2) += k2.col(1) % (t_steps / 2.0);
        arma::mat k3 = calc_drift(pos_k3, beta_mat, border);

        arma::mat pos_k4 = curPos;
        pos_k4.col(0) += t_steps;
        pos_k4.col(1) += k3.col(0) % t_steps;
        pos_k4.col(2) += k3.col(1) % t_steps;
        arma::mat k4 = calc_drift(pos_k4, beta_mat, border);

        arma::mat rk4_drift = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;

        arma::rowvec diff1 = rand_vel_1.row(j);
        arma::rowvec diff2 = rand_vel_2.row(j);

        curPos.col(0) += t_steps;
        curPos.col(1) += rk4_drift.col(0) % t_steps + diff1.t() % sqrt_t_steps;
        curPos.col(2) += rk4_drift.col(1) % t_steps + diff2.t() % sqrt_t_steps;
    }

    // 4. Calculate Final Loss
    arma::vec diff_X = curPos.col(1) - aug_data_ends.col(1);
    arma::vec diff_Y = curPos.col(2) - aug_data_ends.col(2);
    double likelihood_loss = (arma::dot(diff_X, diff_X) + arma::dot(diff_Y, diff_Y)) / (2.0 * N * pos_sd * pos_sd);

    arma::vec d_beta1 = beta_1 - start_beta.subvec(0, M_sq - 1);
    arma::vec d_beta2 = beta_2 - start_beta.subvec(M_sq, 2 * M_sq - 1);

    double prior_loss = arma::dot(d_beta1, prior_beta_sigma * d_beta1) +
                        arma::dot(d_beta2, prior_beta_sigma * d_beta2);

    return likelihood_loss + prior_loss;
}
'

Rcpp::sourceCpp(code = cpp_code)

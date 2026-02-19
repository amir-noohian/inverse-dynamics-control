// compare_torque.cpp  (4-DoF input version, no temp files)
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cassert>
#include <cstdio>

#include "W.hpp"           // calculate_W(...) + calculate_W_eigen(...)
#include "beta_ares_link.hpp"   // initialize_beta()

struct Dataset {
    Eigen::VectorXd t;        // N
    Eigen::MatrixXd q4;       // N x 4
    Eigen::MatrixXd tau4;     // N x 4 (optional)
};

struct Derivs {
    Eigen::MatrixXd dq4;      // N x 4
    Eigen::MatrixXd ddq4;     // N x 4
};

struct Metrics {
    Eigen::RowVector4d rmse;
    Eigen::RowVector4d mae;
    Eigen::RowVector4d r2;
    Eigen::RowVector4d nrmse;
    double rmse_overall;
};

// ------------------ Loader: expects at least t + 4 q; optional 4 tau ------------------
static Dataset load_dat_4dof(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open file: " + path);

    std::vector<double> t;
    std::vector<std::array<double,4>> qv;
    std::vector<std::array<double,4>> tauv;

    std::string line;
    auto parse = [](const std::string& s) {
        std::stringstream ss(s);
        std::vector<double> vals; double v;
        while (ss >> v) vals.push_back(v);
        return vals;
    };

    while (std::getline(fin, line)) {
        // skip blank/comment
        bool blank = true;
        for (char c : line) if (!std::isspace(static_cast<unsigned char>(c))) { blank=false; break; }
        if (blank) continue;
        if (line[0]=='#' || line[0]=='%') continue;

        auto vals = parse(line);
        if (vals.size() < 5) continue; // need t + 4 q
        t.push_back(vals[0]);

        std::array<double,4> qrow{vals[1], vals[2], vals[3], vals[4]};
        qv.push_back(qrow);

        if (vals.size() >= 9) {
            std::array<double,4> taurow{vals[5], vals[6], vals[7], vals[8]};
            tauv.push_back(taurow);
        }
        // if more columns exist, they are ignored
    }
    if (t.empty()) throw std::runtime_error("No valid rows parsed: " + path);

    // shift time so it starts at zero
    double t0 = t.front();
    for (double& ti : t) ti -= t0;

    const int N = static_cast<int>(t.size());
    Dataset D;
    D.t = Eigen::Map<Eigen::VectorXd>(t.data(), N);

    D.q4.resize(N,4);
    for (int i=0;i<N;++i)
        for (int j=0;j<4;++j)
            D.q4(i,j) = qv[i][j];

    if (!tauv.empty() && static_cast<int>(tauv.size()) == N) {
        D.tau4.resize(N,4);
        for (int i=0;i<N;++i)
            for (int j=0;j<4;++j)
                D.tau4(i,j) = tauv[i][j];
    } else {
        D.tau4.resize(0,0);
    }
    return D;
}

// ------------------ LPF & derivatives ------------------
// First-order exponential LPF per column:
// y[n] = alpha*x[n] + (1-alpha)*y[n-1], alpha = dt/(tau+dt), tau = 1/(2π fc)
static Eigen::MatrixXd lowpass_filter(const Eigen::VectorXd& t,
                                      const Eigen::MatrixXd& X,
                                      double fc_hz) {
    const int N = X.rows(), C = X.cols();
    Eigen::MatrixXd Y = X;
    if (N <= 1 || fc_hz <= 0.0) return Y;

    const double tau = 1.0 / (2.0 * M_PI * fc_hz);
    for (int j = 0; j < C; ++j) {
        double yprev = X(0, j);
        Y(0, j) = yprev;
        for (int i = 1; i < N; ++i) {
            double dt = t(i) - t(i - 1);
            double alpha = dt / (tau + dt);
            double y = alpha * X(i, j) + (1.0 - alpha) * yprev;
            Y(i, j) = y;
            yprev = y;
        }
    }
    return Y;
}

// 1) dq from q, 2) LPF(dq) at fc, 3) ddq from filtered dq
static Derivs differentiate4_filtered(const Eigen::VectorXd& t,
                                      const Eigen::MatrixXd& q4,
                                      double fc_hz) {
    const int N = static_cast<int>(t.size());
    Eigen::MatrixXd dq_raw = Eigen::MatrixXd::Zero(N, 4);
    Eigen::MatrixXd ddq    = Eigen::MatrixXd::Zero(N, 4);
    if (N <= 1) return {dq_raw, ddq};

    // velocities from positions (backward/forward at ends)
    for (int i = 1; i < N; ++i) {
        double dt = t(i) - t(i - 1);
        for (int j = 0; j < 4; ++j)
            dq_raw(i, j) = (q4(i, j) - q4(i - 1, j)) / dt;
    }
    dq_raw.row(0) = dq_raw.row(1);

    // filter velocities
    Eigen::MatrixXd dq = lowpass_filter(t, dq_raw, fc_hz);

    // accelerations from filtered velocities
    for (int i = 1; i < N; ++i) {
        double dt = t(i) - t(i - 1);
        for (int j = 0; j < 4; ++j)
            ddq(i, j) = (dq(i, j) - dq(i - 1, j)) / dt;
    }
    ddq.row(0) = ddq.row(1);

    return {dq, ddq};
}

// ------------------ Tau prediction ------------------
static Eigen::MatrixXd predict_tau4(const Eigen::MatrixXd& q4,
                                    const Eigen::MatrixXd& dq4,
                                    const Eigen::MatrixXd& ddq4,
                                    const Eigen::Matrix<double,30,1>& beta) {
    const int N = q4.rows();
    Eigen::MatrixXd tau_pred(N,4);
    for (int i=0;i<N;++i) {
        Eigen::Vector4d qi  = q4.row(i).transpose();
        Eigen::Vector4d dqi = dq4.row(i).transpose();
        Eigen::Vector4d ddqi= ddq4.row(i).transpose();
        auto Wi = calculate_W_eigen(qi, dqi, ddqi);     // 4x30
        tau_pred.row(i) = (Wi * beta).transpose();      // 1x4
    }
    return tau_pred;
}

// ------------------ Metrics ------------------
static Metrics compare_metrics(const Eigen::MatrixXd& y_true,
                               const Eigen::MatrixXd& y_pred) {
    const int N = y_true.rows(), J = y_true.cols();
    Metrics M; M.rmse.setZero(); M.mae.setZero(); M.r2.setZero(); M.nrmse.setZero();
    double sse_total = 0.0;

    for (int j=0;j<J;++j) {
        double sse=0.0, sae=0.0, sst=0.0;
        double ymin = y_true.col(j).minCoeff();
        double ymax = y_true.col(j).maxCoeff();
        double range = ymax - ymin;
        double mean  = y_true.col(j).mean();

        for (int i=0;i<N;++i) {
            double e = y_pred(i,j) - y_true(i,j);
            sse += e*e;
            sae += std::abs(e);
            double d = y_true(i,j) - mean;
            sst += d*d;
        }
        M.rmse(j) = std::sqrt(sse/std::max(1,N));
        M.mae(j)  = sae/std::max(1,N);
        M.r2(j)   = (sst>0.0)?(1.0 - sse/sst):1.0;

        double denom = range;
        if (denom < 1e-9) {
            double var = sst / std::max(1, N-1);
            denom = std::sqrt(std::max(0.0, var));
            if (denom < 1e-9) denom = 1.0;
        }
        M.nrmse(j) = M.rmse(j) / denom;
        sse_total += sse;
    }
    M.rmse_overall = std::sqrt(sse_total / double(std::max(1, N*J)));
    return M;
}

// ------------------ Plotting (gnuplot via pipe; no temp files) ------------------
static void plot_multi4_pipe(const std::vector<double>& t,
                             const std::array<std::vector<double>,4>& Y,
                             const std::string& title,
                             const std::string& ylabel,
                             const std::array<std::string,4>& legend) {
    FILE* gp = popen("gnuplot", "w"); // no -persist so it exits with program
    if (!gp) return;

    fprintf(gp, "set title '%s'\n", title.c_str());
    fprintf(gp, "set xlabel 'time [s]'\n");
    fprintf(gp, "set ylabel '%s'\n", ylabel.c_str());
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '-' w l title '%s', '-' w l title '%s', '-' w l title '%s', '-' w l title '%s'\n",
            legend[0].c_str(), legend[1].c_str(), legend[2].c_str(), legend[3].c_str());

    for (int j=0;j<4;++j) {
        for (size_t i=0;i<t.size();++i)
            fprintf(gp, "%f %f\n", t[i], Y[j][i]);
        fprintf(gp, "e\n");
    }
    fflush(gp);
    pclose(gp); // closes the plot as program ends
}

static void plot_pair_pipe(const std::vector<double>& t,
                           const std::vector<double>& y1,
                           const std::vector<double>& y2,
                           const std::string& title,
                           const std::string& ylabel) {
    FILE* gp = popen("gnuplot", "w"); // no -persist
    if (!gp) return;

    fprintf(gp, "set title '%s'\n", title.c_str());
    fprintf(gp, "set xlabel 'time [s]'\n");
    fprintf(gp, "set ylabel '%s'\n", ylabel.c_str());
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '-' w l title 'measured', '-' w l title 'predicted'\n");

    for (size_t i=0;i<t.size();++i) fprintf(gp, "%f %f\n", t[i], y1[i]);
    fprintf(gp, "e\n");
    for (size_t i=0;i<t.size();++i) fprintf(gp, "%f %f\n", t[i], y2[i]);
    fprintf(gp, "e\n");

    fflush(gp);
    pclose(gp);
}

// ------------------ Main ------------------
int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " data_4dof.dat\n";
            return 1;
        }

        // Load 4-DoF data (t, q1..q4, [tau1..tau4]), shift t so t0=0
        Dataset D = load_dat_4dof(argv[1]);

        // Cutoff from your formula: fc = fc_mult * (wf/(2π)) * L
        const double fc_mult = 10.0;
        const double wf = 0.3141592653589793;  // rad/s
        const double L  = 5.0;                 // dimensionless
        const double fc = fc_mult * (wf / (2.0 * M_PI) * L); // 2.5 Hz for above values

        // Derivatives with velocity LPF at fc
        Derivs derivs = differentiate4_filtered(D.t, D.q4, fc);

        // Predict torques via W * beta
        auto beta = initialize_beta();
        Eigen::MatrixXd tau_pred = predict_tau4(D.q4, derivs.dq4, derivs.ddq4, beta);

        // Metrics (if measured torques present)
        if (D.tau4.size() > 0) {
            auto M = compare_metrics(D.tau4, tau_pred);
            std::cout << "Cutoff fc (Hz): " << fc << "\n";
            std::cout << "RMSE   : " << M.rmse << "\n";
            std::cout << "MAE    : " << M.mae  << "\n";
            std::cout << "R^2    : " << M.r2   << "\n";
            std::cout << "NRMSE  : " << M.nrmse << " (range-normalized)\n";
            std::cout << "Overall RMSE: " << M.rmse_overall << "\n";
        } else {
            std::cout << "Cutoff fc (Hz): " << fc << " (no measured torques found)\n";
        }

        // Prepare vectors for plotting
        const int N = static_cast<int>(D.t.size());
        std::vector<double> tv(N);
        for (int i=0;i<N;++i) tv[i] = D.t(i);

        std::array<std::vector<double>,4> qv, dqv, ddqv;
        for (int j=0;j<4;++j) {
            qv[j].resize(N); dqv[j].resize(N); ddqv[j].resize(N);
            for (int i=0;i<N;++i) {
                qv[j][i]   = D.q4(i,j);
                dqv[j][i]  = derivs.dq4(i,j);
                ddqv[j][i] = derivs.ddq4(i,j);
            }
        }

        // Plots (auto-close with program)
        plot_multi4_pipe(tv, qv,   "Joint positions (4-DoF)",       "position [rad]",
                         {"q1","q2","q3","q4"});
        plot_multi4_pipe(tv, dqv,  "Joint velocities (LPF, 4-DoF)", "velocity [rad/s]",
                         {"dq1","dq2","dq3","dq4"});
        plot_multi4_pipe(tv, ddqv, "Joint accelerations (4-DoF)",   "acc [rad/s^2]",
                         {"ddq1","ddq2","ddq3","ddq4"});

        if (D.tau4.size() > 0) {
            for (int j=0;j<4;++j) {
                std::vector<double> ym(N), yp(N);
                for (int i=0;i<N;++i) { ym[i] = D.tau4(i,j); yp[i] = tau_pred(i,j); }
                plot_pair_pipe(tv, ym, yp, "Torque comparison - J" + std::to_string(j+1),
                               "torque [Nm]");
            }
        } else {
            // predicted-only plots (vs zero) — comment out if not needed
            for (int j=0;j<4;++j) {
                std::vector<double> zero(N, 0.0), yp(N);
                for (int i=0;i<N;++i) yp[i] = tau_pred(i,j);
                plot_pair_pipe(tv, zero, yp, "Torque (pred only) - J" + std::to_string(j+1),
                               "torque [Nm]");
            }
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}

// compare_torque.cpp
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "W.hpp"
#include "beta_slax.hpp"

struct Dataset {
    Eigen::VectorXd t;
    Eigen::MatrixXd q4;
    Eigen::MatrixXd tau4;
};

struct Derivs {
    Eigen::MatrixXd dq4;
    Eigen::MatrixXd ddq4;
};

struct Metrics {
    Eigen::RowVector4d rmse;
    Eigen::RowVector4d mae;
    Eigen::RowVector4d r2;
    Eigen::RowVector4d nrmse;
    double rmse_overall;
};

// ------------------ Loader ------------------
static Dataset load_dat_first4dof(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open file: " + path);

    std::vector<double> t;
    std::vector<std::array<double,7>> q7, tau7;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0]=='#' || line[0]=='%') continue;
        std::stringstream ss(line);
        std::vector<double> vals;
        double v; while (ss >> v) vals.push_back(v);
        if (vals.size() < 8) continue;
        t.push_back(vals[0]);
        std::array<double,7> qrow{};
        for (int j=0;j<7;++j) qrow[j]=vals[1+j];
        q7.push_back(qrow);
        if (vals.size()>=15) {
            std::array<double,7> taurow{};
            for (int j=0;j<7;++j) taurow[j]=vals[8+j];
            tau7.push_back(taurow);
        }
    }

    double t0 = t.front();
    for (auto& ti: t) ti -= t0;

    int N = t.size();
    Dataset D;
    D.t = Eigen::Map<Eigen::VectorXd>(t.data(), N);
    D.q4.resize(N,4);
    for (int i=0;i<N;++i)
        for (int j=0;j<4;++j)
            D.q4(i,j)=q7[i][j];
    if (!tau7.empty() && tau7.size()==t.size()) {
        D.tau4.resize(N,4);
        for (int i=0;i<N;++i)
            for (int j=0;j<4;++j)
                D.tau4(i,j)=tau7[i][j];
    }
    return D;
}

// ------------------ LPF & derivatives ------------------
static Eigen::MatrixXd lowpass_filter(const Eigen::VectorXd& t,
                                      const Eigen::MatrixXd& X,
                                      double fc_hz) {
    int N=X.rows();
    Eigen::MatrixXd Y=X;
    if (N<=1||fc_hz<=0.0) return Y;
    double tau=1.0/(2.0*M_PI*fc_hz);
    for(int j=0;j<X.cols();++j){
        double yprev=X(0,j);
        Y(0,j)=yprev;
        for(int i=1;i<N;++i){
            double dt=t(i)-t(i-1);
            double alpha=dt/(tau+dt);
            double y=alpha*X(i,j)+(1-alpha)*yprev;
            Y(i,j)=y; yprev=y;
        }
    }
    return Y;
}

static Derivs differentiate4_filtered(const Eigen::VectorXd& t,
                                      const Eigen::MatrixXd& q4,
                                      double fc_hz) {
    int N=t.size();
    Eigen::MatrixXd dq_raw=Eigen::MatrixXd::Zero(N,4);
    Eigen::MatrixXd ddq=Eigen::MatrixXd::Zero(N,4);
    if(N<=1) return {dq_raw,ddq};
    for(int i=1;i<N;++i){
        double dt=t(i)-t(i-1);
        for(int j=0;j<4;++j)
            dq_raw(i,j)=(q4(i,j)-q4(i-1,j))/dt;
    }
    dq_raw.row(0)=dq_raw.row(1);
    Eigen::MatrixXd dq=lowpass_filter(t,dq_raw,fc_hz);
    for(int i=1;i<N;++i){
        double dt=t(i)-t(i-1);
        for(int j=0;j<4;++j)
            ddq(i,j)=(dq(i,j)-dq(i-1,j))/dt;
    }
    ddq.row(0)=ddq.row(1);
    return {dq,ddq};
}

// ------------------ Tau prediction ------------------
static Eigen::MatrixXd predict_tau4(const Eigen::MatrixXd& q4,
                                    const Eigen::MatrixXd& dq4,
                                    const Eigen::MatrixXd& ddq4,
                                    const Eigen::Matrix<double,30,1>& beta){
    int N=q4.rows();
    Eigen::MatrixXd tau_pred(N,4);
    for(int i=0;i<N;++i){
        Eigen::Vector4d qi=q4.row(i).transpose();
        Eigen::Vector4d dqi=dq4.row(i).transpose();
        Eigen::Vector4d ddqi=ddq4.row(i).transpose();
        auto Wi=calculate_W_eigen(qi,dqi,ddqi);
        tau_pred.row(i)=(Wi*beta).transpose();
    }
    return tau_pred;
}

// ------------------ Metrics ------------------
static Metrics compare_metrics(const Eigen::MatrixXd& y_true,
                               const Eigen::MatrixXd& y_pred){
    int N=y_true.rows(), J=y_true.cols();
    Metrics M; M.rmse.setZero(); M.mae.setZero(); M.r2.setZero(); M.nrmse.setZero();
    double sse_total=0.0;
    for(int j=0;j<J;++j){
        double sse=0,sae=0,sst=0;
        double ymin=y_true.col(j).minCoeff();
        double ymax=y_true.col(j).maxCoeff();
        double range=ymax-ymin;
        double mean=y_true.col(j).mean();
        for(int i=0;i<N;++i){
            double e=y_pred(i,j)-y_true(i,j);
            sse+=e*e; sae+=std::abs(e);
            double d=y_true(i,j)-mean; sst+=d*d;
        }
        M.rmse(j)=std::sqrt(sse/std::max(1,N));
        M.mae(j)=sae/std::max(1,N);
        M.r2(j)=(sst>0)?(1.0-sse/sst):1.0;
        double denom=range;
        if(denom<1e-9){
            double var=sst/std::max(1,N-1);
            denom=std::sqrt(std::max(0.0,var));
            if(denom<1e-9) denom=1.0;
        }
        M.nrmse(j)=M.rmse(j)/denom;
        sse_total+=sse;
    }
    M.rmse_overall=std::sqrt(sse_total/double(std::max(1,N*J)));
    return M;
}

// ------------------ Plotting (pipes, auto-close) ------------------
static void plot_pair_pipe(const std::vector<double>& t,
                           const std::vector<double>& y1,
                           const std::vector<double>& y2,
                           const std::string& title,
                           const std::string& ylabel){
    FILE* gp=popen("gnuplot -persist","w");
    if(!gp) return;
    fprintf(gp,"set title '%s'\n",title.c_str());
    fprintf(gp,"set xlabel 'time [s]'\n");
    fprintf(gp,"set ylabel '%s'\n",ylabel.c_str());
    fprintf(gp,"set grid\n");
    fprintf(gp,"plot '-' with lines title 'measured', '-' with lines title 'predicted'\n");
    for(size_t i=0;i<t.size();++i) fprintf(gp,"%f %f\n",t[i],y1[i]);
    fprintf(gp,"e\n");
    for(size_t i=0;i<t.size();++i) fprintf(gp,"%f %f\n",t[i],y2[i]);
    fprintf(gp,"e\n");
    fflush(gp);
    pclose(gp); // closes when program ends
}

// ------------------ Main ------------------
int main(int argc,char** argv){
    try{
        if(argc<2){ std::cerr<<"Usage: "<<argv[0]<<" data.dat\n"; return 1; }
        Dataset D=load_dat_first4dof(argv[1]);

        const double wf=0.3141592653589793; // rad/s
        const double L=5.0;
        const double fc_mult=10.0;
        double fc=fc_mult*(wf/(2*M_PI))*L; // =2.5 Hz
        Derivs derivs=differentiate4_filtered(D.t,D.q4,fc);

        auto beta=initialize_beta();
        Eigen::MatrixXd tau_pred=predict_tau4(D.q4,derivs.dq4,derivs.ddq4,beta);

        if(D.tau4.size()>0){
            auto M=compare_metrics(D.tau4,tau_pred);
            std::cout<<"Cutoff fc (Hz): "<<fc<<"\n";
            std::cout<<"RMSE   : "<<M.rmse<<"\n";
            std::cout<<"MAE    : "<<M.mae<<"\n";
            std::cout<<"R^2    : "<<M.r2<<"\n";
            std::cout<<"NRMSE  : "<<M.nrmse<<"\n";
            std::cout<<"Overall RMSE: "<<M.rmse_overall<<"\n";
        }

        int N=D.t.size();
        std::vector<double> tv(N);
        for(int i=0;i<N;++i) tv[i]=D.t(i);

        if(D.tau4.size()>0){
            for(int j=0;j<4;++j){
                std::vector<double> ym(N),yp(N);
                for(int i=0;i<N;++i){ ym[i]=D.tau4(i,j); yp[i]=tau_pred(i,j); }
                plot_pair_pipe(tv,ym,yp,"Torque J"+std::to_string(j+1),"torque [Nm]");
            }
        }

        return 0;
    }catch(const std::exception& e){
        std::cerr<<"Error: "<<e.what()<<"\n";
        return 2;
    }
}

//Estimate P vector, only for non-gravity term

#ifndef PI_VECTOR_HPP_
#define PI_VECTOR_HPP_

#include <Eigen/Dense>
Eigen::Matrix<double, 8, 1> initialize_pi() {
    Eigen::Matrix<double, 8, 1> pi;
    pi << 0.147785256061085,
    0.209451284205178,
    0.0131079080258135,
    1.42981178937137,
    2.43011414756261,
    0.403307288959893,
    0.675223836841805,
    0.587854485146857;

    return pi;
}

#endif
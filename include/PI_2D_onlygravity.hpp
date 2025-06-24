//Estimate P vector, only for non-gravity term

#ifndef PI_VECTOR_HPP_
#define PI_VECTOR_HPP_

#include <Eigen/Dense>
Eigen::Matrix<double, 4, 1> initialize_pi() {
    Eigen::Matrix<double, 4, 1> pi;
    pi << -0.918337877170457,
          -0.103063209829940,
           0.0375697081564328,
          -2.21969748210226;

    return pi;
}

#endif
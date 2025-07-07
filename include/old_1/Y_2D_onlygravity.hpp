#ifndef Y_MATRIX_HPP_
#define Y_MATRIX_HPP_

#include <math.h>
#include <Eigen/Dense>

// Function to calculate the Y matrix
Eigen::Matrix<double, 2, 4> calculate_Y_matrix(const Eigen::Vector4d theta) {

    Eigen::Matrix<double, 2, 4> Y;
    Y.setZero(); // Initialize the matrix with zeros
    double g = 9.81;

    // Using 0-based indexing: theta[1] is theta2, theta[3] is theta4
    double theta2 = theta[1];
    double theta4 = theta[3];

    // First row of Y matrix (Y14, Y15, Y17, Y18)
    Y(0, 0) = g * sin(theta4 + theta2);     // Y14
    Y(0, 1) = -g * cos(theta4 + theta2);    // Y15
    Y(0, 2) = -cos(theta2) * g;             // Y17
    Y(0, 3) = sin(theta2) * g;              // Y18

    // Second row of Y matrix (Y24, Y25, Y27, Y28)
    Y(1, 0) = g * sin(theta4 + theta2);     // Y24
    Y(1, 1) = -g * cos(theta4 + theta2);    // Y25
    Y(1, 2) = 0.0;                          // Y27
    Y(1, 3) = 0.0;                          // Y28

    return Y;
}

#endif // Y_MATRIX_HPP_
#include "Filterable.h"


/// Kalman predict step
void Filterable::predict(VectorXd &x, MatrixXd &P) const {
    x = A * x;
    P = A * P * A.transpose() + Q;
};

/// Kalman update step
void Filterable::update(VectorXd &x, MatrixXd &P, const VectorXd &z) const {
    MatrixXd K = P * H.transpose() * (H * P * H.transpose() + R).inverse();
    MatrixXd I = MatrixXd::Identity(4, 4);

    x = x + K * y;
    P = (I - K * H) * P;
};


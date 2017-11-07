#ifndef FILTER_H_
#define FILTER_H_

#include "Eigen/Dense"


using namespace Eigen;


/// The Kalman filter protocol
struct Filterable {
public:
    /// State transition matrix
    MatrixXd A;

    /// State covariance matrix
    MatrixXd Q;
    
    /// Process measurement matrix
    MatrixXd H;
    
    /// Process covariance matrix
    MatrixXd R;
    
    /// Process measurement error
    VectorXd y;
protected:
    /// Kalman predict step
    void predict(VectorXd &x, MatrixXd &P) const;
    
    /// Kalman update step
    void update(VectorXd &x, MatrixXd &P, const VectorXd &z) const;
};

#endif

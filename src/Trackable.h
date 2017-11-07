#ifndef Trackable_hpp
#define Trackable_hpp


#include "Eigen/Dense"
#include "tools.h"


using namespace Eigen;


/// The Kalman filter
struct Trackable {
public:
    /// Process measurement matrix
    MatrixXd H;
    
    /// Process covariance matrix
    MatrixXd R;
    
    /// Predict state
    void predict(VectorXd &x, MatrixXd &P, VectorXd z, long delta);
    
    /// Update measurement
    virtual void update(VectorXd &x, MatrixXd &P, VectorXd z) = 0;
protected:
    /// State transition matrix
    MatrixXd A(long delta);
    
    /// State covariance matrix
    MatrixXd Q(long delta);
    
    /// Observation Jacobian matrix
    MatrixXd HJ(VectorXd z, VectorXd x);
    
    /// Identity matrix
    MatrixXd I = MatrixXd::Identity(4, 4);
};


void Trackable::predict(VectorXd &x, MatrixXd &P, VectorXd z, long delta) {
    MatrixXd A = this->A(delta);
    MatrixXd Q = this->Q(delta);
    x = A * x;
    P = A * P * A.transpose() + Q;
};


MatrixXd Trackable::A(long delta) {
    MatrixXd M(4, 4);
    M << 1, 0, delta, 0,
         0, 1, 0, delta,
         0, 0, 1, 0,
         0, 0, 0, 1;
    return M;
};


MatrixXd Trackable::Q(long delta) {
    double noise_ax = 9.0;
    double noise_ay = 9.0;
    
    float dt_2 = delta * delta;
    float dt_3 = dt_2 * delta;
    float dt_4 = dt_3 * delta;
    
    MatrixXd M(4, 4);
    M <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
    return M;
};


MatrixXd Trackable::HJ(VectorXd z, VectorXd x) {
    MatrixXd M(3, 4);
    
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);
    
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
    
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return M;
    }
    
    M << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return M;
};


struct Lidar: Trackable {
    Lidar() {
        H = MatrixXd(2, 4);
        R = MatrixXd(2, 2);
        H << 1, 0, 0, 0,
             0, 1, 0, 0;
        R << 0.0225, 0,
             0, 0.0225;
    };
    void update(VectorXd &x, MatrixXd &P, VectorXd z) {
        MatrixXd S = H * P * H.transpose() + R;
        MatrixXd K = P * H.transpose() * S.inverse();
        x = x + K * (z - H * x);
        P = (I * K * H) * P;
    };
};


struct Radar2: Trackable {
    Radar2() {
        H = MatrixXd(3, 4);
        R = MatrixXd(3, 3);
        R << 0.09, 0, 0,
             0, 0.0006, 0,
             0, 0, 0.09;
    };
    void update(VectorXd &x, MatrixXd &P, VectorXd z) {
        MatrixXd H = HJ(z, x);        
        MatrixXd S = H * P * H.transpose() + R;
        MatrixXd K = P * H.transpose() * S.inverse();
        
        VectorXd h(3);
        float px = x[0];
        float py = x[1];
        float vx = x[2];
        float vy = x[3];
        
        if( px == 0. && py == 0. )
            return;
        
        float rho = sqrt( px*px + py*py );
        h << rho, atan2( py, px), (px*vx + py*vy)/rho;
        
        VectorXd y = z - h;
        
        if( y[1] > PI )
            y[1] -= 2.f*PI;
        if( y[1] < -PI )
            y[1] += 2.f*PI;
        
        x = x + K * y;
        P = (I * K * H) * P;
    };
};


struct FusioManager {
public:
    /// Last measurement timestamp
    long long timestamp = 0;
    
    /// System state
    VectorXd x;
    
    /// System state covariance
    MatrixXd P;
    
    /// Predict and update
    void filter(Trackable t, VectorXd z, long long timestamp) {
        long delta = (timestamp - this->timestamp) / 1000000.0;
        this->timestamp = timestamp;
        t.predict(x, P, z, delta);
        t.update(x, P, z);
    };
};


#endif /* Trackable_hpp */

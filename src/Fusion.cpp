#include "Fusion.h"
#include <iostream>


void Fusion::process(const MeasurementPackage m) {
    if (!isInitialized) {
        initialize(m);
        return;
    };
    
    updateStateModel(m.timestamp_);
    predict(x, P);
    
    updateMeasurementModel(x, m.raw_measurements_, m.sensor_type_);
    update(x, P, m.raw_measurements_);
};


void Fusion::initialize(const MeasurementPackage m) {
    
    // Add initial state
    
    switch (m.sensor_type_) {
        case MeasurementPackage::RADAR: {
            
            // Covert polar to cartesian coordinates
            
            float ro = m.raw_measurements_[0];
            float theta = m.raw_measurements_[1];
            x << ro * cos(theta), ro * sin(theta), 0, 0;
            break;
        }
        case MeasurementPackage::LASER: {
            x << m.raw_measurements_[0], m.raw_measurements_[1], 0, 0;
            break;
        }
    }
    
    // Set state covariance model
    
    P << 1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;
    
    // Set state transition model
    
    A = MatrixXd(4, 4);
    A << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
    // Set measurement model
    
    H = MatrixXd(2, 4);
    H << 1, 0, 0, 0,
         0, 1, 0, 0;
    
    // Set state covariance model
    
    Q = MatrixXd(4, 4);
    
    // Define laser measurement covariance model
    
    laserR << 0.0225, 0,
              0, 0.0225;
    
    // Define radar measurement covariance model

    radarR << 0.09, 0, 0,
              0, 0.0006, 0,
              0, 0, 0.09;
    
    // Update initalization variable
    
    isInitialized = true;
    savedTimestamp = m.timestamp_;
}


void Fusion::updateStateModel(long timestamp) {
    
    // Calculate delta time
    
    double dt = (timestamp - savedTimestamp) / 1000000.0;
    savedTimestamp = timestamp;
    
    // Update state model

    A(0, 2) = dt;
    A(1, 3) = dt;
    
    // Update state covariace model
    
    double noise_ax = 9.0;
    double noise_ay = 9.0;
    
    float dt_2 = dt * dt;
    float dt_3 = pow(dt, 3);
    float dt_4 = pow(dt, 4);
    
    Q <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}


void Fusion::updateMeasurementModel(VectorXd x, const VectorXd z, const MeasurementPackage::SensorType type) {
    if (type == MeasurementPackage::RADAR) {
        H = jacobianHMatrix(x);
        R = radarR;
        y = radarObservationModel(x, z);
    } else {
        H = MatrixXd(2, 4);
        H << 1, 0, 0, 0,
        0, 1, 0, 0;
        R = laserR;
        y = z - H * x;
    };
};


VectorXd Fusion::radarObservationModel(const VectorXd x, const VectorXd z) {
    float posx = x(0);
    float posy = x(1);
    float velx = x(2);
    float vely = x(3);
    
    float rho = sqrt((posx * posx) + (posy * posy));
    float phi = atan2(posy, posx);
    float rho_dot = (posx * velx + posy * vely) / rho;
    
    VectorXd h(3);
    h << rho, phi, rho_dot;
    
    VectorXd y = z - h;
    while(y(1) > M_PI || y(1) < -M_PI) {
        if (y(1) > M_PI) {
            y(1) -= M_PI*2;
        } else {
            y(1) += M_PI*2;
        }
    }
    
    return y;
};


MatrixXd Fusion::jacobianHMatrix(const VectorXd z) {
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);
    
    float c1 = pow(px, 2) + pow(py, 2);
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
    
    MatrixXd HJ(3, 4);
    if(fabs(c1) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
        return HJ;
    }
    
    HJ << (px/c2), (py/c2), 0, 0,
    -(py/c1), (px/c1), 0, 0,
    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return HJ;
}

#ifndef Filter_H_
#define Filter_H_


#include "Filterable.h"
#include "measurement_package.h"


/// Fuses sensor measurements
struct Fusion: Filterable {
public:
    /// System state
    VectorXd x = VectorXd(4);
    
    /// System state covariance
    MatrixXd P = MatrixXd(4, 4);
    
    /// Process measurement package
    void process(const MeasurementPackage m);
private:
    /// Indicates Fusion has been initialized
    bool isInitialized = false;
    
    /// Last measurement timestamp
    long long savedTimestamp = 0;
    
    /// Laser sensor state and process matrices
    MatrixXd laserR = MatrixXd(2, 2);
    
    /// Radar sensor process matrices
    MatrixXd radarR = MatrixXd(3, 3);
    
    /// Initializes the sensor fusion with first measurement
    void initialize(const MeasurementPackage m);
    
    /// Updates the Kalman filter state and state covariant matrices
    void updateStateModel(long timestamp);

    /// Updates the Kalman filter measurement and measurement covariant matrices
    void updateMeasurementModel(VectorXd x, const VectorXd z, const MeasurementPackage::SensorType type);

    /// Returns observation matrix for radar measurements
    VectorXd radarObservationModel(const VectorXd x, const VectorXd z);

    /// Returns jacobian observation matrix for 3D measurements
    MatrixXd jacobianHMatrix(const VectorXd x);
};

#endif

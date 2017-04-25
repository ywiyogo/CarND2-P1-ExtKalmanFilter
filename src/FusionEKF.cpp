#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //You'll need to initialize the Kalman Filter, 
  //prepare the Q and F matrices for the prediction step, and call the radar and lidar update functions.
  F = MatrixXd(4, 4);
  Q = MatrixXd(4, 4);
  P = MatrixXd(4, 4);
  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1, 0,
       0, 0, 0, 1;
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    //ekf_.Init(VectorXd(4), MatrixXd &P_in, MatrixXd &F_in,
                        //MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
    Eigen::VectorXd x = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "\n Init measurement with RADAR \n----------------" << endl;
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = -1 * rho * sin(phi);  //-1 due to the positive y axis toward left 
      float vx = rho_dot * cos(phi);
      float vy = -1 * rho_dot * sin(phi);
      
      x << px, py, vx, vy;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "\n Init measurement with LASER \n----------------" << endl;
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    }

    // Initialize the Kalman filter
    ekf_.Init(x, P, F, H_laser_, Hj_, R_laser_, R_radar_, Q);
    //cout << "ekf_.x_ = " << ekf_.x_ << endl;
    // assign the timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  cout << "dt = " << dt << endl;
  // Set the state transition matrix
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0,  dt,
             0, 0, 1,  0,
             0, 0, 0,  1;
  // Set the covariance matrix
  ekf_.Q_ << pow(dt,4)/4*noise_ax, 0, pow(dt,3)/2*noise_ax, 0,
             0, pow(dt,4)/4*noise_ay, 0, pow(dt,3)/2*noise_ay,
             pow(dt,3)/2*noise_ax, 0, pow(dt,2)*noise_ax, 0,
             0, pow(dt,3)/2*noise_ay, 0, pow(dt,2)*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "\n RADAR update \n----------------" << endl;
    ekf_.H_radar_ = tools.CalculateJacobian(ekf_.x_);
    
    // In order to avoid a unnecesarry assignment of H and R matrix ekf_.R_ = R_radar_;
    // I separate the R and H matrix in ekf
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    cout << "\n LASER update \n----------------" << endl;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ =\n" << ekf_.x_ << endl;
  cout << "P_ =\n" << ekf_.P_ << endl;
}

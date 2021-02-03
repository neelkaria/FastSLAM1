#include <iostream>
#include <matplotlibcpp.h>
#include <cmath>
#include <Eigen/Dense>

namespace plt = matplotlibcpp;
using namespace Eigen;
using namespace std;

double toDeg = 180/M_PI, toRad = M_PI/180;

double OFFSET_YAW_RATE_NIOISE = 0.01;
double DT = 0.1;   //time tick [s]
double SIM_TIME = 50.0;  //simulation time [s]
double MAX_RANGE = 20.0;  //maximum observation range
double M_DIST_TH = 2.0;  //Threshold of Mahalanobis distance for data association.
double STATE_SIZE = 3;  //State size [x,y,yaw]
double LM_SIZE = 2;  //LM state size [x,y]
double N_PARTICLE = 100;  //number of particle
double NTH = N_PARTICLE / 1.5;  //Number of particle for re-sampling

// void update(Matrix<double, 3, 1> &x, double &y) {
//      x.a = 1.0;
//      y = 2.0;
// }

// Matrix<double, 3, 1> x;
// double y = 1.0;
// cout << y << endl; //1
// update(x, y);
// cout << y << endl; //2

class Particle
{
     Particle(int n_landmark);

     public:
     double w;
     double x;
     double y;
     double yaw;
     MatrixXd lm =MatrixXd::Zero(8, 2);
};

Particle::Particle(int n_landmark) 
     {
          this->w = 1.0 / N_PARTICLE;
          this->x = 0.0;
          this->y = 0.0;
          this->yaw = 0.0;
          this->lm;

     }

double calc_input(double time)
{    double yaw_rate,v;
     MatrixXd u(1,2);
     if (time <= 3.0)
     {
          v = 0.0;  // [m/s]
          yaw_rate = 0.0;   // [rad/s]
     }
     else
     {
          v = 1.0;  // [m/s]
          yaw_rate = 0.1; // [rad/s]
     }
     

}



int main()
{
     //Fast SLAM Covariance
     Matrix2d Q(2,2);
     Matrix2d R(2,2);

     //Simulation parameter
     Matrix2d Q_sim(2,2);
     Matrix2d R_sim(2,2);

     Q << 3, 0,
          0, 10*toRad; 
     Q=Q*Q;

     R << 1, 0,
          0, 20*toRad;
     R=R*R;

     Q_sim << 0.3, 0,
          0, 2*toRad;
     Q_sim=Q_sim*Q_sim;

     R_sim << 0.5, 0,
          0, 10*toRad;
     R_sim = R_sim*R_sim;

     double time = 0.0;

     //RFID Positions [x, y]
     Matrix <double, 8, 2> RFID;
     RFID << 10.0, -2.0,
             15.0, 10.0,
             15.0, 15.0,
             10.0, 20.0,
              3.0, 15.0,
             -5.0, 20.0,
             -5.0,  5.0,
            -10.0, 15.0;

     int n_landmark = RFID.rows();

     //State Vector [x y yaw v]
     Matrix <double, 3, 1> xEst; //SLAM Estimation
     xEst << 0,0,0;

     Matrix <double, 3, 1> xTrue; //True State
     xTrue << 0,0,0; 

     Matrix <double, 3, 1> xDr; //Dead Reckoning
     xDr << 0,0,0;   

     //history
     Matrix <double, 3, 1> hxEst = xEst;
     Matrix <double, 3, 1> hxTrue = xTrue;
     Matrix<double, 3, 1> hxDR = xTrue;  


     while (SIM_TIME >= time)
     {
          time += DT;
          auto u = calc_input(time);
     }
     

     // Cholesky inverse
}

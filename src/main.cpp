#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <random>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;
using namespace Eigen;
using namespace std;

random_device rd;  
mt19937 gen(rd()); 

double toDeg = 180/M_PI, toRad = M_PI/180;
double OFFSET_YAW_RATE_NOISE = 0.01;
double DT = 0.1;   //time tick [s]
double SIM_TIME = 50.0;  //simulation time [s]
double MAX_RANGE = 20.0;  //maximum observation range
double M_DIST_TH = 2.0;  //Threshold of Mahalanobis distance for data association.
double STATE_SIZE = 3;  //State size [x,y,yaw]
double LM_SIZE = 2;  //LM state size [x,y]
double N_PARTICLE = 100;  //number of particle
double NTH = N_PARTICLE / 1.5;  //Number of particle for re-sampling
//Fast SLAM Covariance
Matrix2d Q(2,2);
Matrix2d R(2,2);

//Simulation parameter
Matrix2d Q_sim(2,2);
Matrix2d R_sim(2,2);

class Particle
{
     public:
          Particle(int n_landmark);
          double w;
          double x;
          double y;
          double yaw;
          vector<Vector2d> lm;
          vector<Matrix2d> lmP;

};

class Jacobians
{
     public:
          Vector2d zp;
          Matrix<double, 2, 3> Hv;
          Matrix2d Hf;
          MatrixXd Sf;
};

Particle::Particle(int n_landmark) 
     {
          this->w = 1.0 / N_PARTICLE;
          this->x = 0.0;
          this->y = 0.0;
          this->yaw = 0.0;
          for(int i = 0; i < n_landmark; i++)
          {
               Matrix2d abc = Matrix2d::Zero(LM_SIZE,LM_SIZE);
               Vector2d def = Vector2d::Zero(LM_SIZE);
               lm.push_back(def);
               lmP.push_back(abc);
          }
     }

VectorXd cumsum(VectorXd &inp)
{
     VectorXd res = VectorXd::Zero(inp.size());
     res(0) = inp(0);
	for (int i = 1; i < inp.size(); i++) 
     {
          res(i) = res(i-1) + inp(i);
     }
     return res;
}


MatrixXd hstack(const MatrixXd &m1, const MatrixXd &m2)
{
    if(m1.cols() == 0){
      return m2;
    }
    else if(m2.cols() == 0){
      return m1;
    }

    uint32_t nrow = m1.rows();
    if(nrow == 0){
      nrow = m2.rows();
    }

    MatrixXd rm(nrow, m1.cols()+m2.cols());
    rm << m1, m2;
  
    return rm;
}


void normalize_weight(vector<Particle> &particles)
{
     double sum_w = 0.0;

     for (int i = 0; i <= particles.size(); i++)
     {
          sum_w += particles[i].w;
     }

     try
     {
          for(int i = 0; i< N_PARTICLE; i++)
          {
               particles[i].w /= sum_w;

          }
     }
     catch(...)
     {
          for(int i = 0; i< N_PARTICLE; i++)
          {
               particles[i].w = 1.0 / N_PARTICLE;
          }
     }
     
}


MatrixXd calc_input(double time)
{
     double yaw_rate,v;
     MatrixXd u(2,1);
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
     u << v, yaw_rate;
     return u; 
}


double pi_2_pi(double angle)
{
    return fmod(fmod((angle)+M_PI, 2*M_PI)-2*M_PI, 2*M_PI)+M_PI;
}

void add_new_landmark(Particle *particle,MatrixXd z)
{
     auto r = z(0);
     auto b = z(1);
     auto lm_id = z(2);

     double s = sin(pi_2_pi(particle->yaw + b));
     double c = cos(pi_2_pi(particle->yaw + b));
     
     particle->lm[lm_id][0] = particle->x + (r * c);
     particle->lm[lm_id][1] = particle->y + (r * s);
     
     double dx = r * c;
     double dy = r * s;
     double d2 = pow(dx,2) + pow(dy,2);
     
     double d = sqrt(d2);

     Matrix2d Gz;
     Gz <<   dx/d , dy/d,
           -dy/d2 , dx/d2;

     particle->lmP[lm_id] = Gz.inverse() * Q * Gz.transpose().inverse(); 
}


Jacobians compute_jacobians(Particle *particle, Vector2d xf, Matrix2d Pf)
{
     Jacobians j;

     double dx = xf(0) - particle->x;
     double dy = xf(1) - particle->y;
     double d2 = pow(dx,2) + pow(dy,2);
     double d = sqrt(d2);
     
     j.zp << d, pi_2_pi(atan2(dy,dx) - particle->yaw);

     j.Hv << -dx/d , -dy/d ,  0.0,
            dy/d2, -dx/d2, -1.0;
     
     j.Hf <<  dx/d , dy/d ,
           -dy/d2, dx/d2;
     
     j.Sf = j.Hf * Pf * j.Hf.transpose() + Q;
     
     return j;
}

MatrixXd motion_model(MatrixXd x, MatrixXd u)
{
     
     MatrixXd F = Matrix<double,3,3>::Identity();
     MatrixXd B(3,2);
     B << DT * cos(x(2,0)), 0,
          DT * sin(x(2,0)), 0,
          0.0, DT;
     
     x = F * x + B * u;
     x(2,0) = pi_2_pi(x(2,0));
     
     return x;
}

void predict_particles(vector<Particle> &particles,MatrixXd u, MatrixXd ud)
{
     uniform_real_distribution<> dis(-1.0, 1.0);
     MatrixXd rand_num(1,2);
     for(int i=0; i<N_PARTICLE; i++)
     {
          rand_num << dis(gen), dis(gen);
          MatrixXd px = MatrixXd::Zero(STATE_SIZE, 1);
          px(0,0) = particles[i].x;
          px(1,0) = particles[i].y;
          px(2,0) = particles[i].yaw;
          ud = u +  (rand_num * R.cwiseSqrt()).transpose();
          px = motion_model(px,ud);
          particles[i].x = px(0,0);
          particles[i].y = px(1,0);
          particles[i].yaw = px(2,0);
     }

}

double compute_weight (Particle *particle, MatrixXd z)
{
     auto lm_id = int(z(2));
     Vector2d xf = particle->lm[lm_id];
     Matrix2d Pf = particle->lmP[lm_id];

     Jacobians j = compute_jacobians(particle,xf,Pf);
     
     Vector2d dx;
     dx << z(0) - j.zp(0), z(1) - j.zp(1);
     dx(1) = pi_2_pi(dx(1));

     Matrix2d invS;
     try
     {
           invS = j.Sf.inverse();
     }
     catch(...)
     {
          return 1.0;
     }
    
     double num = exp(-0.5 * dx.transpose() * invS * dx);
     double den = 2.0 * M_PI * sqrt(j.Sf.determinant());
     double w = num / den;
     
     return w;
}

void update_kf_with_Cholesky(Jacobians &j, Particle *particle,Vector2d &xf, Matrix2d &Pf, Vector2d dz)
{
     Matrix2d S;
     auto PHt = Pf * j.Hf.transpose();
     S = (j.Hf * PHt) + Q;
     S = (S + S.transpose()) * 0.5;
     MatrixXd s_chol(S.llt().matrixL().transpose());
     auto s_chol_inv = s_chol.inverse();
     auto W1 = PHt * s_chol_inv;
     auto W = W1 * (s_chol_inv.transpose());
     xf += W * dz;
     Pf -= W1 * (W1.transpose());
}

void update_landmark(Particle *particle, MatrixXd z) 
{    
     auto lm_id = int(z(2));
     Vector2d xf = particle->lm[lm_id];
     Matrix2d Pf = particle->lmP[lm_id];
     Vector2d dz;
     
     Jacobians j = compute_jacobians(particle,xf,Pf);

     dz << z(0) - j.zp(0), z(1) - j.zp(1);
     dz(1,0) = pi_2_pi(dz(1,0));

     update_kf_with_Cholesky(j,particle,xf,Pf,dz);

     particle->lm[lm_id] = xf.transpose();
     particle->lmP[lm_id] = Pf;
}


void update_with_observation(vector<Particle> &particles, MatrixXd z)
{
     for (int iz = 0; iz < z.cols(); ++iz)
     {
          auto landmark_id = z(2,iz);
     
          for(int ip = 0; ip < N_PARTICLE; ip++)
          {
               
               if (abs(particles[ip].lm[landmark_id][0]) <= 0.01)
               {    
                     add_new_landmark(&particles[ip],z.col(iz)); 
               }
               else
               {    
                    double w = compute_weight(&particles[ip],z.col(iz));
                    particles[ip].w *= w;
                    update_landmark(&particles[ip], z.col(iz));
               }    
          }
     }
}

void resampling (vector <Particle> &particles)
{
     uniform_real_distribution<> dis(0.0, 1.0);

     normalize_weight(particles);

     VectorXd pw(particles.size());
     VectorXd w_cum(particles.size());
     for (int i = 0; i < particles.size(); i++)
     {
          pw(i) = particles[i].w;
     }

     double n_eff = 1.0 / (pw.dot(pw.transpose()));

     if(n_eff < NTH)
     {
          w_cum = cumsum(pw);
          VectorXd b = pw * 0.0 + (1/N_PARTICLE * VectorXd::Ones(pw.size()));
          VectorXd base = cumsum(b) - (1/N_PARTICLE * VectorXd::Ones(b.size()));
          VectorXd resample_id = base+ (dis(gen) * VectorXd::Ones(base.size())/N_PARTICLE);
          vector<int> inds;
          vector<Particle> temp;
          int ind = 0;
          for (int ip = 0;ip < N_PARTICLE; ip++)
          {
               while ((ind < w_cum.size() - 1) & resample_id(ip) > w_cum(ind))
               {
                    ind += 1;
               }
               inds.push_back(ind);
          }
          temp = particles;
          for (int i = 0; i < inds.size(); i++ )
          {
            particles[i].x = temp[inds[i]].x;
            particles[i].y = temp[inds[i]].y;
            particles[i].yaw = temp[inds[i]].yaw;
            particles[i].lm = temp[inds[i]].lm;
            particles[i].lmP = temp[inds[i]].lmP;
            particles[i].w = 1.0 / N_PARTICLE;
          }

     }
 
}

void fast_slam(vector<Particle> &particles, MatrixXd u, MatrixXd ud, MatrixXd z)
{
     predict_particles(particles, u, ud);
     update_with_observation(particles, z);
     resampling(particles);
}


MatrixXd calc_final_state (vector<Particle> &particles, Vector3d xEst)
{
     xEst <<0,0,0;
     normalize_weight(particles);
     
     for (int i = 0; i < N_PARTICLE; i++)
     {
          xEst(0) += particles[i].w * particles[i].x;
          xEst(1) += particles[i].w * particles[i].y;
          xEst(2) += particles[i].w * particles[i].yaw;
     }
     xEst(2) = pi_2_pi(xEst(2));
     return xEst;
}

void observation (Vector3d &xTrue, Ref<Vector3d> xd, Ref <MatrixXd> u, Ref <MatrixXd> RFID, MatrixXd& z,Vector2d& ud)
{
     uniform_real_distribution<> dis(-1.0, 1.0);
     double dx,dy,d,angle,dn;
     xTrue = motion_model(xTrue, u);
     
     for (int i=0; i<RFID.rows(); i++)
     {
          dx = RFID(i,0) - xTrue(0,0);
          dy = RFID(i,1) - xTrue(1,0);
          d = hypot(dx,dy);
          angle = pi_2_pi(atan2(dy,dx) - xTrue(2,0));
          if (d<= MAX_RANGE)
          { 
               dn = d + dis(gen) * pow(Q_sim(0,0),0.5);
               double angle_with_noise = angle + dis(gen) * pow(Q_sim(1,1),0.5);
               MatrixXd zi(3,1);
               zi << dn, pi_2_pi(angle_with_noise), i;
               z = hstack(z,zi);
          }
     }

     auto ud1 = u(0,0) + (dis(gen) * pow(R_sim(0,0),0.5));
     auto ud2 = u(1,0) + (dis(gen) * pow(R_sim(1,1),0.5)) + OFFSET_YAW_RATE_NOISE;
          
     ud << ud1,ud2;
     xd = motion_model(xd, ud);
}

int main()
{    
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

     Particle P(n_landmark);

     //State Vector [x y yaw v]
     Vector3d xEst; //SLAM Estimation
     xEst << 0,0,0;
     Vector3d xTrue; //True State
     xTrue << 0,0,0; 

     Vector3d xDR; //Dead Reckoning
     xDR << 0,0,0;   

     //history
     // MatrixXd hxEst = xEst;
     // MatrixXd hxTrue = xTrue;
     // MatrixXd hxDR = xTrue;  
     Vector3d x_state;
     vector<Particle> particles;
     for(int i=0; i<N_PARTICLE; i++)
     {
          Particle p(n_landmark);
          particles.push_back(p);
     }

     int cnt = 0;
     while (SIM_TIME >= time)
     {
          cout<<"iter no."<<cnt<<endl;
          Vector2d ud;
          MatrixXd z = MatrixXd::Zero(3,0);
         
          time += DT;
          MatrixXd u = calc_input(time);

          observation(xTrue, xDR, u, RFID, z,ud);

          fast_slam(particles, u, ud, z);
          
          xEst = calc_final_state(particles, xEst);
          
          //x_state << xEst(0), xEst(1), xEst(2);
          
          // hxEst = hstack(hxEst,x_state);
          // hxDR = hstack(hxDR,xDR);
          // hxTrue = hstack(hxTrue, xTrue);
          
          //Vectors for Plotting 
          vector<double> part_x,part_y;
          vector<double> x,y;
          vector<double> x_j,y_j;

          for(int i = 0;i < RFID.rows(); i++)
          {
               x.push_back(RFID(i,0));
               y.push_back(RFID(i,1));  
          }
          
          for(int i = 0; i < N_PARTICLE; i++)
          {
               part_x.push_back(particles[i].x);
               part_y.push_back(particles[i].y);
               for (int j = 0; j < particles[i].lm.size(); j++)
               {    
                    x_j.push_back(particles[i].lm[j][0]);
                    y_j.push_back(particles[i].lm[j][1]); 
               }
          }
          plt::plot(x_j, y_j, "xb");
          plt::plot(part_x, part_y, ".r");
          plt::plot(x,y,"xb");
          plt::axis("Equal");
          plt::grid(true);
          plt::pause(0.001);
          cnt +=1;
     }  
          
}

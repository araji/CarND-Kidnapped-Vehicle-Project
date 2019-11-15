/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // Set the number of particles (size of vector<Particle>)
  num_particles = 100 ;  
  weights.resize(num_particles);
  particles.resize(num_particles);
  // init rnd engine 
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

  // assign  random Gaussian noise to each particle , set w =1 .
  for (int i =0 ; i < num_particles ; i++) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0 ;
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
  std::default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {
    if (abs(yaw_rate) != 0 ) {
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate *delta_t;
    } else {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }

    std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (auto pred: predicted) {
      double dist_min = std::numeric_limits<double>::max();
      for(auto observation : observations){
        double distance = dist(observation.x, observation.y, pred.x, pred.y);
        if(distance < dist_min){
          observation.id = pred.id;
        }
        dist_min = distance;
      }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**

   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  double total_weight = 0 ;

  for (int i = 0; i< num_particles; i++)  {
    //needed for each particle
    vector<LandmarkObs> predictions; 
    double weight = 1.0;
    //info we need from each particle
    //double px = particles[i].x;
    //double py = particles[i].y;
    //double ptheta = particles[i].theta;

    //iterate through the observations , apply Homogeneous Transformation    
    for ( uint j = 0; j < observations.size(); j ++) {
      //observation from vehicle in map coordinates
      double map_obs_x, map_obs_y ;
      map_obs_x = particles[i].x + (observations[j].x * cos(particles[i].theta)) - (observations[j].y * sin(particles[i].theta) ) ;
      map_obs_y = particles[i].y + (observations[j].x * sin(particles[i].theta)) - (observations[j].y * cos(particles[i].theta) );

      //initialize threshold to largest value possible before we start iterating
      double distance_min = std::numeric_limits<double>::max(); 
      // will hold the closest landmark based on all observations
      Map::single_landmark_s landmark ;

      for (uint k =0 ; k < map_landmarks.landmark_list.size() ; ++k ) {
        Map::single_landmark_s current_landmark = map_landmarks.landmark_list[k];
        double distance = dist(map_obs_x, map_obs_y,current_landmark.x_f, current_landmark.y_f);
        if (distance < distance_min) {
          distance_min = distance;
          landmark = current_landmark ;
        }
      }

    /**
    N = 2 * pi * sigmax * sigmay
    X =   pow((x - mux),2) / ( 2* sigmax *sigmax )
    Y =   pow((y - muy),2) / ( 2* sigmay *sigmay )
    return (1/N )* exp( - (X+Y) )
    **/
      double X = pow((map_obs_x - landmark.x_f), 2) / pow(std_landmark[0], 2) ;
      double Y = pow((map_obs_y - landmark.y_f), 2) / pow(std_landmark[1], 2) ;
      double N = 2 * M_PI * std_landmark[0]  * std_landmark[1] ;
      weight *= (1/N) * exp( - (X + Y) ) ;

    }
    total_weight += weight ;
    particles[i].weight = weight;
  }

  for (int i = 0; i< num_particles; i++)  {
    particles[i].weight /= total_weight;
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   **/
   
  std::default_random_engine gen;
  std::discrete_distribution<int> dist(weights.begin(), weights.end()); 

  vector<Particle> resampled_particles;
  for (int i = 0; i < num_particles; i++){
        resampled_particles.push_back(particles[dist(gen)]);
  }
  particles = resampled_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
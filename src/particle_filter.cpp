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

//macros and globals
#define MIN_YAW 0.00001
std::default_random_engine gen;


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
  weights.resize(num_particles, 1.0);

  //normal distribution centered around values
  std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

  // assign  random Gaussian noise to each particle , set w =1 .
  for (int i =0 ; i < num_particles ; i++) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0 ;
    particles.push_back(p);
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

  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {
  
    //avoid dividing by zero
    if (abs(yaw_rate) > MIN_YAW  ) {
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate *delta_t;
    } else {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }

    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
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
  for (int i = 0; i < observations.size(); i++)
    {
        LandmarkObs objectO = observations[i];
        double min_dist = std::numeric_limits<double>::max();
        int map_id = -1;
        for (int j = 0; j < predicted.size(); j++)
        {
            LandmarkObs objectP = predicted[j];
            double cur_dist = dist(objectO.x, objectO.y, objectP.x, objectP.y);
            if (cur_dist < min_dist)
            {
                min_dist = cur_dist;
                map_id = objectP.id;
            }
        }
        observations[i].id = map_id;
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
  
  double max_reach = sensor_range*sensor_range ;
  double multivariate_term1 = 1 / ( 2 * M_PI * std_landmark[0] * std_landmark[1]);
  
  for (int i = 0; i< num_particles; i++)  {
    
    vector<LandmarkObs> landmarksInRange; 
    //double weight = 1.0;
    double px = particles[i].x;
    double py = particles[i].y;
    double ptheta = particles[i].theta;
    //keep landmarks that are within sensor_radar range , push into landmarksInRange
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float lx = map_landmarks.landmark_list[j].x_f;
			float ly = map_landmarks.landmark_list[j].y_f;
			int lid = map_landmarks.landmark_list[j].id_i;
			double dx = px - lx;
			double dy = py - ly;

			if (dx*dx + dy*dy <= max_reach ) {
				landmarksInRange.push_back(LandmarkObs{ lid, lx, ly });
			}
		}
    
    //std::cout << "landmarks in sensor range  " << landmarksInRange.size() << std::endl;

    //iterate through the observations , apply homogeneous transformations push result into transformedObservations
    vector<LandmarkObs> transformedObservations;
    for ( uint j = 0; j < observations.size(); j ++) {
      double mox, moy ;
      mox = px + (observations[j].x * cos(ptheta)) - (observations[j].y * sin(ptheta) ) ;
      moy = py + (observations[j].x * sin(ptheta)) + (observations[j].y * cos(ptheta) );
      transformedObservations.push_back(LandmarkObs{ observations[j].id, mox, moy });
    }
    //update observations with closest landmarks in range using helper 
    dataAssociation(landmarksInRange, transformedObservations);

    //reset and compute weigths
    particles[i].weight = 1.0;

    for (uint j = 0; j < transformedObservations.size(); j++) {
			double transformedObsX  = transformedObservations[j].x;
			double transformedObsY  = transformedObservations[j].y;
			int    TransformedObsId = transformedObservations[j].id;
      
      double landmarkX, landmarkY ;
      int landmarksInRangeSize = landmarksInRange.size();

      // 
      int k = 0;
			bool found = false;
			while (!found && k < landmarksInRangeSize) {
				if (landmarksInRange[k].id == TransformedObsId) {
					found = true;
					landmarkX = landmarksInRange[k].x;
					landmarkY = landmarksInRange[k].y;
				}
				k++;
			}

      // get multivariate gaussian 
      //dx == x - MUx ; dy = y-MUy
      double dX = transformedObsX - landmarkX;
			double dY = transformedObsY - landmarkY;
      //std_landmark[0] * std_landmark[1])
      double multivariate_term2X =  dX * dX / (2 * std_landmark[0] * std_landmark[0]) ;
      double multivariate_term2Y =  dY * dY / (2 * std_landmark[1] * std_landmark[1]) ;
      double multivariate_term2 = exp( - ( multivariate_term2X + multivariate_term2Y ) );
      double weight = multivariate_term1 * multivariate_term2 ;
			particles[i].weight = particles[i].weight * weight;
      weights[i] = particles[i].weight;
		}
  }

}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   **/

  std::discrete_distribution<int> weighted_dist(weights.begin(), weights.end());

  std::vector<Particle> resampled_particles;

  for (int i = 0; i < num_particles; ++i) {
    int j = weighted_dist(gen);
    resampled_particles.push_back(particles[j]);
  }
  // exchange old particle list with new resampled list
  particles = std::move(resampled_particles);

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
# Kidnapped Car Project 

Implementation of 2D particle filter as part of the Udacity Self-Driving Car Nanodegree .

## Problem statement :

Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.


## Implementation details:

Overall Flow shown here (Source = Udacity ) :

![ParticleFilter Process Flow] (https://https://github.com/araji/CarND-Kidnapped-Vehicle-Project/udacity_car_nd_localization_process.png)

Given a map with a set of known landmarks (global coordinates ) ,  Particle filter consists of populating the map with a number of particles that could all be the position of the robot .
Through the robot mouvements and sensor readings we update those beliefs narrow the choices for each particles by removing any ladmarks that are out of it range they weigh each particle based on how close that measurement to the belief assocaited wit the particle .

The resampling phase draws from that set of particles (with replacement ) based on the particle weight and therefore converging towards the robot likely postion.



## Running the Code
1. ./clean.sh
2. ./build.sh
3. ./run.sh


![Kidnapped Car Project] (https://github.com/araji/CarND-Kidnapped-Vehicle-Project/particle_filter_at_Work.png)


/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace std;

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 5;

	// noise generationa
	default_random_engine gen;
	double std_x = std[0];  // These are sigma not std
	double std_y = std[1];  // These are sigma not std
  double std_theta = std[2];  // These are sigma not std

	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; i++)
	{
		Particle new_particle;
		new_particle.x = dist_x(gen);
		new_particle.y = dist_y(gen);
		new_particle.theta = dist_theta(gen);
		new_particle.weight = 1.0;

		particles.push_back(new_particle);
	}

	is_initialized = true;

	cout << "Vector size:" << particles.size() << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// noise generation
	default_random_engine gen;
	double std_x = std_pos[0];  // These are sigma not std
	double std_y = std_pos[1];  // These are sigma not std
	double std_theta = std_pos[2];  // These are sigma not std

	std::normal_distribution<double> dist_x(0, std_x);
	std::normal_distribution<double> dist_y(0, std_y);
	std::normal_distribution<double> dist_theta(0, std_theta);

	// Motion equations
	for (int i = 0; i < num_particles; i++)
	{
		if (fabs(yaw_rate) > 1e-5)
		{
			particles[i].x = particles[i].x + (velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)));
			particles[i].y = particles[i].y + (velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)));
			particles[i].theta = particles[i].theta + yaw_rate * delta_t;
		} else {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}

		// add gaussian noise
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	for (int part_index=0; part_index < particles.size(); part_index++)
	{	
		double part_x = particles[part_index].x;
		double part_y = particles[part_index].y;
		double part_theta = particles[part_index].theta;
		particles[part_index].weight = 1.0;

		for (int obs_index=0; obs_index < observations.size(); obs_index++)
		{
			double obs_x = observations[obs_index].x;
			double obs_y = observations[obs_index].y;

			// 1) Transform Coordinates
			// transform each observation from car's coordinates to map coordinates
			double obs_x_map = part_x + obs_x * cos(part_theta) - obs_y * sin(part_theta);
			double obs_y_map = part_y + obs_x * sin(part_theta) + obs_y * cos(part_theta);

			// 2) Associate with the nearest landmark
			double nearest_lm_dist = 2.0 * sensor_range;
			int nearest_lm_id = 0;
			int nearest_lm_index = 0;

			for (int lm_index = 0; lm_index < map_landmarks.landmark_list.size(); lm_index++)
			{
				int lm_id = map_landmarks.landmark_list[lm_index].id_i;
				double lm_x = map_landmarks.landmark_list[lm_index].x_f;
				double lm_y = map_landmarks.landmark_list[lm_index].y_f;

				double diff_x_sq = (lm_x - obs_x_map) * (lm_x - obs_x_map);
				double diff_y_sq = (lm_y - obs_y_map) * (lm_y - obs_y_map);

				double landmark_dist = sqrt(diff_x_sq + diff_y_sq);

				if (landmark_dist < nearest_lm_dist)
				{
					nearest_lm_dist = landmark_dist;
					nearest_lm_id = lm_id;
					nearest_lm_index = lm_index;
					//cout << lm_index << ": nearest_lm_dist:" << nearest_lm_dist << ", nearest_lm_index:" << nearest_lm_index << ", nearest_lm_id:" << nearest_lm_id <<  endl;
				}	
				else
				{
					//cout << "\t" << lm_index << ": nearest_lm_dist:" << nearest_lm_dist << ", nearest_lm_index:" << nearest_lm_index << ", nearest_lm_id:" << nearest_lm_id <<  endl;					
				}	
			}

			observations[obs_index].id = nearest_lm_id;

			// 3) Update Weights
			double x_mean = map_landmarks.landmark_list[nearest_lm_index].x_f;
			double y_mean = map_landmarks.landmark_list[nearest_lm_index].y_f;
			double x = obs_x_map;
			double y = obs_y_map;
			//cout << "Obs_index:" << obs_index << ", Landmark_id:" << nearest_lm_id << ", X_mean:" << x_mean << ", X:" << x << ", Y_mean:" << y_mean << ", Y:" << y << endl;
			double x_diff_sq = (x - x_mean)*(x - x_mean);
			double y_diff_sq = (y - y_mean)*(y - y_mean);
			double sigma_landmark_x = std_landmark[0];
			double sigma_landmark_y = std_landmark[1];
			double numerator = exp(-1.0 * (x_diff_sq/(2.0 * sigma_landmark_x * sigma_landmark_x) + y_diff_sq/(2.0 * sigma_landmark_y * sigma_landmark_y)));
			double denominator = 1.0/2.0 * M_PI * sigma_landmark_x * sigma_landmark_y;
	
			particles[part_index].weight *= numerator*denominator;
			//cout << "\tParticle:" << part_index << ", numerator:" << numerator<< ", denominator:" << denominator << ", weight:" << particles[part_index].weight << endl;
		
			//cout << "Landmark:" << nearest_lm_id << " x:" << x_mean << " y:" << y_mean << endl;
			//cout << "Observation_id:" << observations[obs_index].id << " x:" << x << " y:" << y << endl;

		}
		cout << "Particle:" << part_index << ", weight:" << particles[part_index].weight << endl;

	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<Particle> particles_new;

	for (int i=0; i < particles.size(); i++)
	{
		weights.push_back(particles[i].weight);
	}

  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  //std::map<int, int> m;
  for(int n=0; n<num_particles; ++n) {
	  particles_new.push_back(particles[d(gen)]);
  }

  // what about the old memory of particles that didn't get pick??
  // don't we need to free them?

  /*
  for (int i=0; i < particles.size(); i++)
  {
  	cout << "Old Particles:" << particles[i].x << "," << particles[i].y << "," << particles[i].theta << ", weights:" << particles[i].weight << endl;
  	cout << "New Particles:" << particles_new[i].x << "," << particles_new[i].y << "," << particles_new[i].theta << ", weights:" << particles_new[i].weight << endl;
  }
  */
  particles = particles_new;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

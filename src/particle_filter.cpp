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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles=20;
	default_random_engine gen;
	
	normal_distribution<double> dist_x(0,std[0]);
	normal_distribution<double> dist_y(0,std[1]);
	normal_distribution<double> dist_theta(0,std[2]);
	for(int i=0;i<num_particles;i++){
		Particle particle;
		particle.id=i;
		particle.x=x+dist_x(gen);
		particle.y=y+dist_y(gen);
		particle.theta=theta+dist_theta(gen);
		particles.push_back(particle);
		particle.weight=1;
		weights.push_back(1);
	}
	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	fill(weights.begin(),weights.end(),1);
	for(int i=0;i<num_particles;i++){
		
		
		particles[i].weight=1;
		normal_distribution<double> dist_x(0,std_pos[0]);
		normal_distribution<double> dist_y(0,std_pos[1]);
		normal_distribution<double> dist_theta(0,std_pos[2]);
		if(fabs(yaw_rate)<0.001){
			particles[i].x=	particles[i].x+velocity*delta_t*cos(particles[i].theta)+dist_x(gen);
			particles[i].y=	particles[i].y+velocity*delta_t*sin(particles[i].theta)+dist_y(gen);
		}
		else{
		particles[i].x=particles[i].x+(velocity/yaw_rate)*(sin(particles[i].theta+delta_t*yaw_rate)-sin(particles[i].theta))+dist_x(gen);
		particles[i].y=particles[i].y+(velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+delta_t*yaw_rate))+dist_y(gen);
		particles[i].theta=particles[i].theta+yaw_rate*delta_t+dist_theta(gen);
	}		
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for(int i=0;i<num_particles;i++)
	{
		
		double x=particles[i].x;
		double y=particles[i].y;
		double theta=particles[i].theta;
		vector<LandmarkObs> map_cordinates_obs;
		for(unsigned int j=0;j<observations.size();j++)
		{
			LandmarkObs obs;
			obs.x=x+cos(theta)*observations[j].x-sin(theta)*observations[j].y;
			obs.y=y+sin(theta)*observations[j].x+cos(theta)*observations[j].y;
			obs.id=-1;
			map_cordinates_obs.push_back(obs);
		}
		
		for(unsigned int j=0;j<map_cordinates_obs.size();j++)
		{
			double obs_x=map_cordinates_obs[j].x;
			double obs_y=map_cordinates_obs[j].y;
			double min_distance=-1;
			int landmark_id_temp=-1;
			int index_in_list=-1;
			for(unsigned int k=0;k<map_landmarks.landmark_list.size();k++)
			{
				double landmark_x=map_landmarks.landmark_list[k].x_f;
				double landmark_y=map_landmarks.landmark_list[k].y_f;
				int landmark_id=map_landmarks.landmark_list[k].id_i;
				
				double dist_from_car=dist(x,y,landmark_x,landmark_y);
				if(dist_from_car>sensor_range) continue;

				double difference=dist(obs_x,obs_y,landmark_x,landmark_y);
				
				if(difference<min_distance)
				{
					min_distance=difference;
					landmark_id_temp=landmark_id;
					index_in_list=k;
				}
				if(landmark_id_temp==-1)
				{
					min_distance=difference;
					landmark_id_temp=landmark_id;
					index_in_list=k;
				}
			}

			map_cordinates_obs[j].id=landmark_id_temp;

			double mu_x=obs_x;
			double mu_y=obs_y;
			double landmark_x=map_landmarks.landmark_list[index_in_list].x_f;
			double landmark_y=map_landmarks.landmark_list[index_in_list].y_f;
			
			
			double gauss_norm=(1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
			double exponent=(pow((landmark_x - mu_x),2))/(2 * std_landmark[0]*std_landmark[0]) + (pow((landmark_y - mu_y),2))/(2 * std_landmark[1]*std_landmark[1]);
	
			weights[i]*=gauss_norm*exp(-exponent);
			 
		}
		particles[i].weight=weights[i];
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(),weights.end());
	vector<Particle> resample_particles;
	for( int i=0;i<num_particles;i++){
		int index=distribution(gen);
		resample_particles.push_back(particles[index]);
	}
	particles=resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
	
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

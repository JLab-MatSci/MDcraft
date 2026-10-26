#pragma once
#include <random>
#include <mutex>
#include <mdcraft/data/atom.h>
#include <mdcraft/solver/thermostat/base.h>

namespace mdcraft { 
namespace solver { 
namespace thermostat {

using ::mdcraft::data::Atoms;
using ::mdcraft::data::vector;

class Langevin : public Base {
public:
	Langevin(
		double beta,
		double temperature,
		double time_step,
		int    heat_x = 1,
		int    heat_y = 1,
		int    heat_z = 1,
		double Ux = 0.0,
		double Uy = 0.0,
		double Uz = 0.0,
		double xmin = -1e+15,
		double xmax =  1e+15
	);

	void apply_one(Atoms::iterator atom);

	void set_friction(double beta);
	void set_temperature(double T);
	void set_average_velocity(vector U);
	void set_heating_axes(
		int heat_x = 1,
		int heat_y = 1,
		int heat_z = 1
	);
	void set_dimensions(double xmin, double xmax);

	double friction();
	double temperature();
	vector average_velocity();
	double xmin();
	double xmax();

private:
	std::mutex state_mutex;

	double sigma, 
	       beta,
	       dt, 
	       T; 

	double bxmin, bxmax;

	vector U, heat_axis;

	std::random_device random_dev;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution;

};

} // thermostat
} // solver
} // mdcraft
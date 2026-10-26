#include <mdcraft/solver/thermostat/langevin.h>
#include <stdexcept>

namespace mdcraft {
namespace solver {
namespace thermostat {

Langevin::Langevin(
	double beta,
	double temperature,
	double time_step,
	int    heat_x,
	int    heat_y,
	int    heat_z,
	double Ux,
	double Uy,
	double Uz,
	double xmin,
	double xmax
) : Base(),
    beta(beta), T(temperature), dt(time_step), bxmin(xmin), bxmax(xmax),
    heat_axis(heat_x, heat_y, heat_z),
    U(Ux, Uy, Uz)
{
    if (!(beta >= 0) || !(T >= 0) || !(dt > 0) ||
        !std::isfinite(beta) || !std::isfinite(T) || !std::isfinite(dt))
        throw std::invalid_argument("Langevin requires finite beta/T >= 0 and dt > 0");
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, 1.0);
}

void Langevin::apply_one(Atoms::iterator atomit) {
    std::lock_guard<std::mutex> lock(state_mutex);
	auto& atom = *atomit;
	if ((atom.r.x() < bxmin) || (atom.r.x() > bxmax)) return;
	vector random_force(
		distribution(generator), 
		distribution(generator), 
		distribution(generator)
	);
	const double sqrt_mass = std::sqrt(atom.m);
	vector therm_force = random_force*(sigma*sqrt_mass);
	therm_force -= (atom.v - U)*beta*atom.m;
	therm_force = therm_force.array() * heat_axis.array();
	atom.f += therm_force;
}

void Langevin::set_average_velocity(vector _U) {
    std::lock_guard<std::mutex> lock(state_mutex); U = _U; }

void Langevin::set_friction(double _beta) {
    std::lock_guard<std::mutex> lock(state_mutex);
    if (!(_beta >= 0) || !std::isfinite(_beta)) throw std::invalid_argument("Invalid friction");
	beta = _beta;
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, 1.0);
}

void Langevin::set_temperature(double _T) {
    std::lock_guard<std::mutex> lock(state_mutex);
    if (!(_T >= 0) || !std::isfinite(_T)) throw std::invalid_argument("Invalid temperature");
	T = _T;
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, 1.0);
}

void Langevin::set_heating_axes(
	int heat_x,
	int heat_y,
	int heat_z
) {
    std::lock_guard<std::mutex> lock(state_mutex);
	heat_axis = vector(heat_x, heat_y, heat_z);
}

void Langevin::set_dimensions(double _xmin, double _xmax) {
    std::lock_guard<std::mutex> lock(state_mutex);
	bxmin = _xmin;
	bxmax = _xmax;
}

double Langevin::friction() { std::lock_guard<std::mutex> lock(state_mutex); return beta;}
double Langevin::temperature() { std::lock_guard<std::mutex> lock(state_mutex); return T;}
vector Langevin::average_velocity() { std::lock_guard<std::mutex> lock(state_mutex); return U;}
double Langevin::xmin() { std::lock_guard<std::mutex> lock(state_mutex); return bxmin; }
double Langevin::xmax() { std::lock_guard<std::mutex> lock(state_mutex); return bxmax; }

} // thermostat
} // solver
} // mdcraft
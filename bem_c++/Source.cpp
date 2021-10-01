#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"

using namespace std;
std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;


std::vector<double> rod(double kappa, size_t num_points,
	const std::function<double(double)>& mu,
	const std::function<double(double)>& rho)
{
	const std::vector<std::function<double(double, const std::vector<double>&)>> equation = {
			[=](double x, const std::vector<double>& y) {return y[1] / mu(x); },
			[=](double x, const std::vector<double>& y) {return -kappa * kappa * rho(x) * y[0]; },
	};
	const boundary_value_problem<double> bvp = { equation, {{0,0.0}}, {{1,1.0}} };
	vector<double> points(num_points);
	const double h = 1.0 / num_points;
	for (size_t i = 0; i < num_points; i++)
	{
		points[i] = (i + 1) * h;
	}
	auto sol = bvp.solve(points);
	std::vector<double> solution;
	solution.reserve(num_points);
	std::for_each(sol.begin(), sol.end(), [&](const auto& x) {solution.push_back(x[0]); });
	return solution;
}

double amplitude(double kappa, const std::function<double(double)>& mu,	
	const std::function<double(double)>& rho)
{
	const std::vector<std::function<double(double, const std::vector<double>&)>> equation = {
			[=](double x, const std::vector<double>& y) {return y[1] / mu(x); },
			[=](double x, const std::vector<double>& y) {return -kappa * kappa * rho(x) * y[0]; },
	};
	const boundary_value_problem<double> bvp = { equation, {{0,0.0}}, {{1,1.0}} };
	return bvp.solve()[0];
}


std::map<double, double> frequency_response(double max_kappa, size_t num_points, const std::function<double(double)> & func)
{
	std::map<double, double> result;
	const auto h = max_kappa / (num_points - 1);
	for (double kappa = 0; kappa <= max_kappa; kappa+=h)
	{
		result[kappa] = func(kappa);
	}
	return result;
}


int main()
{
	const size_t num_points = 40;
	const double h = 1.0 / num_points;
	const auto mu = [](double x) {return 1 + x; };
	const auto rho = [](double x) {return 1.0; };
	//plotTheWaveField({ {"blue", rod(10.0,num_points, [](double x) {return 1 + x; }, [](double x) {return 1.0; })} }, "1.txt", h);
	plot_the_frequency_response(frequency_response(5.0, 100, [=](double x) {return amplitude(x, mu, rho); }), "1.txt");
	system("pause");
}

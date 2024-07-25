/*
    Monte Carlo American Option Pricing
    Roger Li

    Steps (Hull 27.8):
        1. Generate GBM paths.
        2. Get cashflow of the last exercise point (at t_n) in its intrinsic value;
        3. Get the discounted values from 2, they will be V's; and the prices that are
            in the money at t_n-1, they will be S's for linear regression for
                V = a + bS + cS^2
            or, less accurately but more efficiently:
                V = a + bS
        4. Perform regression on the data to get a, b, and c.
        5. With the newly obtained best-fit relationship where V is the continuation
            value and S is the stock price at t_n-1...
        6. Get cashflow at t_n-1 by comparing the continuation value and the immediate
            exercise value. If the immediate exercise value is lower than the
            continuation value, cashflow would be 0.0, otherwise it would be the
            immediate exercise value, and the option's valuation in that simluated path
            would be the cashflow discounted at the time the option is exercised!
        7. Repeat for t_n-2, t_n-3 ... until t = 0;
        8. Discount all cashflow in all paths, get average.
*/

#include <iostream>
#include <cmath>
#include <random>

void simulate_GBM_path(double S0, double mu, double sigma, double T, int N, std::vector<double>& path)
{
    std::random_device rd;
    std::mt19937 gen(rd());     // random seed.

    std::normal_distribution<> d(0.0, 1.0); // normal distribution

    double dt = T / N;

    path[0] = S0;

    // Generate GBM path
    for (int i = 1; i < path.size(); ++i) {
        double z_i = d(gen);  // random number from normal distribution.

        // GBM formula (Hull 21.16):
        path[i] = path[i-1] * exp( (mu - sigma*sigma*0.5)*dt + (sigma*z_i*sqrt(dt)) );
    }
}

double discount_price(double S, double r, double t)
{
    double discounted_price = S * exp( -1*r*t );
    return discounted_price;
}

void simulate_price_paths(double S0, double mu, double sigma, double T, 
    int N, double K, bool is_call, std::vector<std::vector<double>>& paths)
{
    for (auto rowIt = paths.begin(); rowIt != paths.end(); ++rowIt)
    {
        simulate_GBM_path(S0, mu, sigma, T, N, *rowIt);
    }
}

void calculate_intrinsic_values(double K, bool is_call, int exercise_point,
    std::vector<std::vector<double>>& price_matrix, std::vector<std::vector<double>>& cf_matrix)
{
    for (int cur_row = 0; cur_row < price_matrix.size(); cur_row++)
    {
        // call: S - K; put: K - S
        cf_matrix[cur_row][exercise_point] = 
            std::max(price_matrix[cur_row][exercise_point] - K, 0.0) ? is_call :
            std::max(K - price_matrix[cur_row][exercise_point], 0.0);
    }
}

void lsmc_american_option_pricing(double S0, double mu, double sigma, double T, 
    int N, double K, bool is_call, int num_paths)
{
    // initialize matrices.
    std::vector<std::vector<double>> simulated_price_paths(num_paths, std::vector<double>(N+1, 0.0));
    std::vector<std::vector<double>> cash_flow_paths(num_paths, std::vector<double>(N+1, 0.0));

    // perform GBM simulation for underlying asset price.
    simulate_price_paths(S0, mu, sigma, T, N, K, is_call, simulated_price_paths);

    // 
}

int main()
{
    lsmc_american_option_pricing(
        100.00, 
        -0.075/365, 
        0.2/365, 
        30, 
        30, 
        100.00, 
        true,
        5
    );

    // for (auto rowIt = paths.begin(); rowIt != paths.end(); ++rowIt) {
    //     for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
    //         std::cout << *colIt << " "; // Dereferencing colIt to access the value
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}
/*
    Monte Carlo American Option Pricing Calculator
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

// emcc monte_carlo.cpp -o lsmc_calculator.js -s EXPORTED_FUNCTIONS='["_lsmc_american_option_pricing_WASM"]' -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' -s ALLOW_MEMORY_GROWTH=1

#include <iostream>
#include <cmath>
#include <random>
#include <numeric>

extern "C" {
    double lsmc_american_option_pricing(double S0, double mu, double sigma, double T, 
        int N, double K, bool is_call, int num_paths, bool is_european, bool reduce_variance);

    double lsmc_american_option_pricing_WASM(
        double S0, double mu, double sigma, double T, 
        int N, double K, bool is_call, int num_paths, bool is_european, bool reduce_variance
    ) {
        return lsmc_american_option_pricing(S0, mu, sigma, T, N, K, 
            is_call, num_paths, is_european, reduce_variance);
    }
}

void simulate_GBM_path(double S0, double mu, double sigma, 
    double T, int N, int path_index,
    std::vector<std::vector<double>>& paths, int true_path_nums, bool reduce_variance)
{
    std::random_device rd;
    std::mt19937 gen(rd());     // random seed.

    std::normal_distribution<> d(0.0, 1.0); // normal distribution

    double dt = T / N;

    std::vector<double>& path = paths[path_index];

    path[0] = S0;
    if (reduce_variance) paths[path_index + true_path_nums][0] = S0;

    // Generate GBM path
    for (int i = 1; i < path.size(); ++i) {
        double z_i = d(gen);  // random sample from normal distribution.

        // GBM formula (Hull 21.16):
        double deltaWt = (sigma*z_i*sqrt(dt));
        path[i] = path[i-1] * exp( (mu - sigma*sigma*0.5)*dt + deltaWt );
        if (reduce_variance)
        {
            std::vector<double>& antithetic_path = paths[path_index + true_path_nums];
            antithetic_path[i] = antithetic_path[i-1] * exp( (mu - sigma*sigma*0.5)*dt - deltaWt );
        }
    }
}

double discount_price(double S, double r, double t)
{
    return S * exp( -1*r*t );   // t should be in days. Assumes continuous compounding.
}

void simulate_price_paths(double S0, double mu, double sigma, double T, 
    int N, double K, std::vector<std::vector<double>>& paths, bool reduce_variance)
{
    int true_path_nums = reduce_variance ? paths.size()/2 : paths.size();
    for (int path = 0; path < true_path_nums; path++)
    {
        simulate_GBM_path(
            S0,
            mu,
            sigma,
            T,
            N,
            path,
            paths,
            true_path_nums,
            reduce_variance
        );
    }
}

bool get_immediate_payoff_samples(double K, bool is_call, int exercise_point,
    std::vector<std::vector<double>>& price_matrix, std::vector<double>& optimal_cfs)
{
    bool has_in_the_money_paths = false;
    for (int cur_row = 0; cur_row < price_matrix.size(); cur_row++)
    {
        double payoff = 
            is_call ? 
            std::max(price_matrix[cur_row][exercise_point] - K, 0.0) :
            std::max(K - price_matrix[cur_row][exercise_point], 0.0);
        optimal_cfs[cur_row] = payoff;
        if (payoff > 0.0) has_in_the_money_paths = true;
    }
    return has_in_the_money_paths;
}


void perform_linear_regression(std::vector<double>& x, 
    std::vector<double>& y, std::vector<double>& results)
{
    int n = x.size();

    double sumX = std::accumulate(x.begin(), x.end(), 0.0);
    double sumY = std::accumulate(y.begin(), y.end(), 0.0);
    double sumXY = 0.0;
    double sumX2 = 0.0;

    for (int i = 0; i < n; ++i) {
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double denominator = (n * sumX2 - sumX * sumX);
    if (denominator == 0) {
        std::cerr << "warning: trying to divide by zero." << std::endl;
        results[0] = 1;
        results[1] = 0;
        return;
    }

    double m = (n * sumXY - sumX * sumY) / denominator;
    double b = (sumY * sumX2 - sumX * sumXY) / denominator;

    results[0] = m;
    results[1] = b;
}

bool get_continuation_samples(std::vector<int>& optimal_exercise_pt,
    std::vector<double>& optimal_exercise_payoff, double r, int cur_t, double dt, 
    std::vector<double>& discounted_cont_val_samples)
{
    bool has_continuation_values = false;
    for (int cur_path = 0; cur_path < discounted_cont_val_samples.size(); cur_path++)
    {
        if (optimal_exercise_pt[cur_path] == -1)
        {
            discounted_cont_val_samples[cur_path] = 0.0;
            continue;
        }
        discounted_cont_val_samples[cur_path] = 
            discount_price(optimal_exercise_payoff[cur_path],
            r,
            (optimal_exercise_pt[cur_path] - cur_t) * dt
        );
        if (discounted_cont_val_samples[cur_path] > 0.0) {
            has_continuation_values = true;
        }
    }
    return has_continuation_values;
}

double get_estimated_continuation_value(std::vector<double>& regression_results, double x)
{
    return regression_results[0] * x + regression_results[1];
}

void update_optimal_values(int cur_exercise_pt, std::vector<double>& regression_results, 
    std::vector<int>& optimal_exercise_pt, std::vector<double>& optimal_exercise_payoff,
    std::vector<double>& immediate_payoff)
{
    for (int path = 0; path < immediate_payoff.size(); path++)
    {
        if (immediate_payoff[path] <= 0.0) continue;
        double estimated_continuation_value = get_estimated_continuation_value(
            regression_results, immediate_payoff[path]
        );
        if (estimated_continuation_value < immediate_payoff[path])
        {
            optimal_exercise_pt[path] = cur_exercise_pt;
            optimal_exercise_payoff[path] = immediate_payoff[path];
        }
    }
}

double lsmc_american_option_pricing(double S0, double mu, double sigma, double T, 
    int N, double K, bool is_call, int num_paths, bool is_european=false, 
    bool reduce_variance=false)
{

    std::cout << "New calculation:" << std::endl;
    std::cout << "\t S0\t" << S0 << std::endl;
    std::cout << "\t mu\t" << mu << std::endl;
    std::cout << "\t sigma\t" << sigma << std::endl;
    std::cout << "\t T\t" << T << std::endl;
    std::cout << "\t N\t" << N << std::endl;
    std::cout << "\t K\t" << K << std::endl;
    std::cout << "\t is_call\t" << is_call << std::endl;
    std::cout << "\t num_paths\t" << num_paths << std::endl;
    std::cout << "\t is_european\t" << is_european << std::endl;
    std::cout << "\t reduce_variance\t" << reduce_variance << std::endl;
    std::cout << "Starting calculation..." << std::endl;

    // daily-fy data:
    mu = mu / 365;
    sigma = sigma / sqrt(365);

    // initialize matrices.
    if (reduce_variance) num_paths = num_paths * 2;
    std::vector<std::vector<double>> simulated_price_paths(num_paths, std::vector<double>(N+1, 0.0));
    std::vector<int> optimal_exercise_pt(num_paths, -1);
    std::vector<double> optimal_exercise_payoff(num_paths);

    std::cout << "Performing GBM price simulations..." << std::endl;

    // perform GBM simulation for underlying asset price.
    simulate_price_paths(S0, mu, sigma, T, N, K, simulated_price_paths, reduce_variance);

    std::cout << "Performing backward induction to determine optimal exercise points..." << std::endl;

    // handle first iteration.
    get_immediate_payoff_samples(K, is_call, N,
        simulated_price_paths, optimal_exercise_payoff);
    for (int path = 0; path < optimal_exercise_payoff.size(); path++) 
    {
        if (optimal_exercise_payoff[path] > 0.0)   
            optimal_exercise_pt[path] = N;
    }

    // perform the main loop for LSMC
    if (!is_european)
    {
        std::vector<double> immediate_payoff_samples(num_paths);
        std::vector<double> discounted_continuation_value_samples(num_paths);
        std::vector<double> regression_results(2);
        for (int cur_exercise_pt = simulated_price_paths[0].size()-2; 
            cur_exercise_pt > -1; cur_exercise_pt--)
        {
            if (
                !get_immediate_payoff_samples(
                    K, is_call, cur_exercise_pt, 
                    simulated_price_paths, immediate_payoff_samples
                )
            ) continue;
            get_continuation_samples(
                optimal_exercise_pt, optimal_exercise_payoff, 
                mu, cur_exercise_pt, T/N, 
                discounted_continuation_value_samples
            );
            perform_linear_regression(immediate_payoff_samples, 
                discounted_continuation_value_samples, regression_results);
            update_optimal_values(cur_exercise_pt, regression_results, 
                optimal_exercise_pt, optimal_exercise_payoff, immediate_payoff_samples);
        }
    }

    std::cout << "Discounting and averaging optimal payoffs..." << std::endl;

    // add and average.
    double sum = 0.0;
    for (int path = 0; path < optimal_exercise_pt.size(); path++)
    {
        if (optimal_exercise_pt[path] == -1) continue;
        sum += discount_price(optimal_exercise_payoff[path],
            mu, optimal_exercise_pt[path] * (T/N));
    }

    double estimated_price = sum / num_paths;

    std::cout << "Final estimated price: " << estimated_price << std::endl;

    return estimated_price;
}
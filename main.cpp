#include <iostream>
#include "black_scholes.h"

int main() {
    double S = 100.0;   // Spot price
    double K = 100.0;   // Strike price
    double T = 1.0;     // Time to maturity in years
    double r = 0.05;    // Risk-free interest rate
    double sigma = 0.2; // Volatility

    double callPrice = blackScholesCall(S, K, T, r, sigma);
    std::cout << "Black-Scholes Call Price: " << callPrice << std::endl;

    return 0;
}
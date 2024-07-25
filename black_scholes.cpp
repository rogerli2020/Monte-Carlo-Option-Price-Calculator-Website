#include "black_scholes.h"
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

double normalCDF(double value) 
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}

double blackScholesCall(double S, double K, double T, double r, double sigma) 
{
    if (S <= 0.0 || K <= 0.0 || T <= 0.0 || sigma <= 0.0) 
    {
        throw std::invalid_argument("Invalid input to Black-Scholes formula");
    }

    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    return S * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
}

double blackScholesPut(double S, double K, double T, double r, double sigma)
{
    double callPrice = blackScholesCall(S, K, T, r, sigma);
    double putPrice = callPrice + K * exp( -1*r*T ) - S;    // put call parity

    return putPrice;
}
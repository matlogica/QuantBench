# Phoenix Autocallable Note - Reference Implementation

## Overview
This repository provides a reference implementation for pricing and risk calculation of Phoenix Autocallable Notes using Monte Carlo simulation. It focuses on correctly handling the discontinuities present in these structured products' payoffs through smoothing techniques and Algorithmic Adjoint Differentiation (AADC).

## What is a Phoenix Autocallable Note?
A Phoenix Autocallable Note is a structured product linked to multiple underlying assets with three key features:
1. **Contingent Coupon Payments**: Periodic coupons paid only if all underlying assets are above their respective barrier levels
2. **Early Redemption (Autocall)**: Automatic early redemption if all underlying assets are at or above their initial levels on observation dates
3. **Principal Protection with Barrier**: Full principal returned at maturity if worst-performing asset remains above the barrier; otherwise, investor participates in the downside

## Key Features of the Implementation
- **Smoothed Digital Functions**: Implements continuous approximations of discontinuous payoff components
- **Monte Carlo Simulation**: Path-wise simulation with log-normal price models and correlation
- **Sensitivity Calculation**: Implements both bump-and-revalue and adjoint algorithmic differentiation approaches
- **Performance Analysis**: Comparison of convergence rates and computational efficiency

## Model Assumptions
This is a simplistic demonstration model with the following assumptions:
- Log-normal price dynamics for underlying assets
- Constant volatility and correlation
- No dividends or borrowing costs
- Constant risk-free rate
- No credit risk of issuer
- European-style payoff features (no early exercise features beyond autocalls)

## Core Components
- `PhoenixAutocallableNote`: Struct for storing the product parameters
- `pathPayoff`: Calculates the total payoff for a given price path
- `contLess`: Smoothed digital function that handles the discontinuities
- `monteCarloPhoenixNote`: Monte Carlo simulation for pricing and sensitivity calculation
- `calculateSensitivities`: Implements bump-and-revalue approach for Greeks calculation

## Usage
The repository includes several examples:
1. **Basic Pricing**: Demonstrates how to price a simple Phoenix Autocallable Note
2. **Sensitivity Analysis**: Shows how to calculate Greeks (volatility and correlation sensitivities)
3. **Convergence Analysis**: Compares accuracy across different path counts and smoothing parameters

## Focus on Discontinuities
The primary focus of this implementation is properly handling the discontinuities in the payoff function:
- The autocall feature creates a discontinuity when underlying assets cross their threshold levels
- The contingent coupon payments create discontinuities at the barrier levels
- The principal barrier creates a discontinuity at maturity

These discontinuities pose challenges for numerical methods, particularly when calculating sensitivities. The implementation demonstrates how smoothing techniques can improve convergence rates and how AADC provides more accurate and efficient sensitivity calculations.

## Performance Considerations
- Smoothing parameter selection affects both accuracy and convergence rate
- Path count requirements vary based on which risk metrics are most important
- AADC provides significant performance advantages for risk calculations

## Limitations
This is a reference implementation intended for educational purposes. Production systems would require additional features:
- More sophisticated market models (stochastic volatility, local volatility, etc.)
- Additional risk metrics (Greeks, stress tests, etc.)
- Market calibration capabilities
- More rigorous error analysis and control

## License
[MIT License](LICENSE)

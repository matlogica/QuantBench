#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include "aadc/aadc_tools.h"
#include "aadc/idouble.h"

#define AADC_ALLOW_TO_PASSIVE_BOOL

#include "aadc/ibool.h"
#include "aadc/aadc.h"

#include <random>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>

typedef idouble Real;


// digital approximation with two europeans, h should be calibrated in real system
// represents continous less function. If h=0, it is a step function

inline Real fastSigmoid(const Real& a, const Real& b, const Real& h) {
    if (h == 0) {
        // Handle the sharp transition case
        return iIf(a < b, 1.0, 0.0);
    }
    
    // Scale input to sigmoid by h
    Real x = (b - a) / (0.02*h);
    
    // Fast approximation using rational function
    Real sig = 0.5 * (x / (std::sqrt(1.0 + x*x)) + 1.0);
    return sig;
}

inline Real contLess(const Real& a, const Real& b, const Real& h) {
    if (h == 0) {
        return iIf(a < b , 1.0 , 0.0);
    }

    return fastSigmoid(a, b, h);
}

inline Real contNot(const Real& a) {
    return 1.0 - a;
}

// In the case of this specific Phoenix Autocallable Note, the contingent coupon works as follows:

// The coupon amount is 2.55% of the nominal amount ($1,000) per quarter, which equals $25.50 per note per quarter (or 10.20% annually)
// The coupon is only paid on a scheduled payment date if the reference level (price) of the worst-performing underlying asset (either S&P 500 or Apple stock, whichever has performed worse) is at or above its Coupon Barrier Level (which is set at 70% of its initial reference level)
// If the condition is not met (i.e., if either underlying drops below 70% of its initial level on an observation date), no coupon is paid for that period

// This is implemented for 2 assets, but we can generalize it to n assets
inline Real coupon(const Real& a1_init, const Real& a2_init, const Real& a1, const Real& a2, const Real& b1, const Real& b2, const Real& amount, const Real& H) {

    Real h1 = H* a1_init;
    Real h2 = H* a2_init;
    Real a1_perf = a1 / a1_init;
    Real a2_perf = a2 / a2_init;

    Real a1_is_worse = contLess(a1_perf, a2_perf, H * 1.0);

    Real pay_if_a1_is_worse = a1_is_worse * contLess(b1, a1, h1);
    Real pay_if_a2_is_worse = contNot(a1_is_worse) * contLess(b2, a2, h2);

    Real coupon = amount * (pay_if_a1_is_worse + pay_if_a2_is_worse);

    return coupon;
}

// Principal protection in the Phoenix Autocallable Note is a feature designed to protect your initial investment under certain conditions. Here's how it works in this specific product:

// Principal Barrier Level: This note has a principal barrier level set at 70% of the initial reference level for each underlying asset (the S&P 500 Index and Apple stock).
// Final Observation Only: Unlike the coupon and autocall features that are checked at multiple dates, principal protection is only evaluated at maturity (the Final Reference Valuation Date of April 27, 2018).
// How the protection works:

// If the note has not been autocalled earlier and reaches maturity
// AND if the Final Reference Level of the Worst-Performing Basket Constituent (either S&P 500 or Apple, whichever performed worse) is at or above its Principal Barrier Level (70% of initial level)
// THEN investors receive 100% of their initial investment ($1,000 per note)


// When protection fails:

// If the Final Reference Level of the Worst-Performing Basket Constituent is below its Principal Barrier Level (70% of initial)
// THEN investors receive a reduced amount calculated as:

// $1,000 × [Final Reference Level / Initial Reference Level]
// This means the investor participates directly in the downside performance of the worst-performing underlying

inline Real finalPrincipalProtection(const Real& a1_init, const Real& a2_init, const Real& a1, const Real& a2, const Real& b1, const Real& b2, const Real& amount, const Real& H) {

    Real h1 = H* a1_init;
    Real h2 = H* a2_init;

    Real a1_perf = a1 / a1_init;
    Real a2_perf = a2 / a2_init;

    Real pay_if_a1_is_worse = contLess(a1_perf, a2_perf, H * 1.0);
    Real pay_if_a2_is_worse = contNot( pay_if_a1_is_worse);

    Real a1_under = contLess(a1, b1, h1);
    Real a2_under = contLess(a2, b2, h2);

    Real principal_at_end = pay_if_a1_is_worse * (
        contNot(a1_under) * amount + a1_under * amount * a1 / a1_init
    )
    + pay_if_a2_is_worse * (
        contNot(a2_under) * amount + a2_under * amount * a2 / a2_init
    );

    return principal_at_end;

}

struct PhoenixAutocallableNote {
    Real a1_init;
    Real a2_init;
    std::vector<Real> coupon_amounts;
    Real a1_barrier, a2_barrier;
    Real a1_autocall_barrier, a2_autocall_barrier;

    Real nominal;

    Real principal_barrier;
};


inline Real pathPayoff(
    const PhoenixAutocallableNote& note,
    const std::vector<Real>& a1,
    const std::vector<Real>& a2,
    const Real& H=0.0
) {
    assert(a1.size() == a2.size());
    assert(note.coupon_amounts.size() == a1.size());

    Real total_amount = 0.0;
    Real not_autocalled_so_far = 1.0;

    Real h1 = H* note.a1_init;
    Real h2 = H* note.a2_init;

    for (int cp_i = 0; cp_i < note.coupon_amounts.size(); ++cp_i) {
        if (cp_i == 0) {
            // No auto call on first coupon date
            total_amount += not_autocalled_so_far * coupon(
                note.a1_init, note.a2_init,
                a1[cp_i], a2[cp_i],
                note.a1_barrier, note.a2_barrier,
                note.coupon_amounts[cp_i], H
            );
        }
        else if (cp_i < note.coupon_amounts.size()) {
            Real autocalled_on_this_coupon = contLess(note.a1_autocall_barrier, a1[cp_i], h1) * contLess(note.a2_autocall_barrier, a2[cp_i], h2);

            Real autocall_amount = autocalled_on_this_coupon * note.nominal;

            total_amount += not_autocalled_so_far * (
                autocall_amount + 
                coupon(
                    note.a1_init, note.a2_init,
                    a1[cp_i], a2[cp_i],
                    note.a1_barrier, note.a2_barrier,
                    note.coupon_amounts[cp_i], H
                )
            );
            not_autocalled_so_far *= contNot(autocalled_on_this_coupon);
        }
    }

    // Final payoff
    Real final_amount = finalPrincipalProtection(
        note.a1_init, note.a2_init,
        a1.back(), a2.back(),
        note.a1_barrier, note.a2_barrier,
        note.nominal, H
    );
    
    total_amount += not_autocalled_so_far * final_amount;
    return total_amount;
}


// Unit Tests

class PhoenixAutocallableNoteTest : public ::testing::Test {
    protected:
        void SetUp() override {
            // Create a basic Phoenix note based on the term sheet
            note.a1_init = 2095.15;  // S&P 500 initial level
            note.a2_init = 97.82;    // Apple initial level
            
            // Set barrier levels (70% of initial)
            note.a1_barrier = 0.7 * note.a1_init;  // S&P 500 barrier
            note.a2_barrier = 0.7 * note.a2_init;  // Apple barrier
            
            // Set autocall barriers (100% of initial)
            note.a1_autocall_barrier = 1.0 * note.a1_init;
            note.a2_autocall_barrier = 1.0 * note.a2_init;
            
            note.nominal = 1000.0;  // $1000 nominal amount
            note.principal_barrier = 0.7;  // 70% principal barrier
            
            // Set up 8 quarterly coupon amounts (2.55% each = $25.50)
            note.coupon_amounts = std::vector<Real>(8, 25.50);
        }
    
        PhoenixAutocallableNote note;
    };
    
    // Test helper functions first
    TEST(ContLessTest, BasicFunctionality) {
        // Test the sharp digital (h=0)
        EXPECT_EQ(contLess(1.0, 2.0, 0.0), 1.0);
        EXPECT_EQ(contLess(2.0, 1.0, 0.0), 0.0);
        EXPECT_EQ(contLess(1.0, 1.0, 0.0), 0.0);  // Equal case
        
        // Test the smooth digital (h=1.0)
        EXPECT_NEAR(AAD_PASSIVE(contLess(0.0, 30.0, 1.0)), 1.0, 1e-6);     // a << b
        EXPECT_NEAR(AAD_PASSIVE(contLess(30.0, 0.0, 1.0)), 0.0, 1e-6);     // a >> b

        EXPECT_NEAR(AAD_PASSIVE(contLess(0.0, 1.0, 1.0)), 1.0, 1e-3);     // a << b
        EXPECT_NEAR(AAD_PASSIVE(contLess(2.0, 1.0, 1.0)), 0.0, 1e-3);     // a >> b
        EXPECT_NEAR(AAD_PASSIVE(contLess(1.0, 1.0, 1.0)), 0.5, 1e-10);     // a = b
        EXPECT_NEAR(AAD_PASSIVE(contLess(0.5, 1.0, 1.0)), 1.0, 1e-2);    // a < b
        EXPECT_NEAR(AAD_PASSIVE(contLess(1.5, 1.0, 1.0)), 0.0, 1e-2);    // a > b
    }
    
    TEST(CouponTest, BasicFunctionality) {
        // Test coupon payment when Asset 1 is worse and above barrier
        Real a1_init = 100.0, a2_init = 100.0;
        Real a1 = 80.0, a2 = 90.0;  // Asset 1 is worse
        Real barrier1 = 70.0, barrier2 = 70.0;
        Real amount = 25.50;
        
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.0)), 25.50, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.1)), 25.50, 1e-2);
        
        // Test coupon payment when Asset 2 is worse and above barrier
        a1 = 90.0, a2 = 80.0;  // Asset 2 is worse
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.0)), 25.50, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.1)), 25.50, 1e-2);
        
        // Test no coupon when Asset 1 is worse and below barrier
        a1 = 60.0, a2 = 90.0;  // Asset 1 is worse and below barrier
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.0)), 0.0, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.1)), 0.0, 1e-2);
        
        // Test no coupon when Asset 2 is worse and below barrier
        a1 = 90.0, a2 = 60.0;  // Asset 2 is worse and below barrier
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.0)), 0.0, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(coupon(a1_init, a2_init, a1, a2, barrier1, barrier2, amount, 0.1)), 0.0, 1e-2);
    }
    
    TEST(FinalPrincipalProtectionTest, BasicFunctionality) {
        Real a1_init = 100.0, a2_init = 100.0;
        Real barrier1 = 70.0, barrier2 = 70.0;
        Real amount = 1000.0;
        
        // Test full principal when both assets above barrier
        Real a1 = 80.0, a2 = 90.0;  // Both above barrier, Asset 1 is worse
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.0)), 1000.0, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.1)), 1000.0, 0.1);
        
        // Test proportional loss when Asset 1 is worse and below barrier
        a1 = 60.0, a2 = 90.0;  // Asset 1 below barrier
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.0)), 600.0, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.1)), 600.0, 0.1);
        
        // Test proportional loss when Asset 2 is worse and below barrier
        a1 = 90.0, a2 = 50.0;  // Asset 2 below barrier
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.0)), 500.0, 1e-10);
        EXPECT_NEAR(AAD_PASSIVE(finalPrincipalProtection(a1_init, a2_init, a1, a2, barrier1, barrier2, amount,0.1)), 500.0, 0.1);
    }
    
    // Now test the full path payoff
    TEST_F(PhoenixAutocallableNoteTest, NoAutocallFullProtection) {
        // Create a path where no autocall occurs but both assets stay above barrier
        // SPX performs slightly worse but stays above 70% barrier
        std::vector<Real> spx_path = {
            0.95 * note.a1_init,  // 95% of initial
            0.90 * note.a1_init,  // 90% of initial
            0.85 * note.a1_init,  // 85% of initial
            0.80 * note.a1_init,  // 80% of initial
            0.78 * note.a1_init,  // 78% of initial
            0.76 * note.a1_init,  // 76% of initial
            0.74 * note.a1_init,  // 74% of initial
            0.72 * note.a1_init   // 72% of initial - still above 70% barrier
        };
        
        // AAPL performs better
        std::vector<Real> aapl_path = {
            0.98 * note.a2_init,  // 98% of initial
            0.95 * note.a2_init,  // 95% of initial
            0.92 * note.a2_init,  // 92% of initial
            0.89 * note.a2_init,  // 89% of initial
            0.87 * note.a2_init,  // 87% of initial
            0.85 * note.a2_init,  // 85% of initial
            0.83 * note.a2_init,  // 83% of initial
            0.80 * note.a2_init   // 80% of initial
        };
        
        Real payoff = pathPayoff(note, spx_path, aapl_path, 0.0);
        
        // Expected payoff = all 8 coupons + full principal = 8 * 25.50 + 1000 = 1204
        EXPECT_NEAR(AAD_PASSIVE(payoff), 1204.0, 1e-10);

        Real smooth_payoff = pathPayoff(note, spx_path, aapl_path, 0.1);

        // Expected payoff = all 8 coupons + full principal = 8 * 25.50 + 1000 = 1204
        EXPECT_NEAR(AAD_PASSIVE(smooth_payoff), 1204.0, 1.0);
    }
    
    TEST_F(PhoenixAutocallableNoteTest, EarlyAutocall) {
        // Create a path where an autocall occurs on the 3rd observation date
        std::vector<Real> spx_path = {
            0.95 * note.a1_init,   // 1st period: 95% of initial (no autocall check)
            0.98 * note.a1_init,   // 2nd period: 98% of initial (no autocall)
            1.05 * note.a1_init,   // 3rd period: 105% of initial (autocall occurs)
            1.10 * note.a1_init,   // These values shouldn't matter
            1.15 * note.a1_init,
            1.20 * note.a1_init,
            1.25 * note.a1_init,
            1.30 * note.a1_init
        };
        
        std::vector<Real> aapl_path = {
            0.96 * note.a2_init,   // 1st period: 96% of initial
            0.99 * note.a2_init,   // 2nd period: 99% of initial (no autocall)
            1.02 * note.a2_init,   // 3rd period: 102% of initial (autocall occurs)
            1.05 * note.a2_init,   // These values shouldn't matter
            1.08 * note.a2_init,
            1.10 * note.a2_init,
            1.12 * note.a2_init,
            1.15 * note.a2_init
        };
        
        Real payoff = pathPayoff(note, spx_path, aapl_path);
        
        // Expected payoff = 3 coupons + full principal = 3 * 25.50 + 1000 = 1076.5
        EXPECT_NEAR(AAD_PASSIVE(payoff), 1076.5, 1e-10);

        Real smooth_payoff = pathPayoff(note, spx_path, aapl_path, 0.1);
        // Expected payoff = 3 coupons + full principal = 3 * 25.50 + 1000 = 1076.5
        EXPECT_NEAR(AAD_PASSIVE(smooth_payoff), 1076.5, 0.1);
    }
    
    TEST_F(PhoenixAutocallableNoteTest, NoCouponNoAutocall) {
        // Create a path where no coupons are paid because the worst performer
        // is always below the barrier, and no autocall occurs
        std::vector<Real> spx_path = {
            0.65 * note.a1_init,  // 65% of initial (below barrier)
            0.63 * note.a1_init,  // 63% of initial (below barrier)
            0.60 * note.a1_init,  // 60% of initial (below barrier)
            0.58 * note.a1_init,  // 58% of initial (below barrier)
            0.56 * note.a1_init,  // 56% of initial (below barrier)
            0.54 * note.a1_init,  // 54% of initial (below barrier)
            0.52 * note.a1_init,  // 52% of initial (below barrier)
            0.50 * note.a1_init   // 50% of initial (below barrier)
        };
        
        // AAPL performs better but still not enough for autocall
        std::vector<Real> aapl_path = {
            0.90 * note.a2_init,
            0.88 * note.a2_init,
            0.86 * note.a2_init,
            0.84 * note.a2_init,
            0.82 * note.a2_init,
            0.80 * note.a2_init,
            0.78 * note.a2_init,
            0.75 * note.a2_init
        };
        
        Real payoff = pathPayoff(note, spx_path, aapl_path);
        
        // Expected payoff = no coupons + principal loss = 0 + 1000 * 0.5 = 500
        EXPECT_NEAR(AAD_PASSIVE(payoff), 500.0, 1e-10);

        Real smooth_payoff = pathPayoff(note, spx_path, aapl_path, 0.1);
        // Expected payoff = no coupons + principal loss = 0 + 1000 * 0.5 = 500
        EXPECT_NEAR(AAD_PASSIVE(smooth_payoff), 500.0, 0.1);
    }
    
    TEST_F(PhoenixAutocallableNoteTest, PartialCouponsWithPrincipalLoss) {
        // Create a path where some coupons are paid, no autocall, and principal barrier breached
        std::vector<Real> spx_path = {
            0.80 * note.a1_init,  // Above barrier
            0.75 * note.a1_init,  // Above barrier
            0.72 * note.a1_init,  // Above barrier
            0.68 * note.a1_init,  // Below barrier
            0.65 * note.a1_init,  // Below barrier
            0.62 * note.a1_init,  // Below barrier
            0.60 * note.a1_init,  // Below barrier
            0.58 * note.a1_init   // Below barrier at maturity
        };
        
        // AAPL performs better
        std::vector<Real> aapl_path = {
            0.85 * note.a2_init,
            0.82 * note.a2_init,
            0.80 * note.a2_init,
            0.78 * note.a2_init,
            0.76 * note.a2_init,
            0.74 * note.a2_init,
            0.72 * note.a2_init,
            0.70 * note.a2_init
        };
        
        Real payoff = pathPayoff(note, spx_path, aapl_path);
        
        // Expected payoff = 3 coupons + principal loss = 3 * 25.50 + 1000 * 0.58 = 656.5
        EXPECT_NEAR(AAD_PASSIVE(payoff), 656.5, 1e-10);

        Real smooth_payoff = pathPayoff(note, spx_path, aapl_path, 0.1);

        // Expected payoff = 3 coupons + principal loss = 3 * 25.50 + 1000 * 0.58 = 656.5
        EXPECT_NEAR(AAD_PASSIVE(smooth_payoff), 656.5, 0.1);
    }
    
    TEST_F(PhoenixAutocallableNoteTest, PartialCouponsWithAutocall) {
        // Create a path where some coupons are missed, then autocall occurs
        std::vector<Real> spx_path = {
            0.65 * note.a1_init,  // Below barrier (no coupon)
            0.75 * note.a1_init,  // Above barrier (coupon paid)
            0.85 * note.a1_init,  // Above barrier (coupon paid)
            0.95 * note.a1_init,  // Above barrier (coupon paid)
            1.05 * note.a1_init,  // Above barrier and above autocall (autocall occurs)
            1.10 * note.a1_init,  // Doesn't matter
            1.15 * note.a1_init,  // Doesn't matter
            1.20 * note.a1_init   // Doesn't matter
        };
        
        // AAPL performs similarly
        std::vector<Real> aapl_path = {
            0.70 * note.a2_init,  // At barrier
            0.80 * note.a2_init,  // Above barrier
            0.90 * note.a2_init,  // Above barrier
            0.95 * note.a2_init,  // Above barrier
            1.05 * note.a2_init,  // Above barrier and above autocall
            1.10 * note.a2_init,  // Doesn't matter
            1.15 * note.a2_init,  // Doesn't matter
            1.20 * note.a2_init   // Doesn't matter
        };
        
        Real payoff = pathPayoff(note, spx_path, aapl_path);
        
        // Expected payoff = 4 coupons (including the autocall date) + principal = 4 * 25.50 + 1000 = 1102
        EXPECT_NEAR(AAD_PASSIVE(payoff), 1102.0, 1e-10);

        Real smooth_payoff = pathPayoff(note, spx_path, aapl_path, 0.1);

        // Expected payoff = 4 coupons (including the autocall date) + principal = 4 * 25.50 + 1000 = 1102
        EXPECT_NEAR(AAD_PASSIVE(smooth_payoff), 1102.0, 0.1);
    }
    

// Function to generate independent normal random variables for Monte Carlo simulation
// These can be reused for AAD calculations to ensure consistent paths
inline std::vector<std::vector<std::pair<Real, Real>>> generateNormalRandoms(
    int numPaths,
    int numSteps,
    unsigned int seed = 12345
) {
    std::vector<std::vector<std::pair<Real, Real>>> randomVars(numPaths);
    
    std::mt19937 gen(seed);
    std::normal_distribution<double> normalDist(0.0, 1.0);
    
    for (int path = 0; path < numPaths; ++path) {
        randomVars[path].resize(numSteps);
        for (int step = 0; step < numSteps; ++step) {
            randomVars[path][step].first = normalDist(gen);
            randomVars[path][step].second = normalDist(gen);
        }
    }
    
    return randomVars;
}

// Monte Carlo simulation using pre-generated random variables
// This is useful for AAD calculations and variance reduction techniques
inline Real monteCarloPhoenixNoteWithRandoms(
    const PhoenixAutocallableNote& note,
    const std::vector<std::vector<std::pair<Real, Real>>>& randomVars,
    Real a1_spot,
    Real a2_spot,
    Real vol1,
    Real vol2,
    Real correlation,
    Real r,
    Real tMax,
    Real smoothingParam = 0.0,
    bool print_paths = false
) {
    int numPaths = randomVars.size();
    int numSteps = randomVars[0].size();
    
    // Time step
    Real dt = tMax / numSteps;
    Real sqrtDt = std::sqrt(dt);
    
    // Drift terms (risk-neutral measure)
    Real drift1 = (r - 0.5 * vol1 * vol1) * dt;
    Real drift2 = (r - 0.5 * vol2 * vol2) * dt;
    
    // Volatility terms
    Real volSqrtDt1 = vol1 * sqrtDt;
    Real volSqrtDt2 = vol2 * sqrtDt;
    
    // Results accumulator
    Real sumPayoffs = 0.0;
    
    for (int path = 0; path < numPaths; ++path) {
        // Initial asset prices for this path
        Real price1 = a1_spot;
        Real price2 = a2_spot;
        
        // Store asset price paths
        std::vector<Real> path1(numSteps + 1);
        std::vector<Real> path2(numSteps + 1);
        
        path1[0] = price1;
        path2[0] = price2;
        
        // Generate correlated asset paths using pre-generated random numbers
        for (int step = 0; step < numSteps; ++step) {
            // Get pre-generated standard normal random variables
            Real z1 = randomVars[path][step].first;
            Real z2 = randomVars[path][step].second;
            
            // Apply Cholesky decomposition to get correlated random variables
            Real w1 = z1;
            Real w2 = correlation * z1 + std::sqrt(1.0 - correlation * correlation) * z2;
            
            // Update asset prices using log-normal model
            price1 *= std::exp(drift1 + volSqrtDt1 * w1);
            price2 *= std::exp(drift2 + volSqrtDt2 * w2);
            
            path1[step + 1] = price1;
            path2[step + 1] = price2;
        }
        
        // Sample prices at observation dates
        std::vector<Real> observedPrices1(note.coupon_amounts.size());
        std::vector<Real> observedPrices2(note.coupon_amounts.size());
        
        // Map time steps to observation dates
        int stepsPerObs = numSteps / note.coupon_amounts.size();
        for (int i = 0; i < note.coupon_amounts.size(); ++i) {
            int step = (i + 1) * stepsPerObs;
            step = std::min(step, numSteps); // Ensure we don't go past the end
            observedPrices1[i] = path1[step];
            observedPrices2[i] = path2[step];
        }
        
        // Calculate payoff for this path
        Real pathPayoffResult = pathPayoff(note, observedPrices1, observedPrices2, smoothingParam);
        if(0) if (print_paths && path < 4) {
            std::cout << "Path " << path << ": Payoff = " << pathPayoffResult << std::endl;
        }
        sumPayoffs += pathPayoffResult;
    }
    
    // Return average payoff across all paths
    return sumPayoffs / numPaths;
}

// Structure to hold sensitivity results
struct SensitivityResults {

    template<class stream>
    void outputHeaders(stream& str, std::string prefix) {
        str 
            << prefix << "price_time,"
            << prefix << "risk_time,"
            << prefix << "price,"
            << prefix << "dPricedVol1,"
            << prefix << "dPricedVol2,"
            << prefix << "dPricedCorr,"
            << prefix << "dPricedSpot1,"
            << prefix << "dPricedSpot2,"
            << prefix << "gamma_spot1,"
            << prefix << "gamma_spot2"
        ;
    }
    template<class stream>
    void output(stream& str) {
        str
            << price_time << ","
            << risk_time << ","
            << price << ","
            << dPricedVol1 << ","
            << dPricedVol2 << ","
            << dPricedCorr << ","
            << dPricedSpot1 << ","
            << dPricedSpot2 << ","
            << gamma_spot1 << ","
            << gamma_spot2
        ;
    }


    double price_time;
    double risk_time;
    Real price;
    Real gamma_spot1;
    Real gamma_spot2;
    Real dPricedSpot1;
    Real dPricedSpot2;
    
    Real dPricedVol1;
    Real dPricedVol2;
    Real dPricedCorr;
};

// Calculate sensitivities using the same random variables for all evaluations
// This reduces variance in the sensitivity estimates
inline SensitivityResults calculateSensitivitiesWithFixedRandoms(
    const PhoenixAutocallableNote& note,
    const std::vector<std::vector<std::pair<Real, Real>>>& randomVars,
    Real a1_spot,
    Real a2_spot,
    Real vol1,
    Real vol2,
    Real correlation,
    Real r,
    Real tMax,
    Real smoothingParam = 0.0,
    Real vol_bump_size = 0.001,
    Real spot_bump_size = 1.0 // %
) {
    SensitivityResults results;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
       
    // Base price
    results.price = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1, vol2, correlation, r, tMax, smoothingParam, false
    );

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

    results.price_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6;

//    return results;

    Real a1_spot_bump_size = a1_spot * spot_bump_size / 100.0;

    // Calculate dPricedSpot1 using bump-and-revalue with the same random variables
    Real bumpedPriceSpotUp1 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot + a1_spot_bump_size, a2_spot, vol1, vol2, correlation, r, tMax, smoothingParam
    );
    Real bumpedPriceSpotDown1 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot - a1_spot_bump_size, a2_spot, vol1, vol2, correlation, r, tMax, smoothingParam
    );
    results.dPricedSpot1 = (bumpedPriceSpotUp1 - bumpedPriceSpotDown1) / (2.0 * a1_spot_bump_size);
    results.gamma_spot1 = (bumpedPriceSpotUp1 - 2.0 * results.price + bumpedPriceSpotDown1) / (a1_spot_bump_size * a1_spot_bump_size);

    // Calculate dPricedSpot2 using bump-and-revalue with the same random variables

    Real a2_spot_bump_size = a2_spot * spot_bump_size / 100.0;

    Real bumpedPriceSpotUp2 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot + a2_spot_bump_size, vol1, vol2, correlation, r, tMax, smoothingParam
    );
    Real bumpedPriceSpotDown2 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot - a2_spot_bump_size, vol1, vol2, correlation, r, tMax, smoothingParam
    );
    results.dPricedSpot2 = (bumpedPriceSpotUp2 - bumpedPriceSpotDown2) / (2.0 * a2_spot_bump_size);
    results.gamma_spot2 = (bumpedPriceSpotUp2 - 2.0 * results.price + bumpedPriceSpotDown2) / (a2_spot_bump_size * a2_spot_bump_size);

    // Calculate dPricedVol1 using bump-and-revalue with the same random variables
    Real bumpedPriceVolUp1 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1 + vol_bump_size, vol2, correlation, r, tMax, smoothingParam
    );
    Real bumpedPriceVolDown1 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1 - vol_bump_size, vol2, correlation, r, tMax, smoothingParam
    );

    results.dPricedVol1 = (bumpedPriceVolUp1 - bumpedPriceVolDown1) / (2.0 * vol_bump_size);
    
    // Calculate dPricedVol2 using bump-and-revalue with the same random variables
    Real bumpedPriceVolUp2 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1, vol2 + vol_bump_size, correlation, r, tMax, smoothingParam
    );
    Real bumpedPriceVolDown2 = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1, vol2 - vol_bump_size, correlation, r, tMax, smoothingParam
    );
    results.dPricedVol2 = (bumpedPriceVolUp2 - bumpedPriceVolDown2) / (2.0 * vol_bump_size);
    
    // Calculate dPricedCorr using bump-and-revalue with the same random variables
    Real bumpedPriceCorrUp = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1, vol2, correlation + 1. * vol_bump_size, r, tMax, smoothingParam
    );
    Real bumpedPriceCorrDown = monteCarloPhoenixNoteWithRandoms(
        note, randomVars, a1_spot, a2_spot, vol1, vol2, correlation - 1. * vol_bump_size, r, tMax, smoothingParam
    );
    results.dPricedCorr = (bumpedPriceCorrUp - bumpedPriceCorrDown) / (2.0 * 1. * vol_bump_size);

    std::chrono::high_resolution_clock::time_point end2 = std::chrono::high_resolution_clock::now();

    results.risk_time = std::chrono::duration_cast<std::chrono::microseconds>(end2 - end).count() / 1e6;
    
    return results;
}


typedef __m256d mmType;

struct AADCPricingKernel {
    std::shared_ptr<aadc::AADCFunctions<mmType> > aadc_funcs;
    std::vector<std::pair<aadc::AADCArgument, aadc::AADCArgument> > random_normals;

    aadc::AADCArgument spot1_arg;
    aadc::AADCArgument spot2_arg;

    aadc::AADCArgument vol1_arg;
    aadc::AADCArgument vol2_arg;
    aadc::AADCArgument correlation_arg;
    aadc::AADCArgument r_arg;

    aadc::AADCResult path_price_result;
};

std::shared_ptr<AADCPricingKernel> recordAADCKernel(
    const PhoenixAutocallableNote& note,
    const std::vector<std::vector<std::pair<Real, Real>>>& randomVars,
    Real a1_spot,
    Real a2_spot,
    Real vol1,
    Real vol2,
    Real correlation,
    Real r,
    Real tMax,
    Real smoothingParam
) {
    auto res = std::make_shared<AADCPricingKernel>();
    // Create an AADC function object
    res->aadc_funcs = std::make_shared<aadc::AADCFunctions<mmType> >();
    
    // Set up the AADC function with the parameters

    std::vector<std::vector<std::pair<Real, Real>>> onePathRandomVars(1);
    onePathRandomVars[0] = randomVars[0];

    res->random_normals.resize(onePathRandomVars[0].size());

    res->aadc_funcs->startRecording();

    for (int i = 0; i < onePathRandomVars[0].size(); ++i) {
        res->random_normals[i].first = onePathRandomVars[0][i].first.markAsInputNoDiff();
        res->random_normals[i].second = onePathRandomVars[0][i].second.markAsInputNoDiff();
    }

    res->spot1_arg = a1_spot.markAsInput();
    res->spot2_arg = a2_spot.markAsInput();
    res->vol1_arg = vol1.markAsInput();
    res->vol2_arg = vol2.markAsInput();
    res->correlation_arg = correlation.markAsInput();
    res->r_arg = r.markAsInput();

    auto path_payoff = monteCarloPhoenixNoteWithRandoms(
        note, onePathRandomVars, a1_spot, a2_spot, vol1, vol2, correlation, r, tMax, smoothingParam
    );

    res->path_price_result = path_payoff.markAsOutput();

    std::cout << "payoff : " << path_payoff << std::endl;

    res->aadc_funcs->stopRecording();

    assert(res->aadc_funcs->getPassiveWarnings().size() == 0);

    if (res->aadc_funcs->getPassiveWarnings().size() > 0) {
        std::cout << "AADCPricingKernel: Passive warnings: " << std::endl;
        for (const auto& warning : res->aadc_funcs->getPassiveWarnings()) {
            std::cout << warning.location << ":" << warning.line << std::endl;
        }
    }
    
    return res;
}


auto calculateSensitivitiesWithFixedRandomsAADC(
    std::shared_ptr<AADCPricingKernel> kernel,
    const std::vector<std::vector<std::pair<Real, Real>>>& randomVars,
    Real a1_spot,
    Real a2_spot,
    Real vol1,
    Real vol2,
    Real correlation,
    Real r,
    double spot_bump_size = 5.0 // %
) {
    SensitivityResults results;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    auto ws = kernel->aadc_funcs->createWorkSpace();

    auto num_paths = randomVars.size();
    assert(num_paths % aadc::mmSize<mmType>() == 0);
    auto num_avx_paths = num_paths / aadc::mmSize<mmType>();

    ws->setVal(kernel->vol1_arg, vol1.val);
    ws->setVal(kernel->vol2_arg, vol2.val);
    ws->setVal(kernel->correlation_arg, correlation.val);
    ws->setVal(kernel->r_arg, r.val);

    results.price = 0.;
    results.dPricedSpot1 = 0.;
    results.dPricedSpot2 = 0.;
    results.dPricedVol1 = 0.;
    results.dPricedVol2 = 0.;
    results.dPricedCorr = 0.;

    Real dPricedSpot1_Up = 0.;
    Real dPricedSpot2_Up = 0.;
    Real dPricedSpot1_Down = 0.;
    Real dPricedSpot2_Down = 0.;

    Real a1_spot_bump_size = a1_spot * spot_bump_size / 100.0;
    Real a2_spot_bump_size = a2_spot * spot_bump_size / 100.0;

    for (int i = 0; i < num_avx_paths; ++i) {
        // Set normals for 4 paths
        for (int j = 0; j < kernel->random_normals.size(); ++j) {
            for (int k = 0; k < aadc::mmSize<mmType>(); ++k) {
                ws->valp(kernel->random_normals[j].first)[k] = randomVars[i * aadc::mmSize<mmType>() + k][j].first.val;
                ws->valp(kernel->random_normals[j].second)[k] = randomVars[i * aadc::mmSize<mmType>() + k][j].second.val;
            }
        }

        for (int ext_scen_i = 0; ext_scen_i < 5; ++ext_scen_i) { // 0 - normal, 1 - spot1 up, 2 - spot2 up, 3 - spot1 down, 4 - spot2 down
            ws->setVal(kernel->spot1_arg, a1_spot.val);
            ws->setVal(kernel->spot2_arg, a2_spot.val);
            if (ext_scen_i == 1) {
                ws->setVal(kernel->spot1_arg, a1_spot.val + a1_spot_bump_size.val);
            } else if (ext_scen_i == 2) {
                ws->setVal(kernel->spot2_arg, a2_spot.val + a2_spot_bump_size.val);
            } else if (ext_scen_i == 3) {
                ws->setVal(kernel->spot1_arg, a1_spot.val - a1_spot_bump_size.val);
            } else if (ext_scen_i == 4) {
                ws->setVal(kernel->spot2_arg, a2_spot.val - a2_spot_bump_size.val);
            }
            // Calculate the price
            kernel->aadc_funcs->forward(*ws);
            if (ext_scen_i == 0) {
                results.price += aadc::mmSum(ws->val(kernel->path_price_result));
            }

            if(0) if (i == 0) {
                for (int k = 0; k < aadc::mmSize<mmType>(); ++k) {
                    std::cout << "AADC Path " << i * aadc::mmSize<mmType>() + k << " price: " << ws->valp(kernel->path_price_result)[k] << std::endl;
                }
            }

            // Calculate the sensitivities
            ws->setDiff(kernel->path_price_result, 1.0);
            kernel->aadc_funcs->reverse(*ws);
            if (ext_scen_i == 0) {
                results.dPricedSpot1 += aadc::mmSum(ws->diff(kernel->spot1_arg));
                results.dPricedSpot2 += aadc::mmSum(ws->diff(kernel->spot2_arg));
                results.dPricedVol1 += aadc::mmSum(ws->diff(kernel->vol1_arg));
                results.dPricedVol2 += aadc::mmSum(ws->diff(kernel->vol2_arg));
                results.dPricedCorr += aadc::mmSum(ws->diff(kernel->correlation_arg));
            } else if (ext_scen_i == 1) {
                dPricedSpot1_Up += aadc::mmSum(ws->diff(kernel->spot1_arg));
            } else if (ext_scen_i == 2) {
                dPricedSpot2_Up += aadc::mmSum(ws->diff(kernel->spot2_arg));
            } else if (ext_scen_i == 3) {
                dPricedSpot1_Down += aadc::mmSum(ws->diff(kernel->spot1_arg));
            } else if (ext_scen_i == 4) {
                dPricedSpot2_Down += aadc::mmSum(ws->diff(kernel->spot2_arg));
            }
        }
    }

    results.price /= num_paths;
    results.dPricedSpot1 /= num_paths;
    results.dPricedSpot2 /= num_paths;
    results.dPricedVol1 /= num_paths;
    results.dPricedVol2 /= num_paths;
    results.dPricedCorr /= num_paths;

    dPricedSpot1_Up /= num_paths;
    dPricedSpot1_Down /= num_paths;
    dPricedSpot2_Up /= num_paths;
    dPricedSpot2_Down /= num_paths;
    results.gamma_spot1 = (dPricedSpot1_Up - dPricedSpot1_Down) / (2.0 * a1_spot_bump_size);
    results.gamma_spot2 = (dPricedSpot2_Up - dPricedSpot2_Down) / (2.0 * a2_spot_bump_size);

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

    results.price_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6;
    results.risk_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6;
    return results;
}
    

class ProgressBar {
private:
    int total;
    int barWidth;

public:
    ProgressBar(int total, int width = 70) : total(total), barWidth(width) {}
    
    void update(int current) {
        float progress = (float)current / total;
        int pos = barWidth * progress;
        
        std::cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")\r";
        std::cout.flush();
        
        if (current == total) std::cout << std::endl;
    }
};

// Usage example:
// ProgressBar bar(numPaths);
// for (int i = 0; i < numPaths; ++i) {
//     // Run simulation
//     bar.update(i+1);
// }

struct MarketParameters {
    Real a1_spot;
    Real a2_spot;
    Real vol1;
    Real vol2;
    Real correlation;
    Real r;

    Real scenario_val;
};

void ladderScenario(const MarketParameters& base, MarketParameters& scenario, int step, int num_steps, const std::string ladder) {
    scenario.a1_spot = base.a1_spot;
    scenario.a2_spot = base.a2_spot;
    scenario.vol1 = base.vol1;
    scenario.vol2 = base.vol2;
    scenario.correlation = base.correlation;
    scenario.r = base.r;
    if (ladder == "vol1") {
        scenario.scenario_val = scenario.vol1 = base.vol1 * (0.5 + double(step) / num_steps);
    } else if (ladder == "vol2") {
        scenario.scenario_val = scenario.vol2 = base.vol2 * (0.5 + double(step) / num_steps);
    } else if (ladder == "corr") {
        scenario.scenario_val = scenario.correlation = base.correlation * (0.5 + double(step) / num_steps);
    } else if (ladder == "spot1") {
        scenario.scenario_val = scenario.a1_spot = base.a1_spot * (0.5 + double(step) / num_steps);
    } else if (ladder == "spot2") {
        scenario.scenario_val = scenario.a2_spot = base.a2_spot * (0.5 + double(step) / num_steps);
    } else if (ladder == "r") {
        scenario.scenario_val = scenario.r = base.r * (0.5 + double(step) / num_steps);
    } else {
        std::cerr << "Invalid ladder type: " << ladder << std::endl;
    }
}

// Example usage
void runMonteCarloExample(int params_set_start, int params_set_end) {
    // Create a PhoenixAutocallableNote similar to the one in the term sheet
    PhoenixAutocallableNote note;
    note.a1_init = 2095.15;  // S&P 500 initial
    note.a2_init = 97.82;    // Apple initial
    
    // Barriers at 70% of initial
    note.a1_barrier = 0.7 * note.a1_init;
    note.a2_barrier = 0.7 * note.a2_init;
    
    // Autocall barriers at 100% of initial
    note.a1_autocall_barrier = 1.0 * note.a1_init;
    note.a2_autocall_barrier = 1.0 * note.a2_init;
    
    note.nominal = 1000.0;
    note.principal_barrier = 0.7;
    
    // 8 quarterly coupons of 2.55% each
    note.coupon_amounts = std::vector<Real>(8, 25.50);
    
    // Market parameters
    MarketParameters market;
    market.a1_spot = note.a1_init;  // Current S&P 500 price
    market.a2_spot = note.a2_init;  // Current Apple price
    market.vol1 = 0.2;         // 20% annualized volatility for S&P 500
    market.vol2 = 0.3;         // 30% annualized volatility for Apple
    market.correlation = 0.6;  // 60% correlation between S&P 500 and Apple
    market.r = 0.01;           // 1% risk-free rate
    Real tMax = 2.0;         // 2 years to maturity
    
    // Monte Carlo parameters
//    int numPaths = 4;    // Number of paths
    int numPaths = 100000;    // Number of paths
    int numSteps = 96;       // 96 time steps (4 steps per week for 24 months)
    Real smoothingParam = 0.4; // Smoothing parameter for digital functions
    
    int seed = 12345; // Seed for random number generation

    // Generate random variables once and reuse
    auto randomVars = generateNormalRandoms(numPaths, numSteps, seed);

    // Calculate price and sensitivities
    auto results = calculateSensitivitiesWithFixedRandoms(
        note, randomVars, market.a1_spot, market.a2_spot, market.vol1, market.vol2, market.correlation, market.r, tMax, 0.0
    );
    
    std::cout << "Phoenix Autocallable Note - Monte Carlo Results:" << std::endl;
    std::cout << "Price: " << results.price << std::endl;
    std::cout << "dPricedVol1 (S&P 500): " << results.dPricedVol1 << std::endl;
    std::cout << "dPricedVol2 (Apple): " << results.dPricedVol2 << std::endl;
    std::cout << "dPricedSpot1 (S&P 500): " << results.dPricedSpot1 << std::endl;
    std::cout << "dPricedSpot2 (Apple): " << results.dPricedSpot2 << std::endl;
    std::cout << "dPricedSpot1 (S&P 500) gamma: " << results.gamma_spot1 << std::endl;
    std::cout << "dPricedSpot2 (Apple) gamma: " << results.gamma_spot2 << std::endl;
    std::cout << "dPricedCorr: " << results.dPricedCorr << std::endl;

    // Calculate price with different smoothing parameter
    for (int bump_i = 0; bump_i < 5; ++bump_i)        
    {
        double bump_size = 1.0 + bump_i;
        std::cout << "============================== FD_BUMP: " << bump_size << "%" << std::endl;
        auto results = calculateSensitivitiesWithFixedRandoms(
            note, randomVars, market.a1_spot, market.a2_spot, market.vol1, market.vol2, market.correlation, market.r, tMax, smoothingParam, 0.001, bump_size
        );
        
        std::cout << "Phoenix Autocallable Note - Monte Carlo Results:" << std::endl;
        std::cout << "Price: " << results.price << std::endl;
        std::cout << "dPricedVol1 (S&P 500): " << results.dPricedVol1 << std::endl;
        std::cout << "dPricedVol2 (Apple): " << results.dPricedVol2 << std::endl;
        std::cout << "dPricedSpot1 (S&P 500): " << results.dPricedSpot1 << std::endl;
        std::cout << "dPricedSpot2 (Apple): " << results.dPricedSpot2 << std::endl;
        std::cout << "dPricedSpot1 (S&P 500) gamma: " << results.gamma_spot1 << std::endl;
        std::cout << "dPricedSpot2 (Apple) gamma: " << results.gamma_spot2 << std::endl;
        std::cout << "dPricedCorr: " << results.dPricedCorr << std::endl;
    }

    auto aadc_kernel = recordAADCKernel(
        note, randomVars, market.a1_spot, market.a2_spot, market.vol1, market.vol2, market.correlation, market.r, tMax, smoothingParam
    );

    auto aadc_results = calculateSensitivitiesWithFixedRandomsAADC(
        aadc_kernel, randomVars, market.a1_spot, market.a2_spot, market.vol1, market.vol2, market.correlation, market.r
    );

    std::cout << "AADCPricing - Monte Carlo Results:" << std::endl;
    std::cout << "Price: " << aadc_results.price << std::endl;
    std::cout << "dPricedVol1 (S&P 500): " << aadc_results.dPricedVol1 << std::endl;
    std::cout << "dPricedVol2 (Apple): " << aadc_results.dPricedVol2 << std::endl;
    std::cout << "dPricedSpot1 (S&P 500): " << aadc_results.dPricedSpot1 << std::endl;
    std::cout << "dPricedSpot2 (Apple): " << aadc_results.dPricedSpot2 << std::endl;
    std::cout << "dPricedSpot1 (S&P 500) gamma: " << aadc_results.gamma_spot1 << std::endl;
    std::cout << "dPricedSpot2 (Apple) gamma: " << aadc_results.gamma_spot2 << std::endl;
    std::cout << "dPricedCorr: " << aadc_results.dPricedCorr << std::endl;
    std::cout << "==============================" << std::endl;


    // Run pricing for range of vol1 values
    struct ParamSet {
        Real smoothingParam;
        int numPaths;
        std::string ladder;
    };

    std::vector<std::string> ladders = {
        "spot1", "spot2", "vol1", "vol2", "corr"
    };

    std::vector<int> numPathsList = {
        10000, 50000, 100000, 500000, 1000000
    };

    std::vector<double> smoothingParams = {
        0.0, 0.1, 0.4, 0.8, 1.6, 2.4, 3.5
};

    std::vector<ParamSet> params;

    for (const auto& ladder : ladders) {
        for (const auto& numPaths : numPathsList) {
            for (const auto& smoothingParam : smoothingParams) {
                params.push_back({smoothingParam, numPaths, ladder});
            }
        }
    }
    // Add "truth ground with 10m paths"
    for (const auto& ladder : ladders) {
            params.push_back({0.8, 10000000, ladder});
    }
    std::cout << "Number of parameter sets: " << params.size() << std::endl;


    for (int param_set_i = params_set_start; param_set_i < std::min<int>(params.size(), params_set_end); ++param_set_i) {
        const auto& param = params[param_set_i];

        std::stringstream run_name;
        run_name << "run_" << param.smoothingParam << "_" << param.numPaths;
        std::cout << "==============================" << std::endl;
        std::cout << "Running with smoothingParam: " << param.smoothingParam << ", numPaths: " << param.numPaths << std::endl;

        std::string run_name_str = run_name.str();

        Real smoothingParam = param.smoothingParam;
        int numPaths = param.numPaths;

        // Rerecord aadc kernel because smoothingParam is hardcoded
        auto aadc_kernel = recordAADCKernel(
            note, randomVars, market.a1_spot, market.a2_spot, market.vol1, market.vol2, market.correlation, market.r, tMax, smoothingParam
        );
    
        // Generate random variables once and reuse
        auto randomVars = generateNormalRandoms(numPaths, numSteps, seed);

        std::stringstream params_out;

        params_out << "Sensitivity Analysis Parameters:" << std::endl;
        params_out << "smoothingParam: " << smoothingParam << std::endl;
        params_out << "numPaths: " << numPaths << std::endl;


        std::stringstream csv_out;

        csv_out << param.ladder;

        csv_out << ","; results.outputHeaders(csv_out, "");
        csv_out << ","; results.outputHeaders(csv_out, "smooth_");
        csv_out << ","; results.outputHeaders(csv_out, "aadc_");
        csv_out << std::endl;

        int num_steps = 50;

        ProgressBar bar(num_steps);

        for (int scen_i = 0; scen_i < num_steps; ++scen_i) {
            bar.update(scen_i+1);

            MarketParameters scenario;

            ladderScenario(market, scenario, scen_i, num_steps, param.ladder);

            auto results = calculateSensitivitiesWithFixedRandoms(
                note, randomVars, scenario.a1_spot, scenario.a2_spot, scenario.vol1, scenario.vol2, scenario.correlation, scenario.r, tMax, 0.0
            );
            auto results_smooth = calculateSensitivitiesWithFixedRandoms(
                note, randomVars, scenario.a1_spot, scenario.a2_spot, scenario.vol1, scenario.vol2, scenario.correlation, scenario.r, tMax, smoothingParam
            );
            auto aadc_results = calculateSensitivitiesWithFixedRandomsAADC(
                aadc_kernel, randomVars, scenario.a1_spot, scenario.a2_spot, scenario.vol1, scenario.vol2, scenario.correlation, scenario.r
            );
            csv_out << scenario.scenario_val << ",";
            results.output(csv_out);
            csv_out << ",";
            results_smooth.output(csv_out);
            csv_out << ",";
            aadc_results.output(csv_out);
            csv_out << std::endl;
        }
        std::cout << "==============================" << std::endl;

        std::ofstream file(std::string("sensitivity_results_") + param.ladder + "_" + run_name_str + ".csv");
        if (file.is_open()) {
            file << csv_out.str();
            file.close();
            std::cout << "Results saved to sensitivity_results_" << run_name_str << ".csv" << std::endl;
        } else {
            std::cerr << "Error opening file for writing." << std::endl;
        }
        std::ofstream params_file(std::string("sensitivity_results_") + run_name_str + ".txt");
        if (params_file.is_open()) {
            params_file << params_out.str();
            params_file.close();
            std::cout << "Parameters saved to sensitivity_params_" << run_name_str << ".txt" << std::endl;
        } else {
            std::cerr << "Error opening file for writing." << std::endl;
        }

        std::cout << "==============================" << std::endl;
    } // params
}


TEST(AADCPricing, asset_autocall) {


    if (std::getenv("AADCPARAMSET") != nullptr) {
        int params_set_start = std::stoi(std::getenv("AADCPARAMSET"));
        runMonteCarloExample(params_set_start, params_set_start+1);
    } else {
        runMonteCarloExample(0, 10000000);
    }        

}


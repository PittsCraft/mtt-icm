#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <vector>
#include "../include/mtt-icm/monte-carlo-icm.hpp"
#include "../src/monte-carlo-icm-util.hpp"

using boost::format;
using namespace std;

BOOST_AUTO_TEST_SUITE(icm_monte_carlo)

    BOOST_AUTO_TEST_CASE(values)
    {
        vector<double> stacks = {98575., 92424., 91638., 91481., 85173.,
                                 82184., 81072., 80449., 78163., 78102.};
        vector<double> payouts = {50, 20, 10, 4, 2,
                                  1, 0, 0, 0, 0};
        // Values from https://www.holdemresources.net/icmcalculator
        vector<double> expectedResult = {9.83818, 9.29185, 9.22142, 9.20733, 8.63667,
                                         8.36304, 8.2607, 8.20324, 7.99161, 7.98595};
        auto result = icm::monteCarloIcm(stacks, payouts, 500000);
        // Buyin = sum(payouts) / nb_players = 8,7
        // Error threshold = 0.1 = buy_in / 87
        for (int i = 0; i < 10; i++) {
            auto error = fabs(result[i] - expectedResult[i]);
            BOOST_TEST(error <= 1E-1, format("Monte-Carlo ICM value for player %1%, error : %2%") % i % error);
        }
    }

    BOOST_AUTO_TEST_CASE(permutation_total_sort_constant_rng)
    {
        // MC Permutation should sort in reverse order according to weights if RNG is a constant > 1 fn
        const int size = 15;
        int destination[size];
        double weights[size];
        icm::mc_util::RNG rng = []() -> double {
            return 2;
        };
        for (int i = 0; i < size; ++i) {
            weights[i] = i;
        }
        auto sorter = icm::mc_util::defaultSorter();
        icm::mc_util::monteCarloPermutation(weights, destination, size, rng, sorter);
        for (int i = 0; i < size; ++i) {
            BOOST_TEST(destination[i] == size - i - 1);
        }
    }

    BOOST_AUTO_TEST_CASE(permutation_total_sort_constant_rng_reverse)
    {
        // Let's ensure the permutation correctly behaves for descending weights
        const int size = 15;
        int destination[size];
        double weights[size];
        icm::mc_util::RNG rng = []() -> double {
            return 2;
        };
        for (int i = 0; i < size; ++i) {
            weights[i] = size - i;
        }
        auto sorter = icm::mc_util::defaultSorter();
        icm::mc_util::monteCarloPermutation(weights, destination, size, rng, sorter);
        for (int i = 0; i < size; ++i) {
            BOOST_TEST(destination[i] == i);
        }
    }

    BOOST_AUTO_TEST_CASE(permutation_partial_sort_constant_rng)
    {
        // Should work as well with partial sort for the first relevant ranks
        const int size = 15;
        int destination[size];
        double weights[size];
        icm::mc_util::RNG rng = []() -> double {
            return 2;
        };
        for (int i = 0; i < size; ++i) {
            weights[i] = i;
        }
        int relevantRanksCount = 10;
        auto partialSorter = icm::mc_util::partialSorter(relevantRanksCount);
        icm::mc_util::monteCarloPermutation(weights, destination, size, rng, partialSorter);
        for (int i = 0; i < relevantRanksCount; ++i) {
            BOOST_TEST(destination[i] == size - i - 1);
        }
    }

    BOOST_AUTO_TEST_CASE(first_zero_payout)
    {
        double payoutsNormal[] = {1, 2, 3, 4, 5, 0, 0, 0};
        BOOST_TEST(icm::mc_util::firstZeroPayout(payoutsNormal, 8) == 5, "First zero index for payouts with trailing zeros");

        double payoutsWithHole[] = {1, 2, 3, 4, 5, 0, 1, 0};
        BOOST_TEST(icm::mc_util::firstZeroPayout(payoutsWithHole, 8) == 7, "First zero index for payouts with a hole before trailing zeros");

        double payoutsNoZeros[] = {1, 2, 3, 4, 5, 6, 7, 8};
        BOOST_TEST(icm::mc_util::firstZeroPayout(payoutsNoZeros, 8) < 0, "First zero index for payouts with no null payout");

        double payoutsNoZerosExceptHole[] = {1, 2, 0, 0, 5, 6, 7, 8};
        BOOST_TEST(icm::mc_util::firstZeroPayout(payoutsNoZerosExceptHole, 8) < 0, "First zero index for payouts with no null payout except holes");
    }
BOOST_AUTO_TEST_SUITE_END()
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <vector>
#include "../include/mtt-icm/icm.hpp"

using boost::format;
using namespace std;

BOOST_AUTO_TEST_SUITE(icm)

    BOOST_AUTO_TEST_CASE(values)
    {
        vector<double> stacks = {98575., 92424., 91638., 91481., 85173.,
                                 82184., 81072., 80449., 78163., 78102.};
        vector<double> payouts = {50, 20, 10, 4, 2,
                                  1, 0, 0, 0, 0};
        // Values from https://www.holdemresources.net/icmcalculator
        vector<double> expectedResult = {9.83818, 9.29185, 9.22142, 9.20733, 8.63667,
                                         8.36304, 8.2607, 8.20324, 7.99161, 7.98595};
        auto result = icm::icm(stacks, payouts);
        for (int i = 0; i < 10; i++) {
            auto error = fabs(result[i] - expectedResult[i]);
            BOOST_TEST(error <= 1E-5, format("ICM value for player %1%, error : %2%") % i % error);
        }
    }
}
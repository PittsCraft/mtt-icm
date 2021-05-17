#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include "../src/math.hpp"

using boost::format;

BOOST_AUTO_TEST_SUITE(math)

    BOOST_AUTO_TEST_CASE(normal_cdf_inverse_test)
    {
        // See the demo() function from the source : https://www.johndcook.com/blog/cpp_phi_inverse/
        double p[] =
                {
                        0.0000001,
                        0.00001,
                        0.001,
                        0.05,
                        0.15,
                        0.25,
                        0.35,
                        0.45,
                        0.55,
                        0.65,
                        0.75,
                        0.85,
                        0.95,
                        0.999,
                        0.99999,
                        0.9999999
                };

        // Exact values computed by Mathematica.
        double exact[] =
                {
                        -5.199337582187471,
                        -4.264890793922602,
                        -3.090232306167813,
                        -1.6448536269514729,
                        -1.0364333894937896,
                        -0.6744897501960817,
                        -0.38532046640756773,
                        -0.12566134685507402,
                        0.12566134685507402,
                        0.38532046640756773,
                        0.6744897501960817,
                        1.0364333894937896,
                        1.6448536269514729,
                        3.090232306167813,
                        4.264890793922602,
                        5.199337582187471
                };

        int numValues = sizeof(p)/sizeof(double);
        for (int i = 0; i < numValues; ++i)
        {
            double computed = NormalCDFInverse(p[i]);
            double error = fabs(exact[i] - computed);
            BOOST_TEST(error <= 4.5E-4, format("NormalCDFInverse for input %1%, error : %2%") % p[i] % error);
        }
    }

BOOST_AUTO_TEST_SUITE_END()

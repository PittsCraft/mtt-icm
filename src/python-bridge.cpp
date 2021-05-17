#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/format.hpp>
#include <utility>
#include "../include/mtt-icm/icm.hpp"
#include "../include/mtt-icm/monte-carlo-icm.hpp"

using namespace std;
namespace np = boost::python::numpy;
namespace py = boost::python;

vector<double> numpyToDouble(const np::ndarray &array) {
    const int strides = array.strides(0);
    const int size = array.shape(0);
    vector<double> result(size);
    for (int i = 0; i < size; ++i) {
        result[i] = ((double *) (array.get_data() + i * strides))[0];
    }
    return result;
}

void reportDoubleToNumpy(vector<double> src, const np::ndarray &array) {
    const int strides = array.strides(0);
    const int size = array.shape(0);
    for (int i = 0; i < size; i++) {
        ((double *) (array.get_data() + i * strides))[0] = src[i];
    }
}

vector<vector<double>> numpyToDouble2d(const np::ndarray &array) {
    const int firstDimStrides = array.strides(0);
    const int secondDimStrides = array.strides(1);
    const int firstDimSize = array.shape(0);
    const int secondDimSize = array.shape(1);
    vector<vector<double>> result(firstDimSize);
    for (int i = 0; i < firstDimSize; i++) {
        result[i] = vector<double>(secondDimSize);
        for (int j = 0; j < secondDimSize; j++) {
            result[i][j] = ((double *) (array.get_data() + i * firstDimStrides + j * secondDimStrides))[0];
        }
    }
    return result;
}

void reportDouble2dToNumpy(vector<vector<double>> src, np::ndarray array) {
    const int firstDimStrides = array.strides(0);
    const int secondDimStrides = array.strides(1);
    const int firstDimSize = array.shape(0);
    const int secondDimSize = array.shape(1);
    for (int i = 0; i < firstDimSize; i++) {
        for (int j = 0; j < secondDimSize; j++) {
            ((double *) (array.get_data() + i * firstDimStrides + j * secondDimStrides))[0] = src[i][j];
        }
    }
}

void check1dParams(np::ndarray stacks, np::ndarray payouts,
                   np::ndarray output) {
    const int stacksNd = stacks.get_nd();
    if (stacksNd != 1) {
        throw invalid_argument(
                "Expected first argument stacks to be a 1d array with shape (max(nb_players, nb_payouts))");
    }
    const int payoutsNd = payouts.get_nd();
    if (payoutsNd != 1) {
        throw invalid_argument(
                "Expected second argument payouts to be a 1d array with shape (max(nb_players, nb_payouts))");
    }
    const int outputNd = output.get_nd();
    if (outputNd != 1) {
        throw invalid_argument(
                "Expected third argument payouts to be a 1d array with shape (max(nb_players, nb_payouts))");
    }
    const int stacksNbPlayers = stacks.shape(0);
    const int payoutsNbPlayers = payouts.shape(0);
    const int outputNbPlayers = output.shape(0);
    if (stacksNbPlayers != payoutsNbPlayers || stacksNbPlayers != outputNbPlayers) {
        throw invalid_argument("Stacks, payouts and output don't have the same size. Add zeros if needed.");
    }
}

void check1dParams(np::ndarray stacks, np::ndarray output) {
    const int stacksNd = stacks.get_nd();
    if (stacksNd != 1) {
        throw invalid_argument(
                "Expected first argument stacks to be a 1d array with shape (nb_players)");
    }

    const int outputNd = output.get_nd();
    if (outputNd != 1) {
        throw invalid_argument(
                "Expected third argument payouts to be a 1d array with shape (nb_players)");
    }
    const int stacksNbPlayers = stacks.shape(0);
    const int outputNbPlayers = output.shape(0);
    if (stacksNbPlayers != outputNbPlayers) {
        throw invalid_argument("Stacks and output don't have the same size. Add zeros if needed.");
    }
}

void monteCarloIcmNumpy(np::ndarray stacks, np::ndarray payouts, np::ndarray output, long trials) {
    check1dParams(stacks, payouts, output);
    const auto vectorStacks = numpyToDouble(stacks);
    const auto vectorPayouts = numpyToDouble(payouts);
    vector<double> result = icm::monteCarloIcm(vectorStacks, vectorPayouts, trials);
    reportDoubleToNumpy(result, output);
}

void monteCarloIcmMultiThreadNumpy(np::ndarray stacks, np::ndarray payouts,
                                   np::ndarray output, long trials, int nbThreads) {
    const int stacksNd = stacks.get_nd();
    if (stacksNd != 2) {
        throw invalid_argument(
                "Expected first argument stacks to be a 2d array with shape (nb_samples, max(nb_players, nb_payouts))");
    }
    const int payoutsNd = payouts.get_nd();
    if (payoutsNd != 2) {
        throw invalid_argument(
                "Expected second argument payouts to be a 2d array with shape (nb_samples, max(nb_players, nb_payouts))");
    }
    const int outputNd = output.get_nd();
    if (outputNd != 2) {
        throw invalid_argument(
                "Expected third argument output to be a 2d array with shape (nb_samples, max(nb_players, nb_payouts))");
    }
    const int stacksNbSamples = stacks.shape(0);
    const int payoutsNbSamples = payouts.shape(0);
    const int outputNbSamples = output.shape(0);
    if (payoutsNbSamples != stacksNbSamples || outputNbSamples != stacksNbSamples) {
        throw invalid_argument(
                str(boost::format(
                        "Expected to have the same first dimension size for first three arguments (stacks, payouts, output) but got %1%, %2% and %3%") %
                    stacksNbSamples % payoutsNbSamples % outputNbSamples));
    }
    const int stacksNbPlayers = stacks.shape(1);
    const int payoutsNbPlayers = payouts.shape(1);
    const int outputNbPlayers = output.shape(1);
    if (payoutsNbPlayers != stacksNbPlayers || outputNbPlayers != stacksNbPlayers) {
        throw invalid_argument(
                str(boost::format(
                        "Expected to have the same second dimension size for first three arguments (stacks, payouts, output) but got %1%, %2% and %3%\nYou can add zeros if needed.") %
                    stacksNbPlayers % payoutsNbPlayers % outputNbPlayers));
    }
    auto vectorStacks = numpyToDouble2d(stacks);
    auto vectorPayouts = numpyToDouble2d(payouts);
    auto result = icm::monteCarloIcm(vectorStacks, vectorPayouts, trials, nbThreads);
    reportDouble2dToNumpy(result, output);
}

void exactIcmNumpy(np::ndarray stacks, np::ndarray payouts, np::ndarray output) {
    check1dParams(stacks, payouts, output);
    const auto vectorStacks = numpyToDouble(stacks);
    const auto vectorPayouts = numpyToDouble(payouts);
    vector<double> result = icm::icm(vectorStacks,
                                     vectorPayouts);
    reportDoubleToNumpy(result, output);
}

void monteCarloIcmWithStoppingRuleNumpy(np::ndarray stacks, np::ndarray payouts,
                                             np::ndarray results, double stoppingAlphaProbability,
                                             double stoppingConfidenceIntervalLength, long long stoppingEvalLag,
                                             // Working around my specific setup issue (SIGSEGV when returning any scalar)
                                             // https://stackoverflow.com/questions/66985776/boost-python-basic-function-returning-int-sigsegv
                                             np::ndarray nbIterationsResultContainer) {
    check1dParams(stacks, payouts, results);
    const auto vectorStacks = numpyToDouble(stacks);
    const auto vectorPayouts = numpyToDouble(payouts);
    auto vectorResults = std::make_shared<vector<double>>();
    auto nbTrials = icm::monteCarloIcm(vectorStacks, vectorPayouts, vectorResults,
                                       stoppingAlphaProbability,
                                       stoppingConfidenceIntervalLength, stoppingEvalLag);
    reportDoubleToNumpy(*vectorResults, results);
    *((long long *)nbIterationsResultContainer.get_data()) = nbTrials;
}

// Simple test functions for this issue : https://stackoverflow.com/questions/66985776/boost-python-basic-function-returning-int-sigsegv
//int testFunction() {
//    std::cout << "Test Function 1\n";
//    return 1;
//}
//
//void testFunction2() {
//    std::cout << "Test Function 2\n";
//}

BOOST_PYTHON_MODULE (icm) {
    np::initialize();
    py::def("monte_carlo_icm", monteCarloIcmNumpy);
    py::def("monte_carlo_icm_multithread", monteCarloIcmMultiThreadNumpy);
    py::def("exact_icm", exactIcmNumpy);
    py::def("monte_carlo_icm_with_stopping_rule", monteCarloIcmWithStoppingRuleNumpy);
//    py::def("test_function", testFunction);
//    py::def("test_function_2", testFunction2);
}
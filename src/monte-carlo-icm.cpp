#include "../include/mtt-icm/monte-carlo-icm.hpp"
#include "monte-carlo-icm-util.hpp"
#include <random>
#import <boost/asio.hpp>
#include <sfmt.hpp>
#include "math.hpp"

using namespace std;

namespace icm {

    /**
     * Low level utilities
     */
    namespace mc_util {
        void monteCarloPermutation(const double *weights, int *destination, const int size, const RNG &rng,
                                   const Sorter &sorter) {
            double draw[size];
            for (int i = 0; i < size; i++) {
                draw[i] = log2(rng()) * weights[i];
                destination[i] = i;
            }
            sorter(draw, destination, size);
        }

        int firstZeroPayout(const double *payouts, int size) {
            int firstZeroPayout = -1;
            for (int i = 0; i < size; ++i) {
                if (firstZeroPayout < 0) {
                    if (payouts[i] == 0) {
                        firstZeroPayout = i;
                    }
                } else if (payouts[i] != 0) {
                    firstZeroPayout = -1;
                }
            }
            return firstZeroPayout;
        }

        vector<double> mean(vector<double *> samples, const int nbPlayers) {
            const long long nbSamples = samples.size();
            vector<double> result(nbPlayers, 0);
            for (int i = 0; i < nbSamples; i++) {
                double *sample = samples[i];
                for (int j = 0; j < nbPlayers; ++j) {
                    result[j] += sample[j];
                }
            }
            for (int i = 0; i < nbPlayers; ++i) {
                result[i] /= (double) nbSamples;
            }
            return result;
        }

        double variance(vector<double *> samples, const int nbPlayers) {
            const long long nbSamples = samples.size();
            vector<double> meanSample = mean(samples, nbPlayers);
            double result = 0;
            for (int i = 0; i < nbSamples; i++) {
                double *sample = samples[i];
                for (int j = 0; j < nbPlayers; ++j) {
                    double diff = sample[j] - meanSample[j];
                    result += diff * diff;
                }
            }
            // Usually, the result should be divided by nbSamples and not nbSamples - 1, but I wrote the latest
            // to respect the paper's stopping rule algorithm ( http://www.lib.ncsu.edu/resolver/1840.4/5244 )
            return result / (double) (nbSamples - 1);
        }

        double stoppingRuleLimit(double alphaProbability, double confidenceIntervalLength, long long nbSamples) {
            double zAlpha = NormalCDFInverse(alphaProbability + (1 - alphaProbability) / 2); // Symmetrical distribution
            return (double)nbSamples * confidenceIntervalLength * confidenceIntervalLength / zAlpha;
        }

        Sorter defaultSorter() {
            return [](const double *draw, int *destination, int size) {
                sort(destination, destination + size,
                     [&](const int &a, const int &b) {
                         // Reverse order
                         return (draw[a] > draw[b]);
                     });
            };
        }

        Sorter partialSorter(int relevantRanksCount) {
            return [=](const double *draw, int *destination, int size) {
                partial_sort(destination, destination + relevantRanksCount, destination + size,
                             [&](const int &a, const int &b) {
                                 // Reverse order
                                 return (draw[a] > draw[b]);
                             }
                );
            };
        }

        Sorter chooseSorter(const double *payouts, int nbPlayers) {
            // Partial sort vs sort : interesting only if we have many zero payouts at the end
            const int relevantRanksCount = firstZeroPayout(payouts, nbPlayers);
            const float partialSortThreshold = 0.2; // Magic number <3
            const bool partialSort = (float) relevantRanksCount / (float) nbPlayers <= partialSortThreshold;
            if (partialSort) {
                return partialSorter(relevantRanksCount);
            }
            return defaultSorter();
        }

        RNG defaultRng() {
            // Initialize random distribution
            // https://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
            // SFMT variant
            random_device rd;
            wtl::sfmt19937 mt(rd());
            shared_ptr<wtl::sfmt19937> mtPtr = std::make_shared<wtl::sfmt19937>(mt);
            uniform_real_distribution<double> dist(0, 1);
            return [mtPtr = make_shared<wtl::sfmt19937>(mt),
                    distPtr = make_shared<uniform_real_distribution<double>>(dist)]() -> double {
                return (*distPtr)((*mtPtr));
            };
        }

        bool canStop(const double alphaProbability, const double confidenceIntervalLength, const vector<double *> &samples,
                     const int nbPlayers) {
            const long long nbSamples = samples.size();
            if (nbSamples < 2) {
                return false;
            }
            double limit = stoppingRuleLimit(alphaProbability, confidenceIntervalLength, nbSamples);
            return variance(samples, nbPlayers) <= limit;
        }

        void fillWeights(const double *stacks, int nbPlayers, double *weights) {
            const double stacksAvg = accumulate(stacks, stacks + nbPlayers, 0.) / (double) nbPlayers;
            for (int i = 0; i < nbPlayers; i++) {
                weights[i] = stacksAvg / stacks[i];
            }
        }
    }

    /**
     * Perform one MC trial
     *
     * @param results destination array to add trial contribution into
     * @param payouts the contribution to add to the results per rank
     * @param weights something proportional to the inverse of the stacks of the players
     * @param permutation an array that will be used to store the permutation (permutation[i] = player that reaches rank i)
     * @param nbPlayers the number of players
     * @param rng random number generator
     * @param sorter function to sort the permutation
     */
    void monteCarloIcmTrial(double *results, const double *payouts, double *weights, int *permutation, int nbPlayers,
                            const mc_util::RNG &rng, const mc_util::Sorter &sorter) {
        mc_util::monteCarloPermutation(weights, permutation, nbPlayers, rng, sorter);
        // Cumulate payouts according to the random permutation
        for (int j = 0; j < nbPlayers; j++) {
            results[permutation[j]] += payouts[j];
        }
    }

    /**
     * Monte-Carlo ICM Ranking algorithm.
     *
     * @param stacks the stacks
     * @param payouts the payouts
     * @param nbPlayers the number of players
     * @param trials the number of trials to perform
     * @param results destination array for ICM EV
     */
    void monteCarloIcm(const double *stacks, const double *payouts, const int nbPlayers, const long long trials,
                       double *results) {
        // Algorithm : https://forumserver.twoplustwo.com/15/poker-theory/new-algorithm-calculate-icm-large-tournaments-1098489/
        mc_util::RNG rng = mc_util::defaultRng();
        int permutation[nbPlayers];

        double weights[nbPlayers];
        mc_util::fillWeights(stacks, nbPlayers, weights);

        double contrib[nbPlayers];
        // Prepare each trial ranking contribution to the total
        for (int i = 0; i < nbPlayers; i++) {
            contrib[i] = payouts[i] / (double) trials;
        }

        const mc_util::Sorter sorter = mc_util::chooseSorter(payouts, nbPlayers);
        for (long i = 0; i < trials; i++) {
            monteCarloIcmTrial(results, contrib, weights, permutation, nbPlayers, rng, sorter);
        }
    }

    /**
     * Monte-Carlo ICM Ranking algorithm with stopping rule.
     *
     * @param stacks the stacks
     * @param payouts the payouts
     * @param nbPlayers the number of players
     * @param results destination array for ICM EV
     * @param stoppingAlphaProbability the confidence probability
     * @param stoppingConfidenceIntervalLength confidence interval length (root mean square)
     * @param stoppingEvalLag the number of trials that are performed before each evaluation of the stopping rule
     */
    long long monteCarloIcm(const double *stacks, const double *payouts, const int nbPlayers,
                            double *results, const double stoppingAlphaProbability,
                            const double stoppingConfidenceIntervalLength, long long stoppingEvalLag) {
        // Algorithm : https://forumserver.twoplustwo.com/15/poker-theory/new-algorithm-calculate-icm-large-tournaments-1098489/
        // Stopping rule : http://www.lib.ncsu.edu/resolver/1840.4/5244 (first method for independent samples)
        mc_util::RNG rng = mc_util::defaultRng();
        int permutation[nbPlayers];

        double weights[nbPlayers];
        mc_util::fillWeights(stacks, nbPlayers, weights);

        // For each trial we'll contribute the real payouts because we want to use them for the stopping rule
        // computation

        const mc_util::Sorter sorter = mc_util::chooseSorter(payouts, nbPlayers);

        vector<double *> samples;
        vector<double *> allocated;
        long long count = 0;
        while (true) {
            // While we don't stop according to the stopping rule, we create a batch of samples
            auto *memory = (double *) calloc(nbPlayers * stoppingEvalLag, sizeof(double));
            // Store the allocated memory area pointer
            allocated.push_back(memory);
            for (long i = 0; i < stoppingEvalLag; i++) {
                double *sampleMemory = memory + i * nbPlayers;
                monteCarloIcmTrial(sampleMemory, const_cast<double *>(payouts), weights, permutation, nbPlayers, rng,
                                   sorter);
                // Store the pointer to the sample
                samples.push_back(sampleMemory);
            }
            count += stoppingEvalLag;
            // Break if the stopping rule says it's enough
            if (mc_util::canStop(stoppingAlphaProbability, stoppingConfidenceIntervalLength, samples, nbPlayers)) {
                break;
            }
        }
        // Sum the samples
        for (double *sample : samples) {
            for (int i = 0; i < nbPlayers; i++) {
                results[i] += sample[i];
            }
        }
        // Divide the sum to get the mean
        for (int i = 0; i < nbPlayers; i++) {
            results[i] /= (double) count;
        }
        // Free the allocated memory areas
        for (double *memory : allocated) {
            free(memory);
        }
        // Return the number of samples that were generated
        return count;
    }

    long long monteCarloIcm(vector<double> stacks, vector<double> payouts,
                            shared_ptr<vector<double>> results, double stoppingAlphaProbability,
                            double stoppingConfidenceIntervalLength, long long stoppingEvalLag) {
        const int nbPlayers = stacks.size();
        double stacksArray[nbPlayers];
        double payoutsArray[nbPlayers];
        double resultsArray[nbPlayers];
        for (int i = 0; i < nbPlayers; ++i) {
            stacksArray[i] = stacks[i];
            payoutsArray[i] = payouts[i];
        }
        long long nbTrials = monteCarloIcm(stacksArray, payoutsArray, nbPlayers, resultsArray,
                                           stoppingAlphaProbability, stoppingConfidenceIntervalLength,
                                           stoppingEvalLag);
        auto resultsVector = vector<double>(nbPlayers);
        for (int i = 0; i < nbPlayers; i++) {
            resultsVector[i] = resultsArray[i];
        }
        *results = resultsVector;
        return nbTrials;
    }

    /**
     * Vector wrapper of MC ICM
     * @param stacks
     * @param payouts
     * @param trials
     * @return ICM EV
     */
    vector<double> monteCarloIcm(const vector<double> &stacks, const vector<double> &payouts, const long trials) {
        const int nbPlayers = stacks.size();
        vector<double> results(nbPlayers, 0);
        double stacksArray[nbPlayers];
        double payoutsArray[nbPlayers];
        double resultsArray[nbPlayers];
        for (int i = 0; i < nbPlayers; ++i) {
            stacksArray[i] = stacks[i];
            payoutsArray[i] = payouts[i];
            resultsArray[i] = 0;
        }
        monteCarloIcm(stacksArray, payoutsArray, nbPlayers, trials, resultsArray);
        for (int i = 0; i < nbPlayers; i++) {
            results[i] = resultsArray[i];
        }
        return results;
    }

    /**
     * Multithread implementation for batch execution
     * @param stacks the stacks
     * @param payouts the payouts
     * @param trials the number of trials
     * @param nbThreads the number of threads to use
     * @return a vector of ICM EV vectors
     */
    vector<vector<double>> monteCarloIcm(vector<vector<double>> &stacks, vector<vector<double>> &payouts,
                                         long trials, int nbThreads) {
        boost::asio::static_thread_pool pool(nbThreads);
        auto exec = pool.executor();
        auto alloc = allocator<void>();

        const int nbSituations = stacks.size();
        vector<vector<double>> result(nbSituations);
        for (int i = 0; i < nbSituations; i++) {
            exec.post([=, &result]() {
                result[i] = monteCarloIcm(stacks[i], payouts[i], trials);
            }, alloc);
        }
        pool.join();
        pool.stop();
        return result;
    }
}
#ifndef MTT_ICM_MONTE_CARLO_ICM_UTIL_HPP
#define MTT_ICM_MONTE_CARLO_ICM_UTIL_HPP

/**
 * Low level utilities
 */
namespace icm::mc_util {

    /**
    * The Random Number Generator type : functions denoted by this alias are expected to draw doubles between 0 and 1.
    * Abstracted in order to plug a proper RNG for quasi-Monte-Carlo in the future.
    */
    using RNG = std::function<double()>;
    /**
     * The sorter functions are expected to sort the indexes of the draw array in reverse order according to the values
     * of the array.
     */
    using Sorter = std::function<void(double *draw, int *destination, int size)>;

    /**
     * Compute a permutation according to the sorting of random values put to the power of given weights.
     * In fact, logarithm is used for performance, and logarithm tables could certainly be used to fasten the
     * computation.
     *
     * @param weights
     * @param destination where the permutation will be written to
     * @param size size of weights and destination array
     * @param rng random number generator
     * @param sorter function to sort the permutation
     */
    void monteCarloPermutation(const double *weights, int *destination, const int size, const RNG &rng,
                               const Sorter &sorter);

    /**
     * Get the index of the first zero payout after which all payouts are zero as well
     * @param payouts the payouts array
     * @param size the size of the array
     * @return the index of the first zero payout
     */
    int firstZeroPayout(const double *payouts, int size);

    /**
     * Compute the mean of a batch of samples
     * @param samples the samples
     * @param nbPlayers the size of each sample aka the number of players
     * @return a vector containing the mean of the samples
     */
    vector<double> mean(vector<double *> samples, const int nbPlayers);

    /**
     * Compute the variance of samples
     * @param samples the samples
     * @param nbPlayers the size of each sample aka the number of players
     * @return the variance
     */
    double variance(vector<double *> samples, const int nbPlayers);

    /**
     * Compute the threshold for the variance under which the stopping limit is considered reached
     * @param alphaProbability confidence probability
     * @param confidenceIntervalLength confidence interval length (root mean square)
     * @param nbSamples number of samples
     * @return the variance threshold
     */
    double stoppingRuleLimit(double alphaProbability, double confidenceIntervalLength, long long nbSamples);

    /**
     * Default sorting using std::sort
     * @return the default sorter
     */
    Sorter defaultSorter();

    /**
     * Sorter using partial sort. Relevant when relevant ranks count is low relatively to the number of players
     * @param relevantRanksCount the number of top ranks that have a non-null price
     * @return a partial sorter
     */
    Sorter partialSorter(int relevantRanksCount);

    /**
     * Pick the appropriate sorter depending on payouts. We can use a partial sort if most of the payouts are zero
     * @param payouts the payouts
     * @param nbPlayers the size of the payouts array aka the number of players
     * @return a sorter
     */
    Sorter chooseSorter(const double *payouts, int nbPlayers);

    /**
     * Create a RNG instance (using SFMT)
     * @return a RNG instance
     */
    RNG defaultRng();

    /**
     * Check if we can stop the trials according to the confidence probability and interval
     * @param alphaProbability the confidence probability
     * @param confidenceIntervalLength confidence interval length (root mean square)
     * @param samples the samples
     * @param nbPlayers the number of players
     * @return true if the stopping rule criteria are met
     */
    bool canStop(const double alphaProbability, const double confidenceIntervalLength, const vector<double *> &samples,
                 const int nbPlayers);

    /**
     * Normalize stacks (not really necessary but cheap) and invert to weights
     * @param stacks the stacks
     * @param nbPlayers the number of players
     * @param weights the destination array
     */
    void fillWeights(const double *stacks, int nbPlayers, double *weights);
}
#endif //MTT_ICM_MONTE_CARLO_ICM_UTIL_HPP

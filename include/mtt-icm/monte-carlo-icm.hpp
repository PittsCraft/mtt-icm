#ifndef MTT_ICM_MONTE_CARLO_ICM_HPP
#define MTT_ICM_MONTE_CARLO_ICM_HPP
#import "vector"
using namespace std;

/**
 * Monte-Carlo ICM estimations.
 * See https://pittscraft.com/posts/poker_mtt_icm_monte_carlo/
 * for detailed explanations
 */
namespace icm {

    /**
     * Vector version of Monte-Carlo ICM algorithm
     * @param stacks the stacks
     * @param payouts the payouts
     * @param trials the number of trials to run
     * @return the ICM estimation
     */
    vector<double> monteCarloIcm(const vector<double> &stacks, const vector<double> &payouts, long trials);

    /**
     * Raw Monte-Carlo ICM with stopping rule algorithm.
     * Stopping rule : http://www.lib.ncsu.edu/resolver/1840.4/5244 (first method for independent samples)
     *
     * @param stacks the stacks (expected size = number of players)
     * @param payouts the payouts (expected size = number of players)
     * @param nbPlayers the number of players
     * @param results the destination pointer where the ICM estimation will be written (expected size = number of players)
     * @param stoppingAlphaProbability the target probability
     * @param stoppingConfidenceIntervalLength the target confidence interval length
     * @param stoppingEvalLag how many trials to run before evaluating the stopping conditions
     * @return the number of trials that were run before reaching the stopping conditions
     */
    long long monteCarloIcm(const double *stacks, const double *payouts, int nbPlayers,
                            double *results, double stoppingAlphaProbability,
                            double stoppingConfidenceIntervalLength, long long stoppingEvalLag);

    /**
     * Vector interface of Monte-Carlo ICM with stopping rule algorithm
     * Stopping rule : http://www.lib.ncsu.edu/resolver/1840.4/5244 (first method for independent samples)
     *
     * @param stacks the stacks (expected size = number of players)
     * @param payouts the payouts (expected size = number of players)
     * @param nbPlayers the number of players
     * @param results the destination pointer where the ICM estimation will be written (expected to be empty)
     * @param stoppingAlphaProbability the target probability
     * @param stoppingConfidenceIntervalLength the target confidence interval length
     * @param stoppingEvalLag how many trials to run before evaluating the stopping conditions
     * @return the number of trials that were run before reaching the stopping conditions
     */
    long long monteCarloIcm(vector<double> stacks, vector<double> payouts,
                            shared_ptr<vector<double>> results, double stoppingAlphaProbability,
                            double stoppingConfidenceIntervalLength, long long stoppingEvalLag);

    /**
     * Batch multi-thread Monte-Carlo ICM computation
     * @param stacks the list of stacks
     * @param payouts the list of payouts
     * @param trials the number of trials to apply for each situation
     * @param nbThreads the number of threads to use
     * @return the list of ICM estimations
     */
    vector<vector<double>> monteCarloIcm(vector<vector<double>> &stacks, vector<vector<double>> &payouts,
                                         long trials, int nbThreads);

}
#endif //MTT_ICM_MONTE_CARLO_ICM_HPP

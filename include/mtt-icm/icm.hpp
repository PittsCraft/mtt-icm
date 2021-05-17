#ifndef MTT_ICM_ICM_HPP
#define MTT_ICM_ICM_HPP
using namespace std;

namespace icm {

    /**
     * Exact ICM algorithm
     * See https://pittscraft.com/posts/poker_mtt_icm_calcul/ (french) for explanations
     *
     * @param stacks the stacks
     * @param payouts the payouts. Expected to be sorted in descending order, or at least having the first zero payout followed only by zero payouts
     * @return the ICM given the stacks and payouts
     */
    vector<double> icm(const vector<double> &stacks, const vector<double> &payouts);

}
#endif //MTT_ICM_ICM_HPP
#include <random>
#include <iostream>

using namespace std;

namespace icm {
    namespace {

        void recursiveIcm(const double *stacks, const double *payouts,
                          double *result,
                          const long long usedPlayers, const int rank,
                          const int size, const double factor, const double stacksSum,
                          double *cacheValues, bool *cacheFlags) {
            if (rank == size || payouts[rank] == 0) {
                return;
            }
            // Point subResult to the cache entry
            double *subResult = cacheValues + usedPlayers * size;
            // If the subResult wasn't computed previously, let's do it
            if (!*(cacheFlags + usedPlayers)) {
                for (int i = 0; i < size; i++) {
                    const long long iMask = 0x1ll << i;
                    if (usedPlayers & iMask) { continue; }
                    // Probability of ith player to have this rank *among remaining players*
                    // (ignoring previously ranked ones)
                    const double newFactor = stacks[i] / stacksSum;
                    // Fill the EV of this subresult for this player
                    subResult[i] += newFactor * payouts[rank];
                    recursiveIcm(stacks, payouts, subResult, usedPlayers | iMask, rank + 1,
                                 size, newFactor, stacksSum - stacks[i], cacheValues, cacheFlags);
                }
                // Mark the cache entry as filled
                *(cacheFlags + usedPlayers) = true;
            }
            // Apply the factor accounting for previously ranked players and add to the destination array
            for (int i = 0; i < size; ++i) {
                result[i] += factor * subResult[i];
            }
        }
    }

    vector<double> icm(const vector<double> &stacks, const vector<double> &payouts) {
        const int size = stacks.size();
        // Prepare the array that will receive the results from the recursive function
        double result[size];
        fill_n(result, size, 0);
        // Just translate vectors to C arrays
        double stacksArray[size];
        double payoutsArray[size];
        for (int i = 0; i < size; ++i) {
            stacksArray[i] = stacks[i];
            payoutsArray[i] = payouts[i];
        }
        // Compute the total stack sum
        const double stacksSum = accumulate(stacks.begin(), stacks.end(), 0.);
        // ## Cache ##
        const auto cacheEntrySize = size * sizeof(double); // One double for each player
        // We'll store all cached subresult in this memory area.
        auto cacheValues = (double *) calloc(pow(2, size), cacheEntrySize);
        // And we'll store the information of whether each one is filled or not.
        auto cacheFlags = (bool *) calloc(pow(2, size), sizeof(bool));

        // Ok let's go
        recursiveIcm(stacksArray, payoutsArray, result, 0, 0, size, 1, stacksSum,
                     cacheValues, cacheFlags);
        // Just translate back to vector
        vector<double> toReturn(size);
        for (int i = 0; i < size; ++i) {
            toReturn[i] = result[i];
        }
        // And free the allocated memory
        free(cacheValues);
        free(cacheFlags);
        return toReturn;
    }
}
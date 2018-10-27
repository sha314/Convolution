//
// Created by shahnoor on 2/1/2018.
//

#include "../include/binomial.h"

using namespace std;


/**
 *
 * @param size
 * @param center
 * @param probability
 * @return
 */
vector<double> binomial_distribution_v1(size_t size, size_t center, double probability) {
    vector<double> binomial(size);
    double sm{};
    binomial[center] = 1;
    double factor;
    for (size_t j{center + 1}; j != size; ++j) {  // for j > center
        factor = (size - j + 1) * probability / (j * (1 - probability));
        binomial[j] = binomial[j - 1] * factor;
    }

    for (size_t j{center - 1}; j != -1; --j) {  // for j < center
        factor = (j + 1) * (1 - probability) / ((size - j) * probability);
        binomial[j] = binomial[j + 1] * factor;
    }

    for (size_t j{}; j != size; ++j) {
        sm += binomial[j];
    }
    for (size_t j{}; j != size; ++j) {
        binomial[j] /= sm;  // normalizing
    }
    return binomial;
}

/**
 *
 * @param size
 * @param center
 * @param probability
 * @return
 */
vector<double> binomial_distribution_v2(size_t size, size_t center, double probability) {
    vector<double> binomial(size);
    binomial[center] = 1;
    double factor;
    double normalization_const={1};
    for (size_t j{center + 1}; j != size; ++j) {  // for j > center
        factor = (size - j + 1) * probability / (j * (1 - probability));
        binomial[j] = binomial[j - 1] * factor;
        normalization_const += binomial[j];
    }

    for (size_t j{center - 1}; j != -1; --j) {  // for j < center
        factor = (j + 1) * (1 - probability) / ((size - j) * probability);
        binomial[j] = binomial[j + 1] * factor;
        normalization_const += binomial[j];
    }

    for (size_t j{}; j != size; ++j) {
        binomial[j] /= normalization_const;  // normalizing
    }
    return binomial;
}



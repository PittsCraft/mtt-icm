#ifndef MTT_ICM_MATH_HPP
#define MTT_ICM_MATH_HPP

/**
 * Compute the inverse of the normal CDF (cumulative distribution function)
 * See https://www.johndcook.com/blog/cpp_phi_inverse/
 *
 * @param p the target normal distribution cumulated value
 * @return the argument for which the normal CDF will return p
 */
double NormalCDFInverse(double p);

#endif //MTT_ICM_MATH_HPP

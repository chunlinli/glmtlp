/**********
    C++ utility functions for GLMTLP.

    Copyright (C) 2021-2022 Chunlin Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
**********/

#pragma once

// #ifndef UTILS_HPP
// #define UTILS_HPP

//#include "glmtlp.hpp"


#include <vector>
#include <string.h>
#include <string>
#include <cmath>
#include <queue>

#include "glmtlp.h"

inline double soft_thresh(double init, double thresh)
{
    if (init > thresh)
        init -= thresh;
    else if (init < -thresh)
        init += thresh;
    else
        init = 0.0;
    return init;
}

inline double link(double mu, int family)
{
    switch (family)
    {
    case 1: // Family::Gaussian:
        return mu;
    case 2: // Family::Binomial:
        return std::log(mu / (1.0 - mu));
    case 3: // Family::Poisson:
        return std::log(mu);
    default:
        return NAN;
    }
}

inline double compute_deviance(const Eigen::VectorXd &y,
                               const Eigen::VectorXd &eta,
                               const Eigen::VectorXd &w,
                               int family)
{
    switch (family)
    {
    case 1: // Family::Gaussian:
    {
        return (y - eta).array().square().matrix().dot(w);
    }
    case 2: // Family::Binomial:
    {
        Eigen::ArrayXd mu = 1.0 / (1.0 + exp(-eta.array()));
        return -2.0 * w.dot((log(mu) * y.array() + log(1.0 - mu) * (1.0 - y.array())).matrix());
    }
    case 3: // Family::Poisson:
    {
        Eigen::ArrayXd mu = exp(eta.array());
        return 2.0 * (w.array() * (log(y.array() / mu + 0.000000001) * y.array())).sum();
    }

    default:
    {
        return NAN;
    }
    }
}

inline void check_user_interrupt()
{
    Rcpp::checkUserInterrupt();
}

inline void glmtlp_warning(const std::string &msg)
{
    Rcpp::warning("[GLMTLP] " + msg);
}


//typedef Eigen::Triplet<double> T;


// optimizer

// template

// void check_user_interrupt();
// void glmtlp_warning(const std::string& msg);
//
// inline double soft_thresh(double init, double thresh);
// inline double link(double mu, int family);
// inline double compute_deviance(const Eigen::VectorXd &y,
//                                const Eigen::VectorXd &eta,
//                                const Eigen::VectorXd &w,
//                                int family);

// double objective(Eigen::VectorXd &v, double bound, double lambda)
// {

//     double obj = 0.0;
//     for (int j = 0; j < v.size(); ++j)
//     {
//         double change = std::abs(v(j)) - lambda;
//         if (change > 0)
//         {
//             obj += change;
//         }
//     }
//     return obj - bound;
// }

//#endif // UTILS_HPP

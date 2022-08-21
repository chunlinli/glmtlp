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

#include "macros.h"
#include "glmtlp.h"

enum class Family
{
    Gaussian,
    Binomial,
    Poisson,
    Noncanonical
};

enum class Method
{
    Lasso,
    RTLP,
    CTLP
};

STRONG_INLINE
double soft_thresh(double init, double thresh)
{
    if (init > thresh)
        init -= thresh;
    else if (init < -thresh)
        init += thresh;
    else
        init = 0.0;
    return init;
}

STRONG_INLINE
double link(double mu, Family family)
{
    switch (family)
    {
    case Family::Gaussian:
        return mu;
    case Family::Binomial:
        return std::log(mu / (1.0 - mu));
    case Family::Poisson:
        return std::log(mu);
    default:
        return NAN;
    }
}

STRONG_INLINE
double compute_deviance(const Eigen::VectorXd &y,
                               const Eigen::VectorXd &eta,
                               const Eigen::VectorXd &w,
                               Family family)
{
    switch (family)
    {
    case Family::Gaussian:
    {
        return (y - eta).array().square().matrix().dot(w);
    }
    case Family::Binomial:
    {
        Eigen::ArrayXd mu = 1.0 / (1.0 + exp(-eta.array()));
        return -2.0 * w.dot((log(mu) * y.array() + log(1.0 - mu) * (1.0 - y.array())).matrix());
    }
    case Family::Poisson:
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

STRONG_INLINE
void check_user_interrupt()
{
    Rcpp::checkUserInterrupt();
}

STRONG_INLINE
void glmtlp_warning(const std::string &msg)
{
    Rcpp::warning("[GLMTLP] " + msg);
}


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

#ifndef UTILS_HPP
#define UTILS_HPP

//#include "glmtlp.hpp"

#include <vector>
#include <string.h>
#include <string>
#include <cmath>
#include <queue>


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

#endif // UTILS_HPP

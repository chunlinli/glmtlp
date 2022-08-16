/**********
    R Interface.

    Copyright (C) 2020-2021  Chunlin Li

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

#include "glmtlp.h"
#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"

extern "C"
{

    SEXP gaussian_l1(SEXP b0,
                     SEXP b,
                     SEXP r,
                     SEXP X,
                     SEXP w,
                     SEXP rho,
                     SEXP lambda,
                     SEXP nlambda,
                     SEXP n,
                     SEXP p,
                     SEXP delta,
                     SEXP tol,
                     SEXP cd_maxit)
    {
        linreg_l1_ssr(REAL(b0),
                      REAL(b),
                      REAL(r),
                      REAL(X),
                      REAL(w),
                      REAL(rho),
                      REAL(lambda),
                      *INTEGER(nlambda),
                      *INTEGER(n),
                      *INTEGER(p),
                      *REAL(delta),
                      *REAL(tol),
                      *INTEGER(cd_maxit));

        return R_NilValue;
    }

    SEXP gaussian_tlp(SEXP b0,
                      SEXP b,
                      SEXP r,
                      SEXP X,
                      SEXP w,
                      SEXP rho,
                      SEXP lambda,
                      SEXP nlambda,
                      SEXP tau,
                      SEXP n,
                      SEXP p,
                      SEXP delta,
                      SEXP tol,
                      SEXP cd_maxit,
                      SEXP dc_maxit)
    {
        linreg_tlp_ssr(REAL(b0),
                       REAL(b),
                       REAL(r),
                       REAL(X),
                       REAL(w),
                       REAL(rho),
                       REAL(lambda),
                       *INTEGER(nlambda),
                       *REAL(tau),
                       *INTEGER(n),
                       *INTEGER(p),
                       *REAL(delta),
                       *REAL(tol),
                       *INTEGER(cd_maxit),
                       *INTEGER(dc_maxit));

        return R_NilValue;
    }

    SEXP gaussian_l0(SEXP b0,
                     SEXP b,
                     SEXP r,
                     SEXP X,
                     SEXP w,
                     SEXP rho,
                     SEXP s,
                     SEXP ns,
                     SEXP lambda,
                     SEXP nlambda,
                     SEXP tau,
                     SEXP n,
                     SEXP p,
                     SEXP delta,
                     SEXP tol,
                     SEXP cd_maxit,
                     SEXP dc_maxit)
    {
        linreg_l0_ssr(REAL(b0),
                      REAL(b),
                      REAL(r),
                      REAL(X),
                      REAL(w),
                      REAL(rho),
                      INTEGER(s),
                      *INTEGER(ns),
                      REAL(lambda),
                      *INTEGER(nlambda),
                      *REAL(tau),
                      *INTEGER(n),
                      *INTEGER(p),
                      *REAL(delta),
                      *REAL(tol),
                      *INTEGER(cd_maxit),
                      *INTEGER(dc_maxit));

        return R_NilValue;
    }

    SEXP logistic_l1(SEXP b0,
                     SEXP b,
                     SEXP y,
                     SEXP X,
                     SEXP rho,
                     SEXP lambda,
                     SEXP nlambda,
                     SEXP n,
                     SEXP p,
                     SEXP delta,
                     SEXP tol,
                     SEXP nr_maxit,
                     SEXP cd_maxit)
    {
        logistic_l1_ssr(REAL(b0),
                        REAL(b),
                        REAL(y),
                        REAL(X),
                        REAL(rho),
                        REAL(lambda),
                        *INTEGER(nlambda),
                        *INTEGER(n),
                        *INTEGER(p),
                        *REAL(delta),
                        *REAL(tol),
                        *INTEGER(nr_maxit),
                        *INTEGER(cd_maxit));

        return R_NilValue;
    }

} // extern "C"

static const R_CallMethodDef CallEntries[] = {
    {"gaussian_l0", (DL_FUNC)&gaussian_l0, 17},
    {"gaussian_l1", (DL_FUNC)&gaussian_l1, 13},
    {"gaussian_tlp", (DL_FUNC)&gaussian_tlp, 15},
    {"logistic_l1", (DL_FUNC)&logistic_l1, 13},
    {NULL, NULL, 0}};

void R_init_glmtlp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

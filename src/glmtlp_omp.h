// this file is adapted from "biglasso_omp.h"
// see https://github.com/YaohuiZeng/biglasso

#ifndef GLMTLP_OMP_H_
#define GLMTLP_OMP_H_
#if defined(_OPENMP)
#include <omp.h>
#else
#ifndef DISABLE_OPENMP
#pragma message("OpenMP is not available... Use single thread...")
#endif
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_num_procs() { return 1; }
inline void omp_set_num_threads(int nthread) {}
inline void omp_set_dynamic(int flag) {}
#endif
#endif // GLMTLP_OMP_H_
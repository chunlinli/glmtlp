# Version 3 Test Results

## Unsolved Issues
- [ ] why the R processes still running after being cut off.
  - temporary solution: run `killall -9 R` to manually cut those R processes.

## Solved Issues
- [x] the very slow progress when p=10000: replace Sigma-based with autoregression, and only output data$X and data$y.


## Result Description and Analysis

- tests01: `test01.csv`, `log01.txt`
    - this test has p=10000 inside and was stopped by Yu due to its slow progress at p=10000.
    - another issue is that after stopping the process via `Ctrl+C`, the R processes are still running. 
- tests02: `test02.csv`, `log02.txt`
    - this test has p=10000 inside and was stopped by Yu due to its slow progress at p=10000.
    - the purpose of this test is to check if the slow speed was due to the machine. Turns out no, it was the problem of p being too large, not the machine's issue.
    - the R processes running issue still exists in this trial.
- tests03: `test03.csv`, `fig03`, `log03.txt`
    - nrep = 10
    - Performance analysis based on the figures
        - in terms of obj, no big difference for all the 12 cases below.
        - in terms of runtime
            1. n=100, p=10: no big difference except logistic has larger sd at large snr level.
            2. n=100, p=100: ncvreg and logistic are similar and glmnet is the smallest, all within 0.15 seconds
            3. n=100, p=1000: ncvreg and logistic are similar and glmnet is the smallest, all within 0.1 seconds
            4. n=200, p=10: logistic < glmnet < ncvreg
            5. n=200, p=100: logistic is the largest, while lasso and ncvreg are close, all within the level of 7 seconds.
            6. n=200, p=1000: ncvreg and logistic are similar and glmnet is the smallest, all within 0.15 seconds
            7. n=500, p=10: logistic < glmnet < ncvreg
            8. n=500, p=100: logistic has larger sd and mean at large snr level.
            9. n=500, p=1000: ncvreg and logistic are similar and glmnet is the smallest, all within 0.8 seconds
            10. n=1000, p=10: logistic < glmnet < ncvreg, all within 0.1 seconds.
            11. n=1000, p=100: logistic and glmnet are comparable for small correlations, for large correlation, glmnet < logistic < ncvreg
            12. n=1000, p=1000: ncvreg and logistic are similar and glmnet is the smallest, all within 4 seconds
- tests07: `test07.csv`, `fig07`, `log07.txt`
  - nrep = 10
  - p = 10, 100, 1000, 10000, 30000
  - this test intends to check the large-scale performance.
  - this test uses the new data generation method: instead of Sigma and Cholesky, autoregression is used to generate X.
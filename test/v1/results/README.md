# Tests Version 1 Results

This folder contains the results for version 1 tests, targeting at `lasso()`.

- `tests.csv`: the test results for `tol = sum((y - mean(y))^2)/2*(1e-7)`, run on Yu's MacBook Pro.
- `tests2.csv`: the test results for `tol = 1e-5`, run on Yu's MacBook Pro.
- `tests3.csv`: the test results for `tol = 1e-4`, run on the DAGS server. (parameter settings have been changed.)
- `tests4.csv`: the test results for `tol = 1e-5`, run on the DAGS server. 
- `tests5.csv`: the test results for `tol = 1e-4`, run on the DAGS server, using the package built for the server. The seed is set as 8053. 
- `tests6.csv`: the test results for `tol = 1e-4`, run on the DAGS server, using the package built for the server. The seed is set as 5451. The goal is to check whether the result is similar to that in Test 5.
- `tests10.csv`: the test result for `tol = 1e-4` (glmtlp), `thresh = 1e-8` (glmnet), and `eps = 1e-4` (ncvreg), run on DAGS server. 
- `tests11.csv`: the same as `tests10.csv`, different seed. 
- `test_analysis.ipynb`: the jupyter notebook to analyze the test results.


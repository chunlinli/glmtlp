# Version 2 Tests

This is the 2nd version of tests, targeting at `cv.l0reg()` and `cv.tlpreg()`, in comparison with SCAD, MCP, and LASSO.


## Scripts
- `test2.R`: the script to systematically compare l0reg, tlpreg, scad, mcp, and lasso.
- `testutils2.R`: the utility functions used in `test2.R`.
- `plot.R`: the script to plot the metrics comparison.
- `extra_metrics.R`: 
    1. rerun the data generation process with the same seed to log down two more measures: "sigma2","risk.null" for computing "risk.rel","err.rel"; 
    2. calculate "risk.rel","err.rel" based on `tests02.csv` and generate a new result file `tests02_appended.csv`.
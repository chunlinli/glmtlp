# Version 2 Test Results

- tests01: `test01.csv`, `fig01`, `log01.txt`
    - nrep = 1
    - `fig01`: the metrics comparison plots for tests01.
- tests02: `test02.csv`, `fig02`, `log02.txt`, `tests02_appended.csv`, `log_extra_metrics.txt`
    - nrep = 5
    - the logging system has been updated, which can now print out parameter settings and the timestamps.
    - `tests02.csv`: only contains 10 measures: "runtime","err","err.test","prop","risk","nzs","fp","fn","F1","opt".
    - `extra_metrics.csv`: the results from `extra_metrics.R`, rerunning the data generation process with the same seed to log down two more measures: "sigma2","risk.null" for computing "risk.rel","err.rel".
    - `tests02_appended.csv`: the results after appending the two more measures "risk.rel","err.rel".
    - `log_extra_metrics.txt`: the log file for `extra_metrics.R`.
- tests03: `test03.csv`, `fig03`, `log03.txt`
    - nrep = 10
    - include more measures ("tp", "risk.rel","err.rel","hd","hd.rel"), now in total 15 measures: "runtime","err","err.test","prop","risk","nzs","tp","fp","fn","F1","opt","risk.rel","err.rel","hd","hd.rel".
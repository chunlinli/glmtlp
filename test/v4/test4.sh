module load R/3.6.1
echo "Running test for logistic ..."
R CMD BATCH test3.R results/test_3.7.Rout
echo "Test finished"

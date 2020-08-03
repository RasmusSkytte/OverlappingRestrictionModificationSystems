# Run the data generators
seq 1 325 | parallel -j 50 "matlab -nodesktop -nosplash -r 'cd code/functions; generateData({}); exit'"

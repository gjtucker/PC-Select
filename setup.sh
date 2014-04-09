mkdir working
mkdir data

matlab -nojvm -nodesktop -nosplash -r "cd src; mex fast_write.cpp; exit;"

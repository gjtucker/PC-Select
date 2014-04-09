#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdio.h>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[])
{
	double* a;
	a = mxGetPr(prhs[1]);
	
	const mwSize* dims = mxGetDimensions(prhs[1]);
	int r = (int) dims[0];
	int c = (int) dims[1];
	
	char* input_buf;
	size_t buflen;

	buflen = mxGetN(prhs[0]) + 1;
	input_buf = mxArrayToString(prhs[0]);

	FILE* f = fopen(input_buf, "w");
	for (int j = 0; j < c; j++) {
		for (int i = 0; i < r; i++) {
			fprintf(f, "%d\t", (int)*(a++));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	mxFree(input_buf);
}

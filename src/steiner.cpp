/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// steiner.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include "spgsolver.h"
//#include "labels.h"

void PickSubset (int n, int k, int s, int *count=NULL) {
	if (!count) fprintf (stderr, "Picking a subset of size %d from %d using seed %d.\n", k, n, s);
	RFWLocalRandom random(s);
	int left = k;
	for (int i=1; i<=n; i++) {
		if (random.GetInteger(0,n-i)<left) {
			left --;
			if (count) {count[i]++;}
			else fprintf (stderr, "%3d: %3d\n", k - left, i); 
		}
	}

	//fprintf (stderr, ".");

	if (!count) {

		int *vcount = new int [n+1];
		for (int i=0; i<=n; i++) vcount[i] = 0;
		int tries = 1000000;
		for (s=1; s<=tries; s++) {
			PickSubset(n,k,s,vcount);
		}
		for (int i=1; i<=n; i++) {
			fprintf (stderr, "%.2f ", (double)n * (double) vcount[i] / ((double)tries * (double)k)); 
		}
		fprintf (stderr, "\n");
		delete [] vcount;
	}


}


int main (int argc, char **argv) {
//int _tmain(int argc, _TCHAR* argv[])
//{
	bool PICK_SUBSET = false;
	if (PICK_SUBSET) {
		PickSubset(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
		exit(-1);
	}

	SPGSolver::Solve(argc, argv);

	return 0;
}


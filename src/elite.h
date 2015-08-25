/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once
#include "solution.h"
#include "rfw_random.h"
#include "spgconfig.h"
#include <cstdio>

class SolutionPool {
private:
	vector<SteinerSolution*> sol;
	int count; //number of solutions currently in the pool
	int capacity; //capacity of the pool
	SteinerConfig *config;
	
	
public:
    SolutionPool(int cap) {
		capacity = cap;
		sol.resize(cap+1);
        //sol = new Solution[capacity + 1];
        for (int i = 0; i <= cap; i++) sol[i] = NULL;
        count = 0;
	}

	inline int GetCapacity() {
		return capacity;
	}

	inline int GetCount() {
		return count;
	}



	/// <summary>
    /// Tries adding a solution to the pool. Returns its new position if successful,
    /// and zero if the solution is not added.
    /// </summary>
    /// <param name="x"></param>
    /// <returns></returns>
    int Add(SteinerSolution *s) {
		bool verbose = false;
		const bool VERBOSE_DIFF = false;

        double bestdiff = 0;

		RFWLocalRandom random;
        //OptRandom random = new OptRandom();

		if (verbose) fprintf(stderr, "Adding solution.\n");
		bool present = false; //is the solution already in the pool
        int wcount = 0; //number of solutions worse than s
        int bestpos = 0;
        bool DETERMINISTIC = false;
		bool SQUARE_SCORE = true;
        double accscore = 0;

        for (int i = 1; i <= count; i++) {
            //any solution that is equal or worse than s is a candidate for removal
            if (!sol[i]->IsBetter(s)) {
				double diff = s->GetDifference(sol[i]); //1 - inter/union (of edges)
				if (diff == 0) {
					present = true;
					break;
				}
				if (VERBOSE_DIFF) fprintf (stderr, " %.3f:%.0f", diff, sol[i]->GetCost());

				wcount ++;

				if (DETERMINISTIC) {
					//Console.Error.Write("[{0}] ", diff);
                    if (bestpos == 0 || diff < bestdiff) {
						bestpos = i;
                        bestdiff = diff;
					} else {
						if (diff == bestdiff) {
							fprintf (stderr, "%.2f ", diff);
                        }
                    } 
				} else {
					double curscore = 1.0 / (double)diff; //(sol[i].GetValue() - s.GetValue());
					if (SQUARE_SCORE) {curscore *= curscore;}
					accscore += curscore;
                    if (bestpos == 0 || (random.GetDouble() * accscore < curscore)) {
						if (VERBOSE_DIFF) fprintf (stderr, "<<<");
						bestpos = i;
                        bestdiff = diff;
					}
				}
			}
		}

		if (VERBOSE_DIFF) fprintf (stderr, "\n");

		//solution is already there
        if (present) return 0;

        if (verbose) fprintf (stderr, "count%d:cap%d ", count, capacity);

		//pool is full and nobody is worse
        if (count == capacity) { //full pool?
			if (wcount == 0) return 0; //nobody worse
		} else { //there's still space
			count++; //add one more element in the last position
            bestpos = count;
		}

        if (verbose) {
			fprintf (stderr, "+");
            fprintf (stderr, "Adding to position %d.\n", bestpos);
		}

        if (sol[bestpos] == NULL) {
			sol[bestpos] = new SteinerSolution(s); //(Solution)s.Clone();
		} else {
			if (verbose) fprintf (stderr, "Replacing %d by %d difference is %.4f, out of %d candidates.", sol[bestpos]->GetCost(), s->GetCost(), bestdiff, wcount);
			sol[bestpos]->CopyFrom(s);
		}
		return bestpos;
	}

    void Pick(int i, SteinerSolution *s) {
		//if (i <= 0 || i > count) throw new Exception("Invalid solution found");
        s->CopyFrom(sol[i]);
	}

	SteinerSolution *GetReference(int i) {
		return sol[i];
	}

	int FindBestPosition () {
		fprintf (stderr, "Here I am (%d,%d)!\n", count, GetCount());
		fflush(stderr);

		if (count == 0) return -1;
		int best = 1;
		for (int i=2; i<=count; i++) {
			fprintf (stderr, "%d ", i);
			fflush(stderr);
			if (sol[i] == NULL) exit(-1);
			if (sol[i]->GetCost() < sol[best]->GetCost()) best = i;
		}

		fprintf (stderr, "Returning %d.\n", best);
		fflush(stderr);

		return best;
	}

	EdgeCost FindBestCost() {
		if (count == 0) return numeric_limits<EdgeCost>::max();
		else return sol[FindBestPosition()]->GetCost();
	}

	void Output(FILE *file, int columns) {
		for (int i=1; i<=count; i++) {
			int column = (i-1) % columns;
			if (column == 0) {
				if (i > 1) fprintf (file, "\n");
				fprintf (file, "%5d :", i);
			}
			fprintf (stderr, " %6lg", sol[i] ? (double)sol[i]->GetCost() : -1.0);
		}
		fprintf (stderr, "\n");
	}


	void Reset() { count = 0; }

	void HardReset() {
		for (int i=0; i<=capacity; i++) {
			if (sol[i]) {
				delete sol[i];
				sol[i] = NULL;
			}
		}
		count = 0;
	}

    int Capacity() { return capacity; }
    int Count() { return count; }

	~SolutionPool() {
		for (int i=0; i<=capacity; i++) {
			if (sol[i]) delete sol[i];
		}
	}
};
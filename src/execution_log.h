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

#include <vector>
#include <omp.h>

using namespace std;

#include "graph.h"
#include "rfw_timer.h"
#include "solution.h"
#include "graph.h"

class ExecutionLog {

public:

	ExecutionLog(Graph *g, RFWTimer *t, double tl) : timerPtr(t), bestSolution(g), timeLimit(tl) {}

	typedef pair<EdgeCost, double> EntryType;

	inline void AddSolution(SteinerSolution &sol) {
		// Obey the time limit.
		if (timeLimit > 0 && timerPtr->getTime() > timeLimit)
			return;

#pragma omp critical (exec_log) 
		{
			if (solCost.empty()) {
				solCost.push_back(EntryType(sol.GetCost(), timerPtr->getTime()));
				bestSolution.CopyFrom(&sol);
			}
			else {
				if (solCost.back().first > sol.GetCost() + EDGE_COST_PRECISION) {
					solCost.push_back(EntryType(sol.GetCost(), timerPtr->getTime()));
					bestSolution.CopyFrom(&sol);
				}
			}
		}
	}

	// These are the cost and running times for each incumbent solution.
	vector<EntryType> solCost;

	// Pointer to a "global" timer.
	RFWTimer *timerPtr = nullptr;

	// Everything is awesome.
	SteinerSolution bestSolution;
	const double timeLimit;

};
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

#include <cstdio>
#include <cstdlib>
#include "uset.h"
#include "graph.h"
#include "simplequeue.h"

class DualAscentData {
	int hidden_vertex; //vertex we'll pretend is not in the graph (useful for strong branching)

public:
	Graph *g;
	UniverseSet *rlist;
    SimpleQueue *queue;
    UniverseSet *component;
	EdgeCost *rcost;
	EdgeCost *dist_from_root;
    EdgeCost *dist_to_terminal; 
	int *unsatlist; //
	int *unsattails; //tails of unsatlist (we could just ask the graph, but keeping an explicit list improves locality)
	int unsatcount; //


    bool *delarc; //delarc[a] is true iff there if (directed) *ARC* can be deleted
	int deledgecount;
	int n;
	int m;
	EdgeCost dualcost; //total dual cost of the current solution
	int root;

    inline bool IsHidden(int v) {return (v ==hidden_vertex);}
	inline void Hide(int v) {hidden_vertex = v;}
	inline void UnhideVertices() {hidden_vertex = 0;}

	DualAscentData (Graph *_g) {
		//fprintf (stderr, "Initializing ascent data.\n");
		g = _g;
        n = g->VertexCount();
        m = g->EdgeCount();
		rlist = new UniverseSet(n);
		queue = new SimpleQueue(n);
        component = new UniverseSet(n);
        rcost = new EdgeCost [2*m+1];
        delarc = new bool [2*m+1];
		unsatlist = new int [2*m];
		unsattails = new int [2*m];
		unsatcount = 0;
		deledgecount = 0;
        UnhideVertices();

		dist_from_root = new EdgeCost [n+1];
		dist_to_terminal = new EdgeCost [n+1];


		ResetReducedCosts();
        for (int e=0; e<=2*m; e++) {delarc[e] = false;}
		dualcost = 0;
	}	

	/*
	bool IsArcDeleted (SPGArc *a) const {
		return deledge[a->label];
	}*/

	~DualAscentData() {
		delete [] dist_to_terminal;
		delete [] dist_from_root;
		delete [] unsatlist;
		delete [] unsattails;
		delete [] delarc;
		delete [] rcost;
		delete component;
		delete queue;
		delete rlist;
	}


	// for 
	void ResetReducedCosts() {
		for (int v=1; v<=n; v++) {
			SPGArc *a, *end;
			//for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
			//int alabel = g.GetOutgoingLabel(v, pa);
			//rcost[alabel] = g.GetArcCost(pa);//a.cost;
            //        }
			for (g->GetBounds(v,a,end); a<end; a++) {
				int alabel = g->GetOutgoingLabel(v,a);
				rcost[alabel] = a->cost;
			}
        }
	}
};

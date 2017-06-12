#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include "binheap.h"
#include "graph.h"
#include "solution.h"
#include "rfw_timer.h"
#include "LSBasics.h"

/*
#include "uf.h"
#include "uset.h"
#include "voronoi.h"
#include "rfw_random.h"
#include "rfw_stack.h"
#include "drawer.h"
#include "dual.h"
#include "stedgelinear.h"
#include "pairheap.h"
#include "buckets.h"
#include <cstring>
#include <cmath>
#include "elite.h"
#include "spgconfig.h"
#include <omp.h>
#include "execution_log.h"
#include "LSVertexInsertion.h"
#include "LSVertexElimination.h"
#include "LSBasics.h"
#include "LSKeyPath.h"
#include "BranchBound.h"
*/

class ConstructiveAlgorithms {
private:
	static void fatal (const string &msg) {
		Basics::fatal(msg);
	}


public:
	static void SPH (Graph &g, SteinerSolution &solution, EdgeCost *pertcost, int root) {
		const bool verbose = false;
        unsigned int v;
        int n = g.VertexCount();
        int m = g.EdgeCount();

		static bool PATH_ERROR = false;

        if (verbose) fprintf (stderr, "Running SPH from %d... ", root);

        BinaryHeap<EdgeCost> heap(n); //vertices to scan
		vector <int> parc(n+1,-1); //parent arc
		vector <EdgeCost> dist(n+1,-1); //current distance

        solution.Reset();
        heap.Insert(root, 0);
        dist[root] = 0; //root already in the tree, but no parent edge
        parc[root] = 0; //already in the tree

        int nterm = g.TerminalCount();
        int tcount = 0;
        if (g.IsTerminal(root)) tcount ++;

		//int scancount = 0;
		//int inscount = 1;
		//int pathsum = 0;

        //while (!heap.IsEmpty()) {
        while (tcount < nterm) {
			if (verbose) {fprintf (stderr, "Here."); fflush(stderr);}
			if (heap.IsEmpty()) {
				for (int x=1; x<=n; x++) {
					if (dist[x] < 0) fprintf (stderr, "Vertex %d has distance %d.\n", x, dist[x]);
				}
				fatal ("heap is empty");
			}
			EdgeCost curcost;
			heap.RemoveFirst(v, curcost); //lowest element in the heap
			if (verbose) {
				fprintf (stderr, "Processing %d with distance %d.\n", v, curcost);
				fflush (stderr);
			}
			if (curcost<0) {
				fprintf (stderr, "Vertex %d has negative cost (%d)", v, curcost);
				fatal ("invalid heap entry");
			}

			//found a terminal: add incoming path to the solution,
			//reinsert its vertices into the heap with zero distance
            if (g.IsTerminal(v) && (parc[v]>0)) {
				int curpath = 0;
				if (verbose) fprintf (stderr, "Adding terminal %d.\n", v);
				//fflush(stderr);
				tcount ++;
                int w = v;
                EdgeCost addedcost = 0;
                do {
					int e = parc[w];
					//if (e <= 0) break;
					curpath ++;
                    if (v!=w) heap.Insert(w, 0); // RECENT CHANGE: NO NEED TO INSERT IF V==W
					//else fprintf (stderr, ".");
					//inscount ++;
                    dist[w] = 0; //its distance label is zero
                    parc[w] = 0; //w now belongs to the tree
                    if (!solution.Insert(e)) fatal ("failed to add edge");
					EdgeCost acost = (pertcost!=NULL) ? pertcost[e] : g.GetCost(e);
                    addedcost += acost;
                    w = g.GetOther(e,w);
				} while (parc[w]>0);
				if (addedcost != curcost) {
					if (!PATH_ERROR) {
						fprintf (stderr, "WARNING (SPH): addedcost:%.20lg curcost:%.20lg\n", addedcost, curcost);
						PATH_ERROR = true;
					}
					//fprintf (stderr, "Inconsistent path in SPH.");
				}
				//pathsum += curpath;
				curcost = 0;
			}

			if (verbose) fprintf (stderr, "Processing a non-terminal.\n");

			// scan outgoing arcs, add anything that improves to the solution
			//scancount ++;
			SPGArc *arc, *end;
			for (g.GetBounds(v,arc,end); arc<end; arc++) {
				int w = arc->head;
				if (parc[w]==0) continue; //already in the solution				
				//if ((w==root) || (solution.GetDegree(w)>0)) continue; //already in solution
                EdgeCost acost = (pertcost) ? pertcost[arc->label] : arc->cost; //g.GetCost(arc->label);                    
                EdgeCost newcost = acost + curcost;

                if (newcost < 0) {
					fprintf (stderr, "Arc %d + cost %lg led to newcost %lg.", acost, (double)curcost, (double)newcost);
					fatal ("Cost overflow.");
				}

				if ((newcost<dist[w]) || (dist[w]==-1)) {
					dist[w] = newcost;
                    if (!heap.Insert(w,newcost)) fatal ("insertion failed");
					//inscount ++;
					parc[w] = arc->label;
					if (verbose) {fprintf (stderr, "Inserted %d with cost %d.\n", w, newcost); fflush (stderr);}
				}
			}
			if (verbose) {fprintf (stderr, "x"); fflush(stderr);}
		}
		//fprintf (stderr, "%.3f ", (double)scancount / (double)n);
		//fprintf (stderr, "%.3f ", (double)pathsum / (double)tcount);
		//fprintf (stderr, "%.3f,%.3f seconds ", (double)inscount / (double)n, timer.getTime()*1000);
		if (verbose) fprintf (stderr, "Found solution %d ", solution.GetCost());
	}


};



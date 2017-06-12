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
#include "uf.h"
#include "uset.h"
//#include "voronoi.h"
#include "rfw_random.h"
#include "rfw_stack.h"
//#include "drawer.h"
#include "dual.h"
//#include "stedgelinear.h"
//#include "pairheap.h"
//#include "buckets.h"
//#include <cstring>
//#include <cmath>
//#include "elite.h"
//#include "spgconfig.h"
//#include <omp.h>
#include "execution_log.h"
#include "constructive.h"
//#include "LSVertexInsertion.h"
//#include "LSVertexElimination.h"
//#include "LSBasics.h"
//#include "LSKeyPath.h"
//#include "BranchBound.h"
#undef INFINITY


class BranchBound {
private:
	static void fatal (const string &msg) {
		Basics::fatal(msg);
	}


public:

	class BBNodeInfo {
		public:
			int lastbranch; //most recent branch 0, 1, or other
			int depth;
			int splits; 
			int laststrict; //last strict branch (0 or 1)

			BBNodeInfo() {
				lastbranch = 2;
				laststrict = 2;
				depth = 0;
				splits = 0;
			}

			BBNodeInfo (BBNodeInfo &info, int b) {
				lastbranch = b;
				bool s = (b>=0 && b<=1); //is this actually only half of the parent?
				laststrict = s ? b : info.laststrict;
				depth = info.depth + 1;
				splits = info.splits + (s ? 1 : 0);
			}

	};

	class BBStats {
	public:
		double lstime;     //local search
		double sphtime;    //running the constructive algorithm
		double dualtime;   //running the dual ascent
		double fixtime;    //fixing nodes by reduced cost
		double branchtime; //branching
		double subtime;    //building the subproblem
		RFWTimer timer;

		double nodecount[3];   //total number of nodes in the graph
		double branchcount[3]; //number of branching nodes that are children of 0, 1, 2
		int maxdepth;
		int maxrealdepth; 
		double fraction_solved;
		GlobalInfo *ginfo;
		
		BBStats() {
			lstime = 0;
			sphtime = 0;
			dualtime = 0;
			fixtime = 0;
			branchtime = 0;
			subtime = 0;
			nodecount[0] = nodecount[1] = nodecount[2] = 0;
			branchcount[0] = branchcount[1] = branchcount[2] = 0;
			maxdepth = 0;
			maxrealdepth = 0;
			fraction_solved = 0;
			timer.start();
			ginfo = NULL;
		}

		void MarkAsSolved (BBNodeInfo &bbinfo) {
			double fraction = 1;
			for (int i=0; i<bbinfo.splits; i++) {
				fraction /= 2;
			}
			fraction_solved += fraction;
		}

		void Output (FILE *file, bool abridged=false) {
			double totaltime = timer.getTime();
			double projtime = totaltime / fraction_solved;
			if (!abridged) {
				fprintf (file, "lstimems %.3f\n", 1000 * lstime);
				fprintf (file, "sphtimems %.3f\n", 1000 * sphtime);
				fprintf (file, "dualtimems %.3f\n", 1000 * dualtime);
				fprintf (file, "fixtimems %.3f\n", 1000 * fixtime);
				fprintf (file, "branchtimems %.3f\n", 1000 * branchtime);
				fprintf (file, "subtimems %.3f\n", 1000 * subtime);

				double othertime = totaltime - lstime - sphtime - dualtime - fixtime - branchtime - subtime;
				fprintf (file, "othertimems %.3f\n", 1000 * othertime);
				fprintf (file, "totaltimems %.3f\n", 1000 * totaltime);
				for (int i=0; i<=2; i++) {
					fprintf (file, "bbcount%d %.0f\n", i, nodecount[i]);
				}
				for (int i=0; i<=2; i++) {
					fprintf (file, "branchcount%d %.0f\n", i, branchcount[i]);
				}
				fprintf (file, "branches %.0f\n", branchcount[0]+branchcount[1]+branchcount[2]);

				double skew = 0.5;
				if (branchcount[0] + branchcount[1] > 0) {
					skew = branchcount[1] / (branchcount[0] + branchcount[1]);
				}
				fprintf (file, "skew %.3f\n", skew);
				fprintf (file, "maxrealdepth %d\n", maxrealdepth);
				fprintf (file, "projseconds %.2f\n", projtime);
				fprintf (file, "projhours %.6f\n", projtime / 3600.0);
				fprintf (file, "projbbnodes %.0f\n", (nodecount[0]+nodecount[1]+nodecount[2]) / fraction_solved);
			}



			fprintf (file, "bbcounttotal %.0f\n", nodecount[0]+nodecount[1]+nodecount[2]);
			fprintf (file, "maxdepth %d\n", maxdepth);
			fprintf (file, "totaltime %.3f\n", totaltime);
			fprintf (file, "pctsolved %.6f\n", 100*fraction_solved); 
			fprintf (file, "projdays %.6f\n", projtime / 86400.0);
			fprintf (file, "projbbnodesm %.6f\n", (nodecount[0]+nodecount[1]+nodecount[2]) / (1000000.0 * fraction_solved));
		}

		double GetTotalNodes() {
			return nodecount[0] + nodecount[1] + nodecount[2];
		}
	};


	static void RunBranchAndBound(Graph &g, int seed, EdgeCost &primal, EdgeCost bestfound, EdgeCost bestknown, SteinerConfig &config, GlobalInfo *ginfo, ExecutionLog *executionLogPtr) {
		fprintf (stderr, "RUNNING BRANCH AND BOUND WITH FIXED %.3f\n", g.GetFixedCost());
		fflush(stderr);
		RFWTimer timer(true);
		int n = g.VertexCount();
		int m = g.EdgeCount();
		vector<int> vfix (n+1, -1);
		vector<int> efix (m+1, -1);
		RFWRandom::randomize(21);
		int root = 0; 
		RFWLocalRandom random(seed);
		if (seed > 0) {
			root = Basics::PickRandomTerminal(g,random); 
			//root = PickBestRoot(g, primal, seed, 200);
		}
		//root = 48;
		fprintf (stderr, "Picked %d as root, with degree %d.\n", root, g.GetDegree(root));
		BBStats bbstats;
		bbstats.ginfo = ginfo;

		fprintf (stderr, "WARNING: IGNORING FIXED COST FOR BB.\n");
		EdgeCost fixed = g.GetFixedCost(); //hack!
		//if (fixed > EDGE_COST_PRECISION) fatal ("BB cannot handle fixed costs properly yet.\n");
		g.SetFixedCost(0);
		primal -= fixed; //WARNING: THIS DOESN'T SEEM RIGHT
		vector<double> priority;
		vector<int> ekeep (g.EdgeCount()+1, 0); // optimistically assume that we get to keep all the edges
		//ScatterBranch(priority, g, 2000, 32, seed, root, primal);
		primal = bestfound;
		int count = BBound(g, vfix, efix, primal, 0, seed, root, &bbstats, BBNodeInfo(), false, config, priority, ekeep.data(), executionLogPtr);
		primal += fixed;
		g.SetFixedCost(fixed);

		double times = timer.getTime();
		fprintf (stderr, "\nAfter %d nodes, solution is %.0f (%.2f ms).\n", count, (double)primal, 1000.0 * timer.getTime());
		fprintf (stderr, "Partial times (ms): const:%.2f ls:%.2f dual:%.2f\n", 1000 * bbstats.sphtime, 1000*bbstats.lstime, 1000*bbstats.dualtime);	
		fprintf (stdout, "seed %d\n", seed);
		bbstats.Output(stdout);
		fprintf (stdout, "bbnodes %d\n", count);

		/*
		fprintf (stdout, "primal %.0f\n", (double)primal);
		fprintf (stdout, "bbtimeus %.3f\n", 1000000.0 * times);
		fprintf (stdout, "bbtimems %.6f\n", 1000.0 * times);
		fprintf (stdout, "bbtimes %.9f\n", times);
		double ratio = (double)primal / (double)bestknown;
		double error = ratio - 1;
		fprintf (stdout, "bbratio %.12f\n", ratio);
		fprintf (stdout, "bberror %.12f\n", error);
		fprintf (stdout, "bbpcterror %.12f\n", 100.0 * error);
		*/
		//EdgeCost actualprimal = primal + (ginfo ? ginfo->fixed : 0);
		Basics::ReportResults(stdout, "bb", times, primal, bestknown);
	}


	// WARNING: WHY DOES THIS RETURN DOUBLE?
	static EdgeCost DualAscent(Graph &g, EdgeCost &primal, vector<int> &efix, GraphMapper *mapper, int &branch, int hide, int seed, int root, BBStats *bbstats, double *priority, bool RUN_PRIMAL, bool FIND_BRANCHING, SteinerConfig *config, ExecutionLog *executionLogPtr) {
		
		//fprintf (stderr, "R%d ", seed);
		RFWLocalRandom random (seed); //((int)mapper); //that's not reproducible...

		

		const bool USE_LOCAL_SEARCH = true;
		const bool FIX_REDUCED_COST = true;
		bool verbose = false;
		branch = -1;
		double best = 0; //why double?
		int bestcount = 0;
		double totalvalue = 0;
		double totalfixed = 0;
		RFWTimer timer(true);
		//fprintf (stderr, "%d ", primal);
		int ASCENT_TRIES = 1;

		//int bestcurprimal = -1; //best primal found in the current round
		//int bestcurbranch = -1;

		for (int i=0; i<ASCENT_TRIES; i++) {
			DualAscentData data(&g); 
			if (hide>0) data.Hide(hide); //used for strong branching
	        //ScoreDualAscent(data, 0, 0, false, random);
			//FastDualAscent(data, 0, 0, false, random);

			//fprintf (stderr, "\n");
			bool TEST_CONNECTEDNESS = false;
			if (TEST_CONNECTEDNESS) {
				SteinerSolution solution(data.g);
				ConstructiveAlgorithms::SPH (g, solution, NULL, Basics::PickRandomTerminal(g,random));
				EdgeCost initcost = solution.GetCost();
				fprintf (stderr, "<%d> ", initcost);
			}


			RFWTimer fdatime(true);
			FastDualAscent(data, root, primal, true, random, config);
			//fprintf (stderr, "DFDA"); fflush(stderr);
			if (bbstats) bbstats->dualtime+=fdatime.getTime();
			double curvalue = (double)data.dualcost;
			totalvalue += curvalue;
			if (curvalue == best) bestcount ++;
			else if (curvalue > best) {
				bestcount = 1;
				best = curvalue;
			}

			// if the dual is better than the best known primal, done
			if (curvalue >= primal - EDGE_COST_PRECISION) {
				//fprintf (stderr, "!");
				break;
			}


			SteinerSolution solution(data.g); //THIS IS TEMPORARY
			if (RUN_PRIMAL) {
				//fprintf (stderr, "RP"); fflush(stderr);
				bool verbose = false;
				bool debug = false;
				//EdgeCost newundirected = 0; //SPHDual (data, solution, true);
				RFWTimer sphtimer(true);
				EdgeCost newprimal = SPHDual (data, solution, false);
				if (bbstats) bbstats->sphtime += sphtimer.getTime();
				//fprintf (stderr, "dual:%.3f d:%.3f pctgap:%.3f ", curvalue, newprimal, 100.0 * (newprimal / curvalue - 1.0)); fflush(stderr);

				//fprintf (stderr, "%d ", newprimal);

				if (USE_LOCAL_SEARCH && newprimal < INFINITE_COST - EDGE_COST_PRECISION) {
					//fprintf (stderr, "Running local search.\n");
					//CheckSolution(g, solution);
					int itcount = 0;

					EdgeCost oldcost = solution.GetCost();

					int maxit = 1;
					if (maxit > 0) {
						//SteinerSolution localsol(data.g);
						RFWTimer lstimer(true);

						for (int i=0; i<maxit; i++) {
						//while (1) {
							itcount ++;
							LSVertexInsertion::VertexInsertion (g, solution, 5*g.VertexCount(), random);
							//MSTPrune(g,solution);
							LSVertexElimination::VertexElimination(g, solution, random);
							//KeyVertexElimination(g, solution, random);
							//MSTPrune(g,solution);
							EdgeCost newcost = solution.GetCost();
							if (newcost > oldcost + EDGE_COST_PRECISION) fatal ("Something wrong.");
							if (newcost > oldcost - EDGE_COST_PRECISION) break;
							oldcost = newcost;
						}
						if (bbstats) bbstats->lstime += lstimer.getTime();
					}

					//fprintf (stderr, "%d ", itcount);

					bool CHECK_LOCAL_OPTIMA = false;
					if (CHECK_LOCAL_OPTIMA) {				
						oldcost = solution.GetCost();
						int moves = LSVertexElimination::VertexElimination(g, solution, random);
						if (solution.GetCost() < oldcost) {
							//fprintf (stderr, "Supposed local minimum is %d, vertex elimination found %d (itcount=%d, moves=%d).\n", oldcost, solution.GetCost(), itcount, moves);	
							fprintf (stderr, "#");
							//fatal ("inconsistent local search");
						} else {
							fprintf (stderr, "_");
						}
					}
				
					//VertexInsertion (g, solution, 5*g.VertexCount(), random);
					//VertexElimination(g, solution, random);
					if (debug) Basics::CheckSolution(g, solution);
					EdgeCost lsprimal = solution.GetCost();
					if (lsprimal > newprimal + EDGE_COST_PRECISION) {
						fprintf (stderr, "Solution changed from %d to %d.\n", newprimal, lsprimal);
						fatal ("invalid local search");
					}
					if (verbose) {
						if (lsprimal < newprimal - EDGE_COST_PRECISION) fprintf (stderr, "#");
						else fprintf (stderr, ".");
						fflush(stderr);
					}
					//fprintf (stderr, "i%d ", newprimal - lsprimal);
					//fprintf (stderr, "improvement:%d (%d->%d)\n", newprimal - lsprimal, newprimal, lsprimal);
					newprimal = lsprimal;
				}
				//fprintf (stderr, "Here!\n");

				//if (bestcurprimal < 0 || newprimal < bestcurprimal) {bestcurprimal = newprimal;}




				//SPH(*data.g, solution, NULL, data.root);
				//newprimal = solution.GetCost();
				//fprintf (stderr, "<PP%.3f>", newprimal);

				if (newprimal < primal - EDGE_COST_PRECISION) {
					primal = newprimal;
					if (executionLogPtr != nullptr) {
						fprintf(stderr, "ADDING SOLUTION FROM BB WITH VALUE %f.\n", primal);
						executionLogPtr->AddSolution(solution);
					}
					const bool OUTPUT_PRIMAL = true;
					if (OUTPUT_PRIMAL) {
						fprintf (stderr, "p%.3f ", primal);
						fflush(stderr);
					}
				}
			}

			//fprintf(stderr, "FRC"); fflush(stderr);

			// fix edges by reduced cost; this is cumulative, and 
			// can be done across multiple calls to B&B

			double totaltime = timer.getTime();
			RFWTimer fixtimer(true);
			int fcount = 0;
			if (curvalue < primal) {
				//fprintf (stderr, "!"); fflush(stderr);
				if (FIX_REDUCED_COST) {fcount = FixReducedCost(data, primal);}
				if (mapper) { // must mark the original edge ids
					int oldm = mapper->oldm;
					for (int e=1; e<=oldm; e++) {
						int f = mapper->e2new[e];
						if (f>0 && data.delarc[f]) {efix[e] = 0;} //CAREFUL: TREATING ARC AS EDGE HERE
					}
				}
			}

			//fprintf (stderr, "<%d> ", fcount);
			totalfixed += fcount;
			if (bbstats) bbstats->fixtime += fixtimer.getTime();
			//return 
			branch = 0;

			//fprintf (stderr, "FB"); fflush(stderr);

			if (FIND_BRANCHING) {
				if (totalfixed > 0.2 * g.EdgeCount()) {
					branch = 0;
					//ASCENT_TRIES ++;
				} else {
					RFWTimer branchtimer(true);
					if (!RUN_PRIMAL) fatal ("Cannot find branching node without primal.\n");
					//if (bestcurprimal == newprimal) 
					branch = GetBranchingNode(data, &solution, priority);
					//fprintf (stderr, "%d ", branch);
					if (bbstats) bbstats->branchtime += branchtimer.getTime();
				}
			}


			if (verbose) {
				fprintf (stderr, "%5d %.0f %.2f %.0f %3d\t\tf%.1f\t\t%.2fms+%.2fms", i, curvalue, totalvalue/(double)(i+1), best, bestcount, totalfixed/(double)(i+1), 1000.0*totaltime / (double)(i+1), 1000.0*fixtimer.getTime());
				fprintf (stderr, " Fixed %d edges. \n\n", fcount);
			}

			//fprintf (stderr, " <<EOL>> \n");

			//fflush(stderr);
		}
		//fprintf (stderr, "\n");
		//fprintf (stderr, "#%d ", bestcurprimal);
		return (EdgeCost)best;
	}


	///----------------------------------------------------------------
	/// Mark all vertices that can reach t. Returns false if this set 
    /// contains the root or any other 'live' terminal (in this case,
    /// the set may not be complete). Otherwise returns true (and 
    /// the complete set).
	///----------------------------------------------------------------

	static bool MarkSaturatedComponent (int t, int root, DualAscentData &data) {
		bool verbose = false;

		// could replace everything by a bool array + stack (though the stack has removal; bool+queue seems better. No timestamps...)

		int v = t;
		Graph *g = data.g;
        UniverseSet *marked = data.component; //vertices that belong to the component
        SimpleQueue *queue = data.queue; //queue for BFS (could be stack/DFS; should check)
        EdgeCost *rcost = data.rcost; //reduced costs

		int *unsatlist = data.unsatlist;
		int *unsattails = data.unsattails;
		data.unsatcount = 0;
		bool KEEPLIST = (unsatlist != NULL);
		if (verbose) {
			fprintf (stderr, "Keeping list: %d.\n", KEEPLIST);
			fflush(stderr);
		}

		// component and queue start with v only

        marked->Reset();
		marked->Insert(v);

		queue->Reset();
		queue->Insert(v);

		// perform a BFS from the initial vertex following only incoming saturated arcs 
		while (!queue->IsEmpty()) {
			v = queue->Remove();

			// look at all incoming arcs 
			SPGArc *a, *end;
			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (marked->Contains(w)) continue; // arc is internal to the component
				if (data.IsHidden(w)) continue;    // w temporarily removed from the graph

				int alabel = g->GetIncomingLabel(v,w,a->label); //label of incoming arc
				if (data.delarc[alabel]) continue; //arc already deleted!
				
				if (rcost[alabel] < EDGE_COST_PRECISION) { //WARNING; THERE MAY BE PRECISION ISSUES HERE
					marked->Insert(w);
					queue->Insert(w);
					
					// if we reach the root or another root component, what we started from is not an actual root component
					if ((w==root) || (w!=t && data.rlist->Contains(w))) {return false;}
				} else if (KEEPLIST) {
					// remember this is a potential incoming arc (to be processed later)
					// (we may later discover a path to w, but that's fine)
					if (verbose) {
						fprintf (stderr, "A%d %d: %d %d\n", data.unsatcount, alabel, w, v);
						fflush(stderr);
					}
					unsattails[data.unsatcount] = w;
					unsatlist[data.unsatcount++] = alabel;
				}
			}
		}
		//fprintf (stderr, "%d ", data.unsatcount);

		return true;
	}


	static void ReplayCut (Graph *g, CutRecorder &recorder, EdgeCost reftotal, int remv) {
		int m = g->EdgeCount();
		vector<EdgeCost> rcost (2*m+1);
		for (int e=1; e<=m; e++) {
			rcost[e] = rcost[e+m] = g->GetCost(e);
		}

		const bool TEST_CONDITIONAL = (remv > 0);


		const bool verbose = true;

		EdgeCost dual = 0;
		EdgeCost largecost = 999999999;
		EdgeCost biggestdiff = 0;

		int p = 0;
		int end = (int)recorder.cutlist.size();
		while (p < end) {
			//if (verbose) fprintf (stderr, "%d/%d ", p, end);
			int i = p;
			EdgeCost mincost = largecost; //WARNING: HARD CODED

			if (recorder.cutlist[p] < 0) fprintf (stderr, "Not sure what happened.\n"); 

			int count = 0;
			// find minimum in range
			while (1) {
				int alabel = recorder.cutlist[i];
				if (alabel < 0) break;
				if (rcost[alabel] < mincost) {
					if (TEST_CONDITIONAL) {
						int v, w;
						g->GetArcEndpoints(alabel, v, w);
						if (v!=remv && w!=remv) {
							mincost = rcost[alabel];
						} else {
							fprintf (stderr, "!");
						}

					} else {
						mincost = rcost[alabel];
					}
				}
				count ++;
				i++;
			}

			// some funny test regarding what would happen if a single node were removed...
			const bool TEST_SECOND = false;
			if (TEST_SECOND) {
				EdgeCost min1, min2;
				min1 = min2 = 999999999;
				for (int j=p; j<i; j++) {
					int alabel = recorder.cutlist[j];
					EdgeCost r = rcost[alabel];

					if (r < min2) {
						// WARNING: SHOULD BE TESTING IF EITHER ENDPOINT IS A NONTERMINAL 
						if (r < min1) {
							min2 = min1;
							min1 = r;
						} else {
							min2 = r;
						}
					}
				}
				EdgeCost diff = min2 - min1;
				if (biggestdiff < diff) biggestdiff = diff;
			}





			//fprintf (stderr, "Count is %d.\n", count);

			fprintf (stderr, "%d:%d ", count, i-p);

			for (int j=p; j<i; j++) {
				int alabel = recorder.cutlist[j];
				rcost[alabel] -= mincost;
			}

			if (mincost == largecost) {
				dual = mincost;
				break;
			}
			
			dual += mincost;
			p = i+1;
		}
		fprintf (stderr, "\n");

		//fprintf (stderr, "%d ", biggestdiff);


		if (verbose) {
			//if (reftotal != dual) {
				fprintf (stderr, "Dual %d -> %d : %d ", reftotal, dual, dual - reftotal);
				if (reftotal != dual) fprintf (stderr, " <<<<<<<\n");
				else fprintf (stderr, "\n");
				fflush (stderr);
			//}
		}
	}


	static void RecordCut(DualAscentData &data, CutRecorder &recorder) {
		int ucount = data.unsatcount;
		int alabel = 0;
		for (int p=0; p<ucount; p++) {
			if ((alabel=data.unsatlist[p]) < 0) continue; //arc previously identified as internal (could save this by defining rcost[-1] as infinity; this is rare, actually)
			recorder.AddArc(alabel);
		}
		recorder.CloseCut();
	}

	//-------------------------------------------------------------------------------
	// Subtract delta from the reduced costs of all arcs listed in data.unsatcount.
	// Fills "newvertices" with the tails of these arcs.
	//------------------------------------------------------------------------------
	static void CommitIncomingList (DualAscentData &data, EdgeCost delta, UniverseSet *newvertices) {
		EdgeCost *rcost = data.rcost;
		//Graph *g = data.g;
		newvertices->Reset();
		const bool DEBUG = false;
		int alabel = 0;

		int ucount = data.unsatcount;
		for (int p=0; p<ucount; p++) {
			//int alabel = data.unsatlist[p]; //this could be faster...
			if ((alabel=data.unsatlist[p]) < 0) continue; //arc previously identified as internal (could save this by defining rcost[-1] as infinity; this is rare, actually)
			if (DEBUG && delta>rcost[alabel]) fatal ("negative reduced cost");

			// decrease reduced cost; if arc is saturated, add its tail to the list of neighbors
			if ((rcost[alabel] -= delta) <= .000001) { //this is very rare
				//int v = ;				
				newvertices->Insert(data.unsattails[p]);
			}
		}
		//fprintf (stderr, "%d:%d/%d  ", delta, newvertices->Count(), data.unsatcount);

		/*
		int *start = &data.unsatlist[0];
		int *end = &data.unsatlist[data.unsatcount];
		for (int *q=start; q<end; q++) {

			//int alabel = data.unsatlist[p]; //this could be faster...
			if ((alabel=*q) < 0) continue; //arc previously identified as internal (could save this by defining rcost[-1] as infinity
			if (DEBUG && delta>rcost[alabel]) fatal ("negative reduced cost");

			// decrease reduced cost; if arc is saturated, add its tail to the list of neighbors
			if ((rcost[alabel] -= delta) <= .000001) { //this is very rare
				//int v = ;				
				int p = (int)(q - start);
				newvertices->Insert(data.unsattails[p]);
			}
		}*/
		
	}


	//-------------------------------------------------------------------------------
	// Process all arcs listed in data.unsatlist; 
	// - replaces those that are internal to "data.component" with -1 (in the list). 
	// - among the remaining arcs, find the one with minimum reduced cost. 
	// - Assign this value to delta.
	//-------------------------------------------------------------------------------

	static void EvaluateIncomingList (DualAscentData &data, EdgeCost &delta) {
		EdgeCost *rcost = data.rcost;
		delta = -1; //placeholder value
		UniverseSet *component = data.component; //vertices in the root component we are processing
		bool verbose = false;
		const bool DEBUG = false;
		//fprintf (stderr, "!");

		for (int p=0; p<data.unsatcount; p++) {
			int alabel = data.unsatlist[p]; //this could be faster...
			int v = data.unsattails[p];

			// discard arc if internal to component
			if (component->Contains(v)) {
				data.unsatlist[p] = -1;
				continue;
			}
			
			if (DEBUG) {
				Graph *g = data.g; //the graph
				int alabel = data.unsatlist[p];
				int v, w;
				g->GetArcEndpoints(alabel, v, w);
				if (verbose) {
					fprintf (stderr, "%d/%d: %d %d\n", p, data.unsatcount, v, w);
					fflush(stderr);
				}
				if (!component->Contains(w)) fatal ("all heads should be in the component");
				if (data.IsHidden(v)) fatal ("there should be no hidden vertices");
				if (data.delarc[alabel]) fatal ("there should be no deleted edge here");
				if (DEBUG && rcost[alabel] == 0) fatal ("inconsistent component");
			}


			EdgeCost c = rcost[alabel];
			if (c<delta || delta==-1) delta = c;
		}
	}

	/// Process the incoming arcs of data.component. If update is false, finds
    /// the smallest incoming reduced cost and writes it to delta. If update is
    /// false, subtracts delta from the reduced cost of every arc.

	static double EvaluateIncomingArcs (DualAscentData &data, EdgeCost &delta, bool update, UniverseSet *newvertices = NULL) {
		UniverseSet *component = data.component;
		EdgeCost *rcost = data.rcost;
		if (!update) {delta = -1;}
		Graph *g = data.g;

		if (update && newvertices) newvertices->Reset();

		int vcount = 0;
		int arccount = 0; //total number of arcs in the cut
		int satcount = 0; //number of arcs that would be saturated
		int smallcount = 0;
		double harmonic = 0;

//            foreach (int v in component.ElementEnumerator()) {
		int p, pend;
		for (component->GetBounds(p,pend); p<pend; p++) {
			int v = component->PickPos(p);
			vcount ++;
			SPGArc *a, *end;
			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head;//g.GetArcHead(pa);
				if (component->Contains(w)) continue; //internal arc
				if (data.IsHidden(w)) continue; //vertex doesn't really exist

				//int alabel = g.GetIncomingLabel(v,w,g.GetArcLabel(pa));
				int alabel = g->GetIncomingLabel(v,w,a->label); //get the label of the arc (w,v)
				if (data.delarc[alabel]) continue; //edge already deleted
				//if (data.IsArcDeleted(a)) continue;
				EdgeCost c = rcost[alabel];
				if (c==0) fatal ("inconsistent component"); //throw new Exception ("Inconsistent component.");
				arccount ++;

                if (update) {
					if (delta > c) fatal ("negative reduced cost"); //throw new Exception ("Cannot deal with negative reduced costs.");
					rcost[alabel] -= delta;
					if (newvertices) {
						if (rcost[alabel] <= 0.000001) {
							newvertices->Insert(w);
						}
					}


					//fprintf (stderr, "Should be saturating arc (%d,%d) with label %d.\n", w, v, a->label);
				} else {
					if (delta==-1 || c<delta) {
						delta = c;
						satcount = 1;
					} else if (delta == c) {
						satcount ++;
					}                         
					//these things are expensive and should not be needed
					/*
					if (c <= 50*delta) smallcount++;
					harmonic += 1.0/(double)(1 + c - delta); //not clear what this measures---changes as the algorithm goes
					*/
				}
			}
		}

		bool verbose = false;
            //double score = 1 / ((double)arccount); //1/(double)smallcount;// / (double)vcount; //(double)delta/(double)arccount; //(double)delta / ((double)arccount*(double)satcount*(double)vcount);
		//double score = 1/((double)(arccount));// * (double)arccount * (double)satcount); //(double)(arccount * smallcount);
		double score = 1/((double)(vcount));// * (double)arccount * (double)satcount); //(double)(arccount * smallcount);

		//score = 1 / (harmonic); //(double)smallcount); //1 / (double)vcount;
        if (verbose) fprintf(stderr, "V%d T%d S%d D%.3f SCORE%.3f\n", component->Count(), arccount, satcount, (double)delta, score);

		//if (!update) fprintf (stderr, "[%d] ", arccount);
		return score;
	}

	static void ComputeDistToTerminal (DualAscentData &data) {
		EdgeCost *rcost = data.rcost;
		EdgeCost *d = data.dist_to_terminal;
		EdgeCost INFINITY = INFINITE_COST; //WARNING: HARDCODED
		Graph *g = data.g;
		int n = g->VertexCount();

		static bool first = true;
		if (first) {
			first = false;
			fprintf (stderr, "Distance to terminal should be computed using implicit reduced costs from distance from root.\n");
		}


		bool verbose = false;

		if (verbose) fprintf (stderr, "\n*********************\nComputing distance to closest terminal.");
		BinaryHeap<EdgeCost> heap = BinaryHeap<EdgeCost>(n);
        int v;
        for (v=1; v<=n; v++) {
			if (g->IsTerminal(v)) {
				d[v] = 0;
				heap.Insert(v, 0);
			} else {
				d[v] = INFINITY;
			}
		}

		//do search in the reverse graph
		while (!heap.IsEmpty()) {
			EdgeCost vdist;
			unsigned int v;
			heap.RemoveFirst(v, vdist);

			//fprintf (stderr, "%d ", d[v]);

			SPGArc *a, *end;
			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (data.IsHidden(w)) continue;

				// look at incoming arc (we are duing a backward search)
				int alabel = g->GetIncomingLabel(v,a);
				if (data.delarc[alabel]) continue;

				// looking at arc wv
				EdgeCost redcost = rcost[alabel];

				const bool verbose = false;
				// difference in potential due to Dijkstra's algorithm from the source
				EdgeCost potdiff = data.dist_from_root[v] - data.dist_from_root[w];
				if (potdiff > 0) {
					redcost -= potdiff;
					if (verbose && redcost > 0) fprintf (stderr, "%d ", redcost);
					if (redcost < 0) fatal ("negative reduced cost is not possible"); //fprintf (stderr, "!");
				} else {
					if (verbose) fprintf (stderr, ".");
				}

				//EdgeCost potdiff = data.dist_from_root

				EdgeCost newdist = vdist + redcost; //rcost[alabel];


				if (newdist < d[w]) {
					d[w] = newdist;
					heap.Insert(w, newdist);
				}
			}
		}
	}



	///----------------------------------------------------------------------
	/// Fix elements by reduced cost. Returns the number of arcs eliminated.
	///----------------------------------------------------------------------
	static int FixReducedCost (DualAscentData &data, EdgeCost primal) {
		const bool USE_TO_TERMINALS = false;
		const bool COUNT_DIRECTED = false;

		Graph *g = data.g;
		int m = g->EdgeCount();
		int n = g->VertexCount();

		//fprintf (stderr, "A"); fflush(stderr);
		ComputeDistFromRoot(data.root, data);
		//fprintf (stderr, "B"); fflush(stderr);
		if (USE_TO_TERMINALS) ComputeDistToTerminal(data);
		EdgeCost *rcost = data.rcost;          //reduced costs
		EdgeCost *droot = data.dist_from_root; //distance from the root
		EdgeCost *dterm = data.dist_to_terminal; //distance to closest terminal WITH UPDATED REDUCED COSTS

		int nedges = 0;
		int nelim = 0;

		EdgeCost gap = primal - data.dualcost;

		int semicount = 0;
		for (int v=1; v<=n; v++) {
			SPGArc *a, *end;
			if (data.IsHidden(v)) continue;
			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g->GetArcHead(pa);
				if (w < v) continue;
				if (data.IsHidden(w)) continue;

				nedges ++;

				//reduced cost of arc (v,w)
                EdgeCost rout = rcost[g->GetOutgoingLabel(v,a)]; //reduced cost of (v,w)
				EdgeCost rin = rcost[g->GetIncomingLabel(v,a)];  //reduced cost of (w,v)
				bool eliminate;
				bool USE_EXTENDED = true;

				if (USE_EXTENDED) {
					//add reduced costs of paths from the root
					if (g->IsTerminal(v) && droot[v]>1.0 - EDGE_COST_PRECISION) { // WARNING: THIS IS HACKY. WE ARE ASSUMING INTEGRAL WEIGHTS WITH PERTURBATION
						fprintf (stderr, "Vertex %d is has distance %.6f from root %d.\n", v, droot[v], data.root);
						fatal ("all terminals should be reachable from root"); 
					}
                        
					EdgeCost xout = droot[v] + rout; //arc(v,w)
					EdgeCost xin = droot[w] + rin;   //arc(w,v)

					if (USE_TO_TERMINALS) {
						if (dterm[w]+dterm[v] > EDGE_COST_PRECISION) fprintf (stderr, "%d:%d ", dterm[w], dterm[v]);
						xout += dterm[w];
						xin += dterm[v];
					}
					eliminate = (xout >= gap - EDGE_COST_PRECISION) && (xin >= gap - EDGE_COST_PRECISION);
                    
					if (COUNT_DIRECTED && !eliminate && ((xout>=gap - EDGE_COST_PRECISION) || (xin>=gap - EDGE_COST_PRECISION))) {semicount ++;}

					/*
					if (USE_TO_TERMINALS && !eliminate) {
						xout = rout + data.dist_to_terminal[w];
                        xin = rin + data.dist_to_terminal[v];
						if (data.dist_to_terminal[v]>0 || data.dist_to_terminal[w]>0) 
						fprintf (stderr, "%d:%d ", data.dist_to_terminal[v], data.dist_to_terminal[w]);
                        //Console.Error.Write("{0}:{1} ", data.dist_to_terminal[w], data.dist_to_terminal[v]);
                        eliminate = (xout>=gap) && (xin>=gap);
						if (eliminate) fprintf (stderr, "p");
                        //if (eliminate) throw new Exception("Finally!");//Console.Error.Write(".");
					}*/
				} else {
					eliminate = (rout >= gap) && (rin >= gap);
				}

				if (eliminate) {
					int alabel = a->label; //g->GetArcLabel(pa);
					if (data.delarc[alabel]==false) { //WARNING: WE ARE ONLY DELETING EDGES, THE SECOND ARC IS STILL THERE. THIS IS NOT CONSISTENT.
						data.delarc[alabel] = true;
						data.delarc[alabel+m] = true;
						data.deledgecount ++;
						nelim++;
					}

					// SHOULD ADD DISTANCE TO TERMINALS HERE
				}
			}
		}
		
		if (COUNT_DIRECTED) fprintf (stderr, "%d:%d:%d ", semicount, nelim, m-nelim);

		return nelim;
	}



    /// Compute the reduced cost distance from the root to each vertex.
    /// The results are stored in data.dist_from_root.
	static void ComputeDistance (int source, vector<EdgeCost> &dist, DualAscentData &data) {
		//EdgeCost *rcost = data.rcost;
		//EdgeCost *d = data.dist_from_root;
		Graph *g = data.g;

		bool verbose = false;

		//if (verbose) fprintf(stderr, "Computing distance from root %d.", dist);
		EdgeCost INFINITY = INFINITE_COST; //WARNING: BAD!
		int n = g->VertexCount();
		unsigned int v;
		for (v=1; v<=(unsigned int)n; v++) {dist[v] = INFINITY;}

		BinaryHeap<EdgeCost> heap(n);
		heap.Insert(source, 0);
		dist[source] = 0;

		if (data.IsHidden(source)) fatal ("root cannot be hidden");

		while (!heap.IsEmpty()) {
			EdgeCost vdist;
			heap.RemoveFirst(v, vdist);
			if (verbose) fprintf (stderr, "%d:%.0f:%d ", v, (double)dist[v], (int)g->IsTerminal(v));
			SPGArc *a, *end;

			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g.GetArcHead(pa);
				if (data.IsHidden(w)) continue;
				int alabel = g->GetOutgoingLabel(v,a);
				if (data.delarc[alabel]) continue; // IS THIS CORRECT? AM I CHECKING EDGES OR VERTICES?

				//ArcCost newdist = vdist + rcost[g.GetOutgoingLabel(v, pa)];
				EdgeCost newdist = vdist + a->cost; //rcost[alabel];
				if (newdist < dist[w]) {
					dist[w] = newdist;
                    heap.Insert(w, newdist);
				}
			}
		}
	}


    /// Compute the reduced cost distance from the root to each vertex.
    /// The results are stored in data.dist_from_root.
	static void ComputeDistFromRoot (int root, DualAscentData &data) {
		EdgeCost *rcost = data.rcost;
		EdgeCost *d = data.dist_from_root;
		Graph *g = data.g;

		bool verbose = false;

		if (verbose) fprintf(stderr, "Computing distance from root %d.", root);
		EdgeCost INFINITY = INFINITE_COST; //WARNING: BAD!
		int n = g->VertexCount();
		unsigned int v;
		for (v=1; v<=(unsigned int)n; v++) {d[v] = INFINITY;}

		BinaryHeap<EdgeCost> heap(n);
		heap.Insert(root, 0);
		d[root] = 0;

		if (data.IsHidden(root)) fatal ("root cannot be hidden");

		while (!heap.IsEmpty()) {
			EdgeCost vdist;
			heap.RemoveFirst(v, vdist);
			if (verbose) fprintf (stderr, "%d:%.0f:%d ", v, (double)d[v], (int)g->IsTerminal(v));
			SPGArc *a, *end;

			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g.GetArcHead(pa);
				if (data.IsHidden(w)) continue;
				int alabel = g->GetOutgoingLabel(v,a);
				if (data.delarc[alabel]) continue; // IS THIS CORRECT? AM I CHECKING EDGES OR VERTICES?

				//ArcCost newdist = vdist + rcost[g.GetOutgoingLabel(v, pa)];
				EdgeCost newdist = vdist + rcost[alabel];
				if (newdist < d[w] - EDGE_COST_PRECISION) {
					d[w] = newdist;
                    heap.Insert(w, newdist);
				}
			}
		}
	}

    /// Compute the reduced cost distance from the root to each vertex.
    /// The results are stored in data.dist_from_root.
	static void ComputeDistToVertex (int target, int root, DualAscentData &data) {
		EdgeCost *rcost = data.rcost;
		EdgeCost *d = data.dist_from_root;
		Graph *g = data.g;

		bool verbose = false;

		if (verbose) fprintf(stderr, "Computing distance to target %d.", target );
		EdgeCost INFINITY = INFINITE_COST; //WARNING: BAD!
		int n = g->VertexCount();
		unsigned int v;
		for (v=1; v<=(unsigned int)n; v++) {d[v] = INFINITY;}

		BinaryHeap<EdgeCost> heap(n);
		heap.Insert(target, 0);
		d[target] = 0;

		if (data.IsHidden(target)) fatal ("target cannot be hidden");

		while (!heap.IsEmpty()) {
			EdgeCost vdist;
			heap.RemoveFirst(v, vdist);
			if (v == root) break;
			if (verbose) fprintf (stderr, "%d:%.0f:%d ", v, (double)d[v], (int)g->IsTerminal(v));
			SPGArc *a, *end;


			// look at incoming arcs!
			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g.GetArcHead(pa);
				if (data.IsHidden(w)) continue;
				int alabel = g->GetIncomingLabel(v,a);
				if (data.delarc[alabel]) continue; // IS THIS CORRECT? AM I CHECKING EDGES OR VERTICES?

				//ArcCost newdist = vdist + rcost[g.GetOutgoingLabel(v, pa)];
				EdgeCost newdist = vdist + rcost[alabel];
				if (newdist < d[w]) {
					d[w] = newdist;
                    heap.Insert(w, newdist);
				}
			}
		}
	}

	static EdgeCost LastComponent(DualAscentData &data, int t, int root, EdgeCost dual, EdgeCost primal, bool useprimal, EdgeCost *newrcost) {
		//ComputeDistFromRoot(root,data); // warning: this may be too expensive---could stop sooner
		//ComputeDistToTerminal(data);
		ComputeDistToVertex(t,root,data); //stores in dist_from_root: terrible!

		//fprintf (stderr, ".");

		EdgeCost dr = data.dist_from_root[root];
		//fprintf (stderr, "%d ", data.dist_from_root[v]);
		//fprintf (stderr, "Distance from root to %d is %d.\n", t, dt);
		Graph *g = data.g;
		int n = g->VertexCount();
		int m = g->EdgeCount();

		if (newrcost) {
			for (int a=1; a<=2*m; a++) newrcost[a] = data.rcost[a];
		}

		// node is already reachable from the root
		if (dr==0) {
			//fprintf (stderr, "p");	
			return dual;
		}

		// dual bound will exceed primal
		if (useprimal && dual+dr>=primal) {return primal+1;} //+1 not needed

		// invariant: at this point, all other terminals are reachable from t or from r
		EdgeCost *d = data.dist_from_root;
		EdgeCost *rcost = data.rcost;

		
		/*
		for (int v=1; v<=n; v++) {
			if (g->IsTerminal(v)) {
				EdgeCost dv = d[v];
				if (dv!=0 && dv!=dt) {
					fprintf (stderr, "v=%d dv=%d dt=%d\n", v, dv, dt);
					//fatal ("invalid terminal distance"); //fprintf (stderr, "%d ", d[v]);
				}
			}
		}
		fprintf (stderr, "\n");
		return 0;
		*/
		for (int v=1; v<=n; v++) {
			SPGArc *a, *end;
			EdgeCost dv = d[v];
			if (dv > dr) {dv = dr;}

			if (data.IsHidden(v)) continue; //make sure we have this check in other places as well

			//fprintf (stderr, "d[%d]=%d ", v, dv);

			for (g->GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g.GetArcHead(pa);
				EdgeCost dw = d[w];
				if (dw > dr) continue; //dw = dr;
				if (data.IsHidden(w)) continue;
				int alabel = g->GetOutgoingLabel(v,a);
				if (data.delarc[alabel]) continue; 

				EdgeCost delta = dv - dw;
				if (delta <= 0) {
					//fprintf (stderr, "%d ", delta);	
					continue;
				}

				//fprintf (stderr, "<%d> ", 
				//fprintf (stderr, "Root %d, arc (%d,%d), delta=%d.\n", root, v, w, delta);
				if (rcost[alabel] < delta - EDGE_COST_PRECISION) {
					fatal ("invalid distances");

				}
				if (newrcost) {
					newrcost[alabel] = rcost[alabel] - delta;
				} else {
					rcost[alabel] -= delta;
				}
			}
		}
		/*
		ComputeDistFromRoot(root, data);
		for (int v=1; v<=n; v++) {
			if (g->IsTerminal(v)) {
				EdgeCost dv = d[v];
				if (dv!=0) fatal ("bad distance");
			}
		}*/
		//if (data.dist_from_root[t]!=0) fatal ("something bad happened");

		//fprintf (stderr, "Should have increased by %d when called at %d.\n", dr, dual);
		return dual + dr;
	}


	
	static EdgeCost GetDistance (UniverseSet &vset, vector<EdgeCost> &dist) {
		const bool verbose = false;
		EdgeCost mind = INFINITE_COST;
		int p, pend;
		for (vset.GetBounds(p,pend); p<pend; p++) {
			int v = vset.PickPos(p);
			if (verbose) fprintf (stderr, "%d ", dist[v]);
			if (dist[v] < mind) mind = dist[v];
		}
		if (verbose) fprintf (stderr, " %d<", mind);
		return mind;
	}

	// c1: original component
	// c2: new vertices to be added
	static double ComputePriority (UniverseSet *c1, UniverseSet *c2, Graph *g, SteinerConfig *config){
		const int DEGREE_COUNT = 2; //will try to remove all boundary arcs
		int mode = config->ROOT_COMP_MODE;
		bool AGGRESSIVE = (mode > 2);
		int onecount = 0;
		int totalcount = 0;
		int totaldegree = 0;

		if (mode == -1) return RFWRandom::getDouble();

		
		if (mode == 0) { // vertex count
			double priority = c1->Count();
			if (c2) priority += c2->Count();
			//fprintf (stderr, "<%.1f> ", priority);
			//fflush(stderr);
			return priority;
		}

		// compute sum of degrees of all vertices in c1 and c2
		for (int i=0; i<2; i++) {
			UniverseSet *c = (i==0) ? c1 : c2;
			if (!c) continue;
			int p, pend;
			for (c->GetBounds(p,pend); p<pend; p++) {
				int deg = g->GetDegree(c->PickPos(p));
				/*
				if (deg == 0) {
					int v = c->PickPos(p);
					fprintf (stderr, "Vertex %d has degree 0; terminal:%d\n", v, g->IsTerminal(v));
				}*/
				if (deg==1) onecount ++;
				totaldegree += deg;
			}
			totalcount += c->Count();
		}
		//if (totaldegree == 0) fprintf (stderr, "Sum of degrees is %d; %d vertices considered.\n", totaldegree, c1->Count()); 

		//fprintf (stderr, "<%d,%d> ", totaldegree, totalcount);
		//double actualp = totaldegree;
		//double actualp = totaldegree - (totalcount - 1); // - DEGREE_COUNT * (totalcount - (onecount+1));
		
		double actualp = totaldegree; // estimated number of outgoing arcs
		if (mode == 1) {
			return actualp; //just sum of degrees
		}

		//if (mode > 2) AGGRESSIVE = true;

		if (AGGRESSIVE) {
			actualp -= DEGREE_COUNT * (totalcount - 1); // - DEGREE_COUNT * (totalcount - (onecount+1));
		} else { //monontone
			//actualp -= DEGREE_COUNT * (totalcount - (onecount+1)); // - DEGREE_COUNT * (totalcount - (onecount+1));
			//int highdegree = totalcount - onecount;

			// idea: we know there are at least two arcs incident to every node except the root
			// (these are the arcs in the tree)
			// To make things monotone, however, we should only account for one arc if the vertex has degree one.
			// Otherwise, when we add a vertex of degree one to the component, our estimate could decrease. The original
			// estaimte
			//int safe_internal = DEGREE_COUNT * (totalcount - (onecount + 1)); //


			int safe_internal = DEGREE_COUNT * (totalcount - 1) - onecount;	
			if (onecount > 0) {
				//safe_internal = DEGREE_COUNT * (totalcount - onecount) + onecount - 1;
			}
			//if (onecount == 0) onecount = 1;
			//int safe_internal = DEGREE_COUNT * totalcount - onecount; //- 1) - onecount;
			
			if (safe_internal < 0) {
				//fprintf (stderr, "%d:%d ", totalcount, safe_internal);
				safe_internal = 0;
			}
			
			//if (safe_internal < 0) fprintf (stderr, "%d", safe_internal);
			//else fprintf (stderr, "-");

			actualp -= safe_internal; //(DEGREE_COUNT * (highdegree - 1) + onecount); 
		} 
		return actualp;
	}


	/// Fast implementation of dual ascent method
	static void FastDualAscent (DualAscentData &data, int root, EdgeCost primal, bool useprimal, RFWLocalRandom &random, SteinerConfig *config) {
		const bool verbose = false;
		Graph *g = data.g;
		int n = g->VertexCount();
		int m = g->EdgeCount();

		bool RANDOM_VERBOSE = false;
		int unsigned seed = random.GetInteger(0,1999999999);
		
		// this is for consitency
		const bool VARIABLE_ROOT = true;

		// pick a consistent seed for a node based on the original seed and the 
		// current set of terminals
		if (VARIABLE_ROOT) {
			if (RANDOM_VERBOSE) fprintf (stderr, "<%u> ", seed);
			for (int v=1; v<=n; v++) {
				if (g->IsTerminal(v)) {
					seed = (unsigned int) (seed + v*v + 1);
				} else {
					if (g->GetDegree(v)>0) {
						seed = (unsigned int) (seed - v*v + 1);
					}
				}
			}
			if (RANDOM_VERBOSE) fprintf (stderr, "[s%u] (%d,%d,%d)\n", seed, n, m, g->TerminalCount());
			random.Randomize(seed); //17); //seed);
		}

		EdgeCost *newrcost = NULL; 

		const bool USE_DEGREES = true;
		const int DEGREE_FACTOR = 2;
		const bool VERBOSE_INDIVIDUAL = false;
		const bool VERBOSE_STATS = false;
		const bool RANDOMIZE_ORDER = true;
		const bool CHECK_WASTE = false;
		const bool USE_LIST = true;
		const bool SPECIAL_LAST = true; //process the last root component more efficiently
		const bool SPECIAL_DEBUG = false;

		const bool RECORD_CUTS = USE_LIST && !SPECIAL_LAST && true;
		if (SPECIAL_DEBUG) newrcost = new EdgeCost[2*m+1];

		CutRecorder recorder;

		BinaryHeap<double> heap(n);

		if (root==0) {
			//root = PickTerminal(*g,random);
			root = Basics::PickRandomTerminal(*g,random);
		}
		//else PickRandomTerminal(*g,random);
		data.root = root;

		//fprintf (stderr, "%d ", root);

		if (RANDOM_VERBOSE) fprintf (stderr, "root%d<%d>", root, random.GetInteger(1,1000));

		UniverseSet *rlist = data.rlist; //tentative root components
		UniverseSet *component = data.component;
		UniverseSet newvertices(n); // vertices added to the component after the fact (this could be combined with 'component')

		EdgeCost total = 0;
		int cuts = 0;
		double wasted_time = 0;
		RFWTimer localtimer(false);
		RFWTimer totaltimer(false);
		if (CHECK_WASTE) totaltimer.start();

		// ALTERNATIVE: COULD HAVE STATES
		// - 0: "free" vertex (may be added to a component)
		// - 1: vertex in the current component
		// - 2: vertex neighboring the current component
		// - 3: root vertex of some tentative component

		rlist->Reset(); // list of currently active roots (this does not need to be a universe set)
		component->Reset(); // list of vertices in the current component

		bool special_called = false;

		int vcounttotal = 0;
		int vcountsingle = 0;
		int vcountwasted = 0;

		bool USE_DISTANCES = false;

		vector<EdgeCost> rootdist;
		
		if (USE_DISTANCES) {
			rootdist.resize(n+1,-1);
			ComputeDistance(root, rootdist, data);
		}

		//fprintf (stderr, "Starting stuff.\n");

		// create one entry for each root component
		for (int v=1; v<=n; v++) {
			//priority[v] = -1;
			if (g->IsTerminal(v) && v!=root) {
				//priority[v] = 1;
				if (verbose) fprintf (stderr, "<i%d", v);
				double perturbation = RANDOMIZE_ORDER ? random.GetDouble()/200 : 0.0;

				double priority = (USE_DISTANCES) ? 1/rootdist[v] : 1;

				if (USE_DEGREES) {
					//priority = g->GetDegree(v);
					component->Insert(v);
					priority = ComputePriority(component, NULL, g, config);
					component->Remove(v);
				}


				heap.Insert(v, priority + perturbation); //randomize initial order
				if (verbose) fprintf (stderr, "> ");
				rlist->Insert(v);
			}
		}

		//fprintf (stderr, "\n\n");

		while (!heap.IsEmpty()) {
			//fprintf (stderr, "!"); fflush(stderr);
			unsigned int v;
			double p;
			heap.RemoveFirst(v,p); //pick component with best tentative score
			if (verbose) fprintf (stderr, "<%d,%.0f> ", v, p);
			if (!rlist->Contains(v)) continue; // vertex is no longer active
			rlist->Remove(v); //will add it later

			//fprintf (stderr, "%d ", v);

			// special case: this is the last root component
			if (!special_called && SPECIAL_LAST && heap.IsEmpty()) {
				//fprintf (stderr, "%d", random.GetInteger(0,1));
				special_called = true;
				//EdgeCost before = total;
				total = LastComponent(data, v, root, total, primal, useprimal, NULL); //newrcost);
				//fprintf (stderr, "%.0f ", total - before);
				break;
			}


			if (CHECK_WASTE) localtimer.start();


			// mark (in data.component) all vertices that can reach v
			bool ok = MarkSaturatedComponent(v, root, data);

			if (VERBOSE_STATS) {
				vcounttotal += component->Count();
				if (!ok) {
					vcountwasted += component->Count();
				} else if (heap.IsEmpty()) {
					vcountsingle += component->Count();
				}
			}

			// was here

			if (VERBOSE_INDIVIDUAL) fprintf (stderr, " %d:%d ", v, component->Count());

			if (!ok) {	// not really a root component: nothing else to do
				if (CHECK_WASTE) wasted_time += localtimer.getTime(); 
				continue;
			}



			// we do have a root component
			double actualp = component->Count(); // list of actual vertices
			if (USE_DISTANCES) {
				actualp = actualp / GetDistance(*component, rootdist); //rootdist[v];
			}
			if (USE_DEGREES) {
				actualp = ComputePriority (component, NULL, g, config);
			}


			if (verbose) fprintf (stderr, " priorities = [%.0f:%.0f] ", p, actualp);
			if (actualp < p - 0.1) {
				if (config->ROOT_COMP_MODE >= 0) {
					fprintf (stderr, "m"); //fflush(stderr);
				//fatal ("invalid priority");
				//fprintf (stderr, "[%.0f:%.0f]", actualp, p);
					int onecount = 0;
					int zerocount = 0;
					int pp, pend;
					for (component->GetBounds(pp,pend); pp<pend; pp++) {
						int v = component->PickPos(pp);
						if (g->GetDegree(v)<=1) {
							if (g->GetDegree(v) == 1) onecount++;
							else zerocount ++;
						}
					}
					if (onecount == 0) {
						fprintf (stderr, "[p was %.3f, actualp is %.3f (zerocount:%d)] ", p, actualp, zerocount);	
						//fatal ("priority function should be monotone");
					}
				}
			} else {
				//fprintf (stderr, "_");
			}

			// should be augmenting priority of one node
			double threshold = 1.25; //25; //25; //1.25; //1.25;
			double additive_threshold = 0;

			// if the actual number of vertices is much bigger than expected,
			// reinsert it into the heap with a different value

			// the current value is good enough if the true value is not much higher than our original estimate *or* 
			// it is still better than whatever is left in the heap
			bool good_enough = heap.IsEmpty() || (actualp <= max (threshold*p+0.1+additive_threshold, heap.PeekTopValue()));

			if (!good_enough) {
				if (verbose) fprintf (stderr, "\n");
				rlist->Insert(v);
				double perturbation = (RANDOMIZE_ORDER ? random.GetDouble()/200 : 0.0);
				
				
				heap.Insert(v,actualp + perturbation);
				if (CHECK_WASTE) wasted_time += localtimer.getTime();
				continue;
			}

			// will perform an augmentation

			//fprintf (stderr, "<<<<<     ");

			// compute bottleneck incoming arcs to this component (cut capacity)
			EdgeCost delta = 0;
			if (USE_LIST) EvaluateIncomingList(data, delta);
			else EvaluateIncomingArcs(data, delta, false); //this could be cheaper

			//fprintf (stderr, "Delta is %d.\n", delta);

			if (delta <= -EDGE_COST_PRECISION) {
				// there's no incoming arc to saturate: this component is not connected to the root
				// disconnected graphs have infinite dual bounds; there are faster ways to detect this, but this
				// is too rare to matter.
				total = INFINITE_COST; //WARNING: HARDCODED. FIX THIS.
				break;
			}
			if (verbose) fprintf (stderr, "Should be augmenting %d's component by %.8f.\n", v, (double)delta);

			// actually perform the augmentation			
			if (USE_LIST) CommitIncomingList(data, delta, &newvertices); //actually augment
			else EvaluateIncomingArcs(data, delta, true, &newvertices); //actually augment

			if (RECORD_CUTS) RecordCut(data, recorder);


			total += delta;
			if (verbose) {fprintf (stderr, "%.3f ", total); fflush(stderr);}

			bool reinsert = true;

			// this is meant to identify early on that a component will be useless
			const bool CHECK_NEW_TERMINALS= false;
			if (CHECK_NEW_TERMINALS) {
				int p, pend;
				for (newvertices.GetBounds(p,pend); p<pend; p++) {
					int w = newvertices.PickPos(p);
					if (w==v) fatal ("something fishy");
					if (w==root || rlist->Contains(w)) {
						//fprintf (stderr, "!");
						reinsert = false;
						break;
					}
				}
			}

			if (verbose) {fprintf (stderr, "R"); fflush(stderr);}

			// add this component back to the heap
			if (reinsert) {
				rlist->Insert(v);
				double newsize = component->Count() + newvertices.Count();
				if (USE_DISTANCES) {
					EdgeCost mind = min(GetDistance(*component,rootdist), GetDistance(newvertices,rootdist));
					//fprintf (stderr, "%d ", mind);
					newsize = newsize / (double)mind;
				}

				double perturbation = RANDOMIZE_ORDER ? random.GetDouble()/200 : 0;
				if (USE_DEGREES) {
					/*
					newsize = 0;
					int p, pend;
					for (component->GetBounds(p,pend); p<pend; p++) {
						newsize += g->GetDegree(component->PickPos(p));
					}
					for (newvertices.GetBounds(p,pend); p<pend; p++) {
						newsize += g->GetDegree(newvertices.PickPos(p));
					}
					newsize -= DEGREE_FACTOR * (component->Count() + newvertices.Count() - 1);
					*/
					newsize = ComputePriority(component, &newvertices, g, config);
				}



				heap.Insert(v,newsize + perturbation);
			}
                    
			//if the dual is at least as good as the primal, we may stop now...
			//(we could maybe stop at primal-1)
			if (useprimal && total>=primal - EDGE_COST_PRECISION) {break;}

			if (verbose) fprintf(stderr, "+%.3f %.3f", (double)delta, (double)total);
		}
		if (VERBOSE_STATS) {
			fprintf (stderr, "DUAL: %.3F <<<s%2.1f w%2.0f f%2.1f>>>\t\t", total, 100.0 * vcountsingle / (double)(vcounttotal), 100.0 * vcountwasted / (double)(vcounttotal), (double)vcounttotal / data.g->VertexCount());
			fflush(stderr);
		}
		if (VERBOSE_INDIVIDUAL) fprintf (stderr, "\n");

		if (verbose) {fprintf (stderr, "DONE WITH MAIN LOOP.\n"); fflush(stderr);}

		//if (verbose) fprintf (stderr, "%d ");
		if (CHECK_WASTE) fprintf (stderr, " w%.1f ", 100.0 * wasted_time / totaltimer.getTime()); 

		if (RECORD_CUTS) {
			//fprintf (stderr, "Recorded cut with %d entries.\n", recorder.cutlist.size());
			ReplayCut (data.g, recorder, total, -1);
			//int bv = GetBranchingNode(data, NULL);
			//ReplayCut (data.g, recorder, total, bv);
		}

		data.dualcost = total;
		//fprintf (stderr, "\n");
		//fprintf (stderr, "<%d> ", data.dualcost);

		if (SPECIAL_DEBUG) {
			if (total < primal) {
				fprintf (stderr, "!");
				int badcount = 0;
				for (int a=1; a<=2*m; a++) {
					if (data.rcost[a]!=newrcost[a]) {
						badcount ++;
						int v, w;
						g->GetArcEndpoints(a,v,w);
						fprintf (stderr, "Arc (%d,%d) has newred=%d and red=%d\n",  v, w, newrcost[a], data.rcost[a]);
					}
				} 
				if (badcount>0) {
					fprintf (stderr, "Primal:%d dual:%d\n", primal, total);	
					fatal("bad stuff");
				}
			} else {
				fprintf (stderr, "#");
			}
		}

		if (newrcost) delete [] newrcost;
	}


	/// Wong's dual ascent algorithm, heuristic version
    static void ScoreDualAscent (DualAscentData &data, int root, EdgeCost primal, bool useprimal, RFWLocalRandom &random) {
		bool verbose = false;
		Graph *g = data.g;

		//RFWLocalRandom random;

		if (verbose) fprintf (stderr, "SD");
		int n = g->VertexCount();
		int m = g->EdgeCount();
		int v;

		if (root==0) {
			root = Basics::PickRandomTerminal(*g, random);
			//fprintf(stderr, "<R%04d>", root);
		}
		data.root = root;

		if (verbose) fprintf(stderr, "Running dual ascent using %d as root.\n", root);

		//list of valid root components
        UniverseSet *rlist = data.rlist;
		SimpleQueue *queue = data.queue;
		UniverseSet *component = data.component;
		//UniverseSet *candlist = new UniverseSet(n);

		rlist->Reset();
		queue->Reset();
		component->Reset();
		//candlist->Reset();

		for (v=1; v<=n; v++) {
			if (g->IsTerminal(v)) {
				if (v==root) continue;
				rlist->Insert(v);
				//candlist->Insert(v);
			}
		}

		if (verbose) fprintf(stderr, "There are %d candidates to start with.\n", rlist->Count());

		EdgeCost total = 0;
		int cuts = 0;
		int recuts = 0;

		RFWTimer localtimer(true);

		//while (!candlist->IsEmpty()) {
		while (!rlist->IsEmpty()) {
			if (verbose) fprintf (stderr, "\n\n***************************\n");
			int t=0;
			bool good = false;

			int bestcount = 0;
			double bestscore = -1;
			int best = 0;

            //look at only a subset of the candidates
			//int candleft = candlist->Count();
			int candleft = rlist->Count();
			int triesleft = (int)ceil(log((double)candleft)/log((double)2));//(int)Math.Ceiling(Math.Log(candleft, 2));
			//triesleft = 20;
			//triesleft = 5;
			//triesleft /= 2;
			//triesleft = 1;
			if (triesleft < 1) triesleft = 1;
			if (triesleft > candleft) triesleft = candleft;

			//if (rlist->Count() != candlist->Count()) fatal ("violated invariant");

			if (verbose) fprintf (stderr, "Will evaluate %d of %d candidates.\n", triesleft, candleft); 

			bool failed = false;
			for (v=1; v<=n; v++) { //this could be slow!
				//if (candlist->Contains(v)) { 
				if (rlist->Contains(v)) {
					candleft --;
					if (random.GetInteger(1,candleft+1) > triesleft) {continue;}

					triesleft --;
					bool ok = MarkSaturatedComponent(v, root, data);
					//if (component->Count()>100) fprintf (stderr, "[%d %.3f]   ", data.component->Count(), (double)localtimer.getTime());
					if (verbose) fprintf (stderr, "Component rooted at %d has %d elements (ok=%d).\n", v, data.component->Count(), (int)ok);
					if (ok) {
						EdgeCost localdelta = 0;
						double score = EvaluateIncomingArcs(data, localdelta, false);  //data, ref localdelta, false);
                            
						score = score * random.GetDouble(); // * random.GetDouble();
						if (score > bestscore) {
							bestcount = 1;
							bestscore = score;
							best = v;
						} else if (score == bestscore) {
							bestcount ++;
							if (random.GetInteger(1,bestcount)==1) {best = v;}
						}
					} else {
						//candlist->Remove(v);
						rlist->Remove(v);
						failed = true;
						break;
					}
				}
				if (failed) break;
			}
			//fprintf (stderr, "Out of loop.\n");
			if (failed) continue;
			if (rlist->IsEmpty()) continue;
			//if (candlist->IsEmpty()) continue;

			t = best;
			v = best;
			//candlist->Remove(v);
			rlist->Remove(v);
			bool ok = MarkSaturatedComponent(v, root, data);
			if (!ok) fatal ("inconsistent saturation");
                
            if (verbose) fprintf (stderr, "<%.3f> ", bestscore);
			good = true;

			if (verbose) fprintf (stderr, "\n%d ", t);
			cuts ++;

			if (verbose) fprintf(stderr, "%d ", component->Count());
                
			if (good) { //could actually perform an augmentation
				EdgeCost delta = 0;
                    //ProcessIncomingArcs(data, ref delta, false); //find out how much we'll augment
				EvaluateIncomingArcs(data, delta, false);
				if (delta < 0) {
					fprintf (stderr, "DUALDISCONNECTED");
					total = 1000000000; //WARNING: HARDCODED. FIX THIS.
					break;
				}

				//ProcessIncomingArcs(data, ref delta, true);  //actually augment 
				if (verbose) fprintf (stderr, "Should be augmenting %d's component by %.3f.\n", best, (double)delta);
				EvaluateIncomingArcs(data, delta, true); //actually augment


				total += delta;
                //int temp = candlist->Count();
				int temp = rlist->Count();
				//candlist->Insert(t);
				rlist->Insert(t);
				//if (candlist->Count() != temp+1) fatal("WHAT?");
				if (rlist->Count() != temp+1) fatal ("what?");
                    
				//if the dual is at least as good as the primal, we may stop now...
				//(we could maybe stop at primal-1)
				if (useprimal && total>=primal) {break;}

                    if (verbose) fprintf(stderr, "+%.3f %.3f", (double)delta, (double)total);
			} else { //reached by another vertex
				if (verbose) fprintf(stderr, "*** ");
				continue;
			} 
		}

		fprintf (stderr, "%d ");

		bool CUTVERBOSE = false;
		if (CUTVERBOSE) fprintf(stderr, "%d:%d ", cuts, recuts);
		//delete candlist;
		data.dualcost = total;
	}



	static void ScatterBranch (vector<double> &priority, Graph &g, int setcount, int setsize, int seed, int root, EdgeCost primal, int fixto, SteinerConfig *config, ExecutionLog *executionLogPtr) {
		int n = g.VertexCount();
		int m = g.EdgeCount();
		vector<int> perm; //permutation
		vector<double> totalbound (n+1, 0);
		vector<double> count (n+1, 0);
		vector<int> vfix(n+1, -1);
		vector<int> efix(m+1, -1);

		fprintf (stderr, "Should be running scatter branch with %d sets of %d nodes.\n", setcount, setsize);
		fflush (stderr);
		// create list of possible branching vertices
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) continue;
			if (g.GetDegree(v) == 0) continue; //make sure this is really something we want to test
			perm.push_back(v);
		}
		int psize = (int) perm.size();
		if (setsize > psize) {
			fprintf (stderr, "Updating set size to %d\n", psize);
			setsize = psize;
		}

		double max_bound = 1.1 * primal;
		//fprintf (stderr, "Maximum allowed bound is %.1f\n", max_bound);

		bool USE_MINIMUM = false;

		for (int s=0; s<setcount; s++) {
			//fprintf (stderr, "%d[%d] ", s, psize); fflush(stderr);
			for (int i=0; i<setsize; i++) {
				int j = RFWRandom::getInteger(i,psize-1);
				int v = perm[j];
				perm[j] = perm[i];
				perm[i] = v;
				vfix[v] = fixto; //RFWRandom::getInteger(0,1);
				//fprintf (stderr, " %d", perm[i]);
			}
			//fprintf (stderr, "\n");

			//fprintf (stderr, "!"); fflush(stderr);
			Graph sub;
			GraphMapper mapper;
			GetSubgraph(g,sub,vfix,efix,mapper);
			// run dual ascent, run primal heuristic, and find branching node
			int bnode = 0;
			//fprintf (stderr, "Should be running dual ascent on graph with %d edges.\n", sub.EdgeCount());
			EdgeCost tempprimal = 5*primal;
			double dual = DualAscent(sub, tempprimal, efix, &mapper, bnode, 0, seed, root, NULL, NULL, false, false, config, executionLogPtr);

			//fprintf (stderr, "D%d ", dual);

			double curvalue = dual; // * (double)dual; //sqrt((double)dual);
			if (curvalue > max_bound) curvalue = max_bound;
			//fprintf (stderr, "!"); fflush(stderr);



			for (int i=0; i<setsize; i++) {
				int v = perm[i];
				count[v] ++;
				if (USE_MINIMUM) {
					if (count[v]==1 || totalbound[v]>curvalue) totalbound[v] = curvalue;
				} else {
					totalbound[v] += curvalue;
				}
				vfix[v] = -1;
			}
		}

		if (USE_MINIMUM) {
			for (int v=1; v<=n; v++) {count[v] = 1;}
		}

		sort (&perm[0], &perm[perm.size()], [&](int x, int y) {return totalbound[x]/count[x] > totalbound[y]/count[y];});
		for (int i=0; i<perm.size(); i++) {
			int v = perm[i];
			if (i < 10) fprintf (stderr, "%d %4.0f %6.2f\n", v, count[v], totalbound[v] / count[v]);
		}


		//fprintf (stderr, "\n");
		fflush(stderr);

		//vector<double> priority (n+1,-1);

		priority.resize(n+1, -1);
		fill(priority.begin(), priority.end(), -1);

		double minpriority = -1;
		for (int v=1; v<=n; v++) {
			//priority[v] = -1;
			if (g.IsTerminal(v)) continue;
			double curpriority = totalbound[v] / count[v];
			priority[v] = curpriority;
			if (minpriority<0 || curpriority<minpriority) minpriority = curpriority;
		}

		double delta = minpriority - 1;
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) continue;
			priority[v] -= delta;
		}

		bool USE_PERMUTATION = true;
		if (USE_PERMUTATION) {
			for (int v=0; v<=n; v++) {
				priority[v] = n+1;
			}

			fprintf (stderr, "Using permutation.\n");
			for (int i=0; i<perm.size(); i++) {
				int v = perm[i];
				priority[v] = (i+1.0);
				//priority[v] *= priority[v];

			}
		}


		fprintf (stderr, "Done scatter branching.\n");
		fflush (stderr);

		/*
		for (i=0; i<m; i++) {elist[i] = i+1;}
		sort(&elist[0], &elist[m],  [&](int x, int y) {return g.GetCost(x)<g.GetCost(y);});
		*/
		//return priority;
	}


	//---------------------------------------------
	// Simple DFS-based branch-and-bound algorithm 
	//---------------------------------------------
	static int BBound(Graph &g, vector<int> &vfix, vector<int> &efix, EdgeCost &primal, int depth, int seed, int root, BBStats *bbstats, BBNodeInfo bbinfo, bool probe, SteinerConfig &config, vector<double> &priority, int *ekeep, ExecutionLog *executionLogPtr = nullptr) {
		//fprintf (stderr, " <%d", depth); fflush(stderr);
		static double maxsubtree = 0; //size of maximum-weight subtree at depth VERBOSE_DEPTH


		//vector<double> priority;
		if (depth > config.DEPTH_LIMIT) {
			int m = g.EdgeCount();
			if (bbstats->ginfo) bbstats->ginfo->bbpruned = 1;
			if (ekeep) {
				int kept = 0;
				for (int e=1; e<=m; e++) {
					if (efix[e] != 0) ekeep[e] = 1;
					if (ekeep[e]) kept ++;
				}
				//fprintf (stderr, "%d/%d ", kept, m);
				if (kept > 9*m/10) config.DEPTH_LIMIT = 0;
			}
			return 0;
		}

		if (bbstats->ginfo) {
			if (bbstats->ginfo->IsSolved()) {
				bbstats->ginfo->bbpruned = 1;
				return 0;
			}
			EdgeCost fixed = bbstats->ginfo->fixed;
			EdgeCost gbound = bbstats->ginfo->UpdateBestFound(primal + fixed) - fixed;
			if (primal > gbound) {
				fprintf (stderr, "Branch-and-bound borrowing primal to %.3f.\n", gbound);	
				primal = gbound;
			}

		}


		Graph sub;
		GraphMapper mapper;
		bool verbose_depth = (depth==config.VERBOSE_DEPTH);
		if (verbose_depth) {
			fprintf (stderr, "B%d ", depth);
			fflush(stderr);
		}

		int depthpenalty[3] = {1,1,1};
		int maxdepth = g.VertexCount() + 5;

		int p = bbinfo.lastbranch;
		bbstats->nodecount[p] ++;
		if (bbstats->maxdepth < depth) bbstats->maxdepth = depth;

		int realdepth = bbinfo.splits;
		if (bbstats->maxrealdepth < realdepth) bbstats->maxrealdepth = realdepth;


		bool SAMPLE_STATS = true;
		if (SAMPLE_STATS && depth>15 && RFWRandom::getInteger(1,10000)==1) {
			bbstats->Output(stderr,true);
			fprintf (stderr, "primal %.0f\n", (double)primal);
			fflush(stderr);

		}

		//maxdepth = 10;
		
		//int depthpenalty[3] = {11,1,0};
		//int maxdepth = 60; //g.VertexCount() + 5;
		if (depth > maxdepth) {
			fprintf (stderr, "MAXDEPTH!\n");
			return 0;
		}
		//fprintf (stderr, "%d ", depth);

		if (depth > g.VertexCount()/2) fprintf (stderr, "d");


		//fprintf (stderr, "Depth %d [%d]: ", depth, bbinfo.lastbranch);
		RFWTimer subtimer(true); //time to build the appropriate subproblem
		GetSubgraph(g,sub,vfix,efix,mapper);
		if (bbstats) bbstats->subtime += subtimer.getTime();
		//fprintf (stderr, "\n\n");
		//if (depth > 45) exit(-1);
		int m = g.EdgeCount();
		vector<int> bkpfix(m+1);
		bkpfix = efix;

		// warning: this is too costly
		const bool DEBUG_FIX = false;
		if (DEBUG_FIX) {
			for (int e=1; e<=m; e++) {
				if (bkpfix[e]!=efix[e]) fatal ("didn't copy");
				//if (efix[e] != -1) fprintf (stderr, "!");
			}
		}

		int curcount = 1;


		// WARNING: this is not a true invariant due to fixing by reduced costs!
		if (depth > g.VertexCount()+1) {
			fprintf (stderr, "Reached depth %d with only %d vertices.\n", depth, g.VertexCount());
			fatal ("internal inconsistency detected");
			//fprintf (stderr, "Something weird happened.\n");
		}


		//fprintf (stderr, "When calling, mapper has %d edges.\n", mapper.oldm);
		int bnode = -1;
		priority.resize(g.VertexCount()+1,1);
		//std::fill(priority.


		int MAX_SCATTER_RUN = -1; //6; //10;
		int MAX_SCATTER_USE = MAX_SCATTER_RUN;
		bool USE_SCATTER = (MAX_SCATTER_RUN >= 0);
		//RFWRandom::randomize(17);
		//warning("z");
		if (USE_SCATTER && depth<=MAX_SCATTER_RUN) { //(depth==0 || depth==2 || depth==4 || depth==6 || depth==8 || depth==10 || depth==12)) { // || depth==6 || depth==9 || depth==12 || depth==15)) { //6 || depth==12) { // || depth == 10) { // || depth == 10) {// || depth == 8) { // || depth==7 || depth==14) {
		//if (USE_SCATTER && (depth==0 || depth==3)) {// || depth==6 || depth==9 || depth==12 || depth==15)) {
			bool last = false;// (depth==3); //(depth==15);

			int EXPZERO = 3;
			int EXPONE = 1;

			if (last) {
				EXPZERO = 2;
				EXPONE = 0;
			}

			fprintf (stderr, "Scatter branching: depth %d, exponents:%d:%d, maxscatter:%d.\n", depth, EXPZERO, EXPONE, MAX_SCATTER_USE);
			fflush(stderr);

			bool USE_AVERAGE= true;
			if (USE_AVERAGE) {
				
				/*
				if (depth == 3) {
					EXPZERO = 1;
					EXPONE = 1;
				}*/

				int n = g.VertexCount();

				if (EXPZERO == 0) {
					priority.resize(n+1);
					std::fill(priority.begin(), priority.end(), 1);
				} else {
					ScatterBranch(priority, sub, 1000, 32 - depth, seed, root, primal, 0, &config, executionLogPtr);
					if (EXPZERO > 1) {
						for (int v=0; v<=n; v++) {
							double p = priority[v];
							for (int i=1; i<EXPZERO; i++) priority[v] *= p;
							//priority[v] = priority[v] * priority[v];  //priority[v] * priority[v]; //p1[v]; // * priority[v] * priority[v] * p1[v]; // * p1[v]; // * p1[v]; // * p1[v];// * p1[v]); // + 2*std::max(priority[v], p1[v]); //priority[v]*priority[v]; // + p1[v]*p1[v]; //make p1 more important
						}
					}
				}

				if (EXPONE > 0) {
					vector<double> p1 (n+1,-1);
					ScatterBranch(p1, sub, 2000, 10, seed, root, primal, 1, &config, executionLogPtr);
					for (int v=0; v<=n; v++) {
						double p = p1[v];
						for (int i=0; i<EXPONE; i++) priority[v] *= p;
						//priority[v] = priority[v] * priority[v];  //priority[v] * priority[v]; //p1[v]; // * priority[v] * priority[v] * p1[v]; // * p1[v]; // * p1[v]; // * p1[v];// * p1[v]); // + 2*std::max(priority[v], p1[v]); //priority[v]*priority[v]; // + p1[v]*p1[v]; //make p1 more important
					}
				}
			} else {
				if (depth>0) {
					ScatterBranch(priority, sub, 1000, 32 - depth, seed, root, primal, 0, &config, executionLogPtr);
				} else {
					ScatterBranch(priority, sub, 1000, 10, seed, root, primal, 1, &config, executionLogPtr);
				}
			}
		}

		// run dual ascent, run primal heuristic, and find branching node
		EdgeCost dual = DualAscent(sub, primal, efix, &mapper, bnode, 0, seed, root, bbstats, depth>MAX_SCATTER_USE ? NULL : &priority[0], true, true, &config, executionLogPtr);

		mapper.Destroy();

		if (depth == 0) {
			printf ("rootprimal %.6f\n", primal);
			printf ("rootdual %.6f\n", dual);
			EdgeCost gap = primal - dual;
			bool solved = false;
			if (gap <= 0.000001) {
				gap = 0;
				solved = true;
			}
			double pctgap = 100.0 * (gap / dual);
			printf ("rootabsgap %.6f\n", gap);
			printf ("rootpctgap %.6f\n", pctgap);
			printf ("rootsolved %d\n", (int)solved);
			fprintf (stderr, "[Root] primal:%.3f dual:%.3f gap:%.3f%%\n", primal, dual, pctgap);
		}



		if (dual >= primal - 0.000001) { // node is solved!
			//warning("b");
			if (verbose_depth) fprintf (stderr, "P%d ", dual);
			bbstats->MarkAsSolved(bbinfo);
			//fprintf (stderr, "%.2f ", 100 * bbstats->fraction_solved);
			//fprintf (stderr, "depth:%d:%d splits:%d lastbranch:%d laststrict:%d\n", depth, bbinfo.depth, bbinfo.splits, bbinfo.lastbranch, bbinfo.laststrict);
		} else {
			//warning("c");
			int n = sub.VertexCount();

			if (verbose_depth) fflush(stderr);
			if (bnode < 0) { // could not find a branching node; this should not have happened!
				//warning("d");
				fprintf (stderr, "No branching vertex available!\n"); //WARNING: SHOULD BE MORE CAREFUL WITH THIS... SHOULD BE COMPUTING MST IF WE ONLY HAVE TERMINALS
				bnode = 0;
				double totalvalue = 0;
				for (int v=0; v<=n; v++) {
					if (sub.IsTerminal(v)) continue;
					if (vfix[v] != -1) continue;
					double curvalue = sub.GetDegree(v)*sub.GetDegree(v)*sub.GetDegree(v);
					totalvalue += curvalue;
					if (bnode==0 || RFWRandom::getDouble()*totalvalue<=curvalue) bnode = v; //sub.GetDegree(v)>sub.GetDegree(bnode)) bnode = v;
				}
				if (bnode == 0) fatal ("ran out of vertices");
			} else if (bnode==0) { // algorithm decide we should not branch because we fixed lots of stuff
				//warning ("e");
				if (verbose_depth) fprintf (stderr, "F ");

				// recursively branch
				curcount += BBound (g, vfix, efix, primal, depth+depthpenalty[2], seed, root, bbstats, BBNodeInfo(bbinfo,2), probe, config, priority, ekeep);


				//fprintf (stderr, "Should be branching on %d.\n", bnode);
			} else { //we should branch
				if (sub.GetDegree(bnode)<3) fprintf (stderr, " %d", sub.GetDegree(bnode));

				EdgeCost gap = primal - dual;
				bool ALLOW_SB = config.STRONG_BRANCHING;
				int smalldepth = config.SB_LOW;
				int largedepth = config.SB_HIGH;
				int refgap = config.SB_GAP;


				bool STRONG_BRANCHING = (depth<smalldepth) || (depth<largedepth && gap>refgap); //(depth<20); //(gap > 400);

				STRONG_BRANCHING = STRONG_BRANCHING && ALLOW_SB;
				if (STRONG_BRANCHING) {
					//warning("g");
					fprintf (stderr, "Gap is %d => strong branch at depth %d (side %d)!\n", gap, depth, bbinfo.lastbranch);
					int newbnode = StrongBranching(sub,primal,efix,root,seed,depth,&config, executionLogPtr);
					fprintf (stderr, "Branching on %d (had %d before).\n", newbnode, bnode);
					bnode = newbnode;
				}


				bool STATIC_PRIORITIES = false; //(depth <= 32);
				if (STATIC_PRIORITIES) {
					double maxpriority = -1;
					int bestv = -1;
					for (int v=1; v<=n; v++) {
						//double curpriority = (double)sub.GetDegree(v) / priority[v]; // * sqrt(sub.GetDegree(v)); 
						double curpriority = 1.0 / priority[v]; // * sqrt(sub.GetDegree(v)); 
						//double curpriority = priority[v]; // * sqrt(sub.GetDegree(v)); 
						if (curpriority < maxpriority) continue;
						if (g.IsTerminal(v)) continue;
						if (vfix[v] != -1) continue;
						bestv = v;
						maxpriority = curpriority;
					}
					if (bestv >= 0) {
						bnode = bestv;
						//fprintf (stderr, "<%d:%.2f>", depth, maxpriority);
						//fprintf (stderr, "%d ", bnode, depth);
					}
				}


				//warning("h");
				//if (depth == 10) fprintf (stdout, "%d %d\n", bnode, g.GetDegree(bnode));
				bool zerofirst;
				int fb = config.FIRST_BRANCH;
				if (fb <= 0) {zerofirst = true;} 
				else if (fb == 1) {zerofirst = false;} 
				else {zerofirst = (RFWRandom::getInteger(1,fb)==1);}
				/*
				switch (config.FIRST_BRANCH) {
					case 0: zerofirst = true; break;
					case 1: zerofirst = false

				bool zerofirst = (config.ZERO_FIRST != 0); //(RFWRandom::getInteger(0,1)>0);
				//fprintf (stderr, "[[%d]] ", zerofirst);
				zerofirst = (RFWRandom::getInteger(0,1)>0);
				*/
				// this is a branching node, so we remember what we did last
				int curlast = bbinfo.laststrict;
				bbstats->branchcount[curlast]++;

				int maxchildren = 2;
				if (probe) {
					zerofirst = (RFWRandom::getInteger(0,1)>0);
					maxchildren = 1;
				}
				//fprintf (stderr, "<%d> ", probe);

				double before = bbstats->GetTotalNodes();

				//warning("i");
				for (int i=0; i<maxchildren; i++) {
					vfix[bnode] = zerofirst ? i : 1-i;
					//fprintf (stderr, " fixing %d to %d.", bnode, vfix[bnode]);
					int newdepth = depth+depthpenalty[i]; //depth+(2-i); //add 2 for 0, 1 for 1.
					curcount += BBound (g, vfix, efix, primal, newdepth, seed, root, bbstats, BBNodeInfo(bbinfo,vfix[bnode]), probe, config, priority, ekeep);
				}


				//vfix[bnode] = 1;
				//curcount += BBound (g, vfix, efix, primal, depth+1);

				vfix[bnode] = -1;
				if (verbose_depth) {
					bool OUTPUT_SUBPROBLEMS = false;
					double after = bbstats->GetTotalNodes();
					double difference = after - before;
					if (difference > maxsubtree) {
						maxsubtree = difference;
						fprintf (stderr, "*");
						if (OUTPUT_SUBPROBLEMS) OutputSubproblem(g, vfix, efix, maxsubtree, primal);
					}
					fprintf (stderr, "<%.0f> ", difference);
				}
			}
			efix = bkpfix; 
		}
		return curcount;
	}

	static void CountSaturated (int v, DualAscentData &data, int &insatcount, int &outsatcount, int &degree, EdgeCost &mincost) {
		insatcount = outsatcount = degree = 0;
		mincost = 999999999;

		int tneighbors = 0;

		Graph *g = data.g;
		SPGArc *a, *end;
		for (g->GetBounds(v,a,end); a<end; a++) {
			if (data.delarc[a->label]) continue;
			if (data.IsHidden(a->head)) continue;
			degree ++;
            int alabel = g->GetIncomingLabel(v, a);
			if (data.rcost[alabel]==0) insatcount ++;
            alabel = g->GetOutgoingLabel(v, a);
            if (data.rcost[alabel]==0) outsatcount ++;
			if (a->cost < mincost) mincost = a->cost;

			if (g->IsTerminal(a->head)) tneighbors ++;
		}

		//fprintf (stderr, "%d ", tneighbors);

		//mincost = tneighbors;

		//fprintf (stderr, "[%d %d %d] ", degree, insatcount, outsatcount);
	}

	//---------------------------------------------
	// deterministic routine to get branching node
	// other things we could consider:
	// - distance from root
	// - actual degree
	// - some pool of solutions
	//---------------------------------------------
	static int GetBranchingNode (DualAscentData &data, SteinerSolution *solution, double *priority = NULL) {
		Graph *g = data.g;
		int n = g->VertexCount();
		//double m = g->EdgeCount();
		const bool USE_PRIMAL = (solution != NULL);
		const bool RANDOMIZE = true;

		/*vector<EdgeCost> dist(n+1);
		ComputeDistance(data.root, dist, data);
		for (int v=1; v<=n; v++) {fprintf (stderr, "%d ", dist[v]);}
		fprintf (stderr, "\n");*/

		//priority = NULL;


		double bestscore = 0;
		int best = 0;
		int nonzero = 0;

		const bool USE_LOWER_ORDER = (!priority); //!USE_PRIMAL; 
		const double DEGREE_ADDITIVE = .1;

		//if (USE_LOWER_ORDER) fprintf (stderr, "L");

		//fprintf (stderr, "%d ", solution->EdgeCount());

		for (int v=1; v<=n; v++) {
			if (g->IsTerminal(v)) continue;
			if (g->GetDegree(v)==0) continue; // OK, this is very hacky; this is meant to represent vertices that are already fixed to zero
			if (data.IsHidden(v)) continue;  //not really needed

			double score = 0;
			if (USE_LOWER_ORDER) {
				int insatcount, outsatcount, degree;
				EdgeCost mincost;
				CountSaturated(v,data, insatcount,outsatcount, degree, mincost);

				//fprintf (stderr, "%d ", mincost);

				// original score...
				//double score = (double)insatcount * (double)outsatcount * (double)degree;

				score =  (double)(insatcount + outsatcount + degree); // / (double)mincost; // * ((double)mincost / (double)1000);/// / ((double)dist[v]); // * (double)degree;
				if (RANDOMIZE) score *= RFWRandom::getDouble();
			
				//score = 1 / (dist[v]); // / 100000);
				// use this for perturbation only
				score = (score / (double)(4*n));
			}

			if (USE_PRIMAL) {
				//score *= solution->GetDegree(v) * solution->GetDegree(v);
				int soldegree = solution->GetDegree(v);
				if (soldegree > 0) nonzero++;

				bool BINARY_SOLUTION = false;
				if (BINARY_SOLUTION && (soldegree>0)) soldegree = soldegree * soldegree;

				double scale = 2.0 * (double)n * (double)n;
				double degterm = (double)(soldegree) + 1.0 / ((double)n * (double)n);


				degterm += DEGREE_ADDITIVE;


				if (priority) {
					double priterm = priority[v]; //((double) (*priority)[v]);
					score += scale * degterm / priterm; // * solution->GetDegree(v);
					//if (soldegree > 0) fprintf (stderr, "d%d:p%f:s%f  ", solution->GetDegree(v), (*priority)[v], score );
				}
				else score += soldegree;
				
				/*
				if (score>0) {
					fprintf (stderr, "%f:%d ", score, solution->GetDegree(v));
				}
				*/

				//score += RFWRandom::getInteger(1,solution->GetDegree(v) * solution->GetDegree(v)  * solution->GetDegree(v)  * solution->GetDegree(v));
				//score += (solution->GetDegree(v) > 2) ? 2 : 1; // * solution->GetDegree(v);
			}

			if (best==0 || score > bestscore) {
				//if (best == 0) fprintf (stderr, "\n");
				bestscore = score;
				best = v;
				//fprintf (stderr, "%.3f[%d,%d] ", bestscore, solution->GetDegree(v), g->TerminalCount());
			}
		}
		//fprintf (stderr, "\n");
		//fprintf (stderr, "Best is %d.\n", best);
		
		//fprintf (stderr, "%d:%.0f:%d    ", best, bestscore, nonzero);

		return best;
	}

	static EdgeCost SPHDual (DualAscentData &data, SteinerSolution &solution, bool UNDIRECTED) {
		const bool verbose = false;
        unsigned int v;
		Graph *g = data.g;
        int n = g->VertexCount();
        int m = g->EdgeCount();
		int root = data.root;
		//bool UNDIRECTED = false;
		static bool first = true;

        if (verbose) fprintf (stderr, "Running SPHDual from %d... ", root);

        BinaryHeap<EdgeCost> heap(n); //vertices to scan
		vector <int> parc(n+1,-1); //parent arc
		vector <EdgeCost> dist(n+1,-1); //current distance

		//int nonterminals_seen = 0;

        solution.Reset();
        heap.Insert(root, 0);
        dist[root] = 0; //root already in the tree, but no parent edge
        parc[root] = 0; //already in the tree

        int nterm = g->TerminalCount();
        int tcount = 0;
        if (g->IsTerminal(root)) tcount ++;
		else fatal ("root must be a terminal");

        while (tcount < nterm) {
			if (verbose) {fprintf (stderr, "Here."); fflush(stderr);}
			if (heap.IsEmpty()) {break;}

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
            if (g->IsTerminal(v) && (parc[v]>0)) {
				int curpath = 0;
				if (verbose) fprintf (stderr, "Adding terminal %d.\n", v);
				//fflush(stderr);
				tcount ++;
                int w = v;
                EdgeCost addedcost = 0;
                do {
					int e = parc[w];
					curpath ++;
                    if (v!=w) heap.Insert(w, 0);
					//inscount ++;
                    dist[w] = 0; //its distance label is zero
                    parc[w] = 0; //w now belongs to the tree
                    if (!solution.Insert(e)) fatal ("failed to add edge");
					EdgeCost acost = g->GetCost(e);
                    addedcost += acost;
                    w = g->GetOther(e,w);
				} while (parc[w]>0);
				if (addedcost != curcost && first) {
					first = false;
					fprintf (stderr, "addedcost:%lg curcost:%lg\n", addedcost, curcost);
					fprintf (stderr, "Inconsistent path in SPH.");
				}
				//pathsum += curpath;
				curcost = 0;
			}

			if (verbose) fprintf (stderr, "Processing a non-terminal.\n");

			// scan outgoing arcs, add anything that improves to the solution
			//scancount ++;
			SPGArc *arc, *end;
			for (g->GetBounds(v,arc,end); arc<end; arc++) {
				int w = arc->head;
				if (data.IsHidden(w)) continue;
				int alabel = g->GetOutgoingLabel(v,arc);

				bool invalid = (data.rcost[alabel]>EDGE_COST_PRECISION || data.delarc[alabel]);
				//if (!data.delarc[alabel]) invalid = false;


				if (UNDIRECTED && invalid) {
					//fprintf (stderr, "u");
					int rlabel = g->GetIncomingLabel(v,arc);
					invalid = (data.rcost[rlabel]>EDGE_COST_PRECISION || data.delarc[rlabel]);
				}

				// we only care about saturated edges
				//if (data.rcost[alabel]>0 || data.delarc[alabel]) continue; 
				if (invalid) continue;
                if ((w==root) || (solution.GetDegree(w)>0)) continue; //already in solution

                EdgeCost acost = arc->cost; //g.GetCost(arc->label);                    
                EdgeCost newcost = acost + curcost;

                if (newcost < -EDGE_COST_PRECISION) {
					fprintf (stderr, "Arc %d + cost %d led to newcost %d.", acost, curcost, newcost);
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
		if (verbose) fprintf (stderr, "Found solution %.3f (%d %d)", solution.GetCost(), tcount, nterm);

		//fprintf (stderr, "r%d:%d/%d ", root, tcount, nterm);

		if (tcount == nterm) return solution.GetCost();
		else return INFINITE_COST; //DISCONNECTED GRAPH
	}

	static int StrongBranching (Graph &g, EdgeCost &primal, vector<int> &efix, int root, int seed, int depth, SteinerConfig *config, ExecutionLog *executionLogPtr) {
		int n = g.VertexCount();

		int bestb = 0;
		double bestvalue = 0;
		bool verbose = true;

		/*
		if (verbose) {
			fprintf (stderr, "Should be performing strong branching.\n");
			fflush(stderr);
		}*/

		//int seed = 17;
		//fprintf (stderr, "SEED MISSING!\n");
		double maxit = 3;
		double zeroweight = max(10,10*depth); 

		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) continue;
			if (g.GetDegree(v) == 0) continue;

			double total[2] = {0.0,0.0};
			int count[2] = {0,0};
			bool onefixed = false;
			bool zerofixed = false;
			for (int t=0; t<maxit; t++) {
				bool zeroside = (t!=0);
				int side = zeroside ? 0 : 1;
				count[side] ++;
				//DualAscentData data(&g);
				//if (t%2 == 0) 
				//if (zeroside) data.Hide(v);
				//else g.MakeTerminal(v);
				int tempb = 0;
				if (!zeroside) g.MakeTerminal(v);
				EdgeCost dual = DualAscent(g, primal, efix, NULL, tempb, zeroside ? v : -1, v*seed*(t+1), root, NULL, NULL, false, false, config, executionLogPtr);
				if (!zeroside) g.UnmakeTerminal(v);

				total[side] += dual;
				if (dual >= 999999999) {
					if (zeroside) zerofixed = true;
					else onefixed = true;
				}
				//data.UnhideVertices();
				//if (total < 0) total = dual;
				//else if (total > dual) total = dual;
				//if (dual > total) total = dual; 
				//total += dual;
			}
			//total = total * g.GetDegree(v);
			bool somefixed = zerofixed || onefixed;
			for (int i=0; i<2; i++) {
				if (count[i] > 0) total[i] /= (double)count[i];
			}

			double curvalue = (zeroweight * total[0] + total[1]) / (1.0 + zeroweight); 
			if (!somefixed && curvalue > bestvalue) {
				bestvalue = curvalue;
				bestb = v;
				double norm = bestvalue; //exp((1/(double)maxit)* log (total)); //root
				fprintf (stderr, " %d:%.1f ", v, norm);
			}
			if (verbose) fflush(stderr);
		}
		//double norm = exp((1/(double)maxit)* log (bestvalue));
		//fprintf (stderr, "%d:%.2f ", bestb, norm);
		return bestb;
	}



	static int PickBestRoot(Graph &g, int primal, int seed, int tries) {

		int n = g.VertexCount();
		vector<double> totaldepth(n+1,0);
		RFWLocalRandom random(seed);
		for (int t=0; t<tries; t++) {
			for (int v=1; v<=n; v++) {
				if (g.IsTerminal(v)) {
					totaldepth[v] += RunProbe(g,primal,random.GetInteger(1,999999999),v);
				}
			}
		}

		int best = -1;
		for (int v=1; v<=n; v++) {
			if (!g.IsTerminal(v)) continue;
			fprintf (stderr, "%d:%.1f ", v, totaldepth[v] / tries);
			if (best<0 || totaldepth[v] < totaldepth[best]) best = v;
		}
		return best;
	}

	static double RunProbe (Graph &g, EdgeCost primal, int seed, int root) {
		int n = g.VertexCount();
		int m = g.EdgeCount();
		vector<int> vfix (n+1, -1);
		vector<int> efix (m+1, -1);
		BBStats bbstats;
		RFWLocalRandom random(seed);
		SteinerConfig config;
		vector<double> priority; 
		fprintf (stderr, "WARNING: PRIORITY NOT SET!\n");
		exit(-1);
		int count = BBound(g, vfix, efix, primal, 0, seed, root, &bbstats, BBNodeInfo(), true, config, priority, NULL);
		return bbstats.GetTotalNodes();
	}


	/*-----------------------------------------------------------------
	 | Build a subgraph of the original graph taking into account all
	 | vertices and edges previously fixed. Changes 'sub'
	 *----------------------------------------------------------------*/
	static void GetSubgraph (Graph &g, Graph &sub, vector<int> &vfix, vector<int> &efix, GraphMapper &mapper) {
		int m = g.EdgeCount();
		int n = g.VertexCount();
		mapper.Reset(n,m);
		//fprintf (stderr, "Mapper has %d edges.\n", mapper.oldm);

		sub.SetVertices(n);
		sub.SetEdges(m); //upper bound...
		//sub.SetFixedCost(g.GetFixedCost());
		//fprintf (stderr, "[Vertex 82 is fixed to %d.]", vfix[82]);
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v) || vfix[v]==1) {sub.MakeTerminal(v);}
		}
		//fprintf (stderr, "\n");
		//bool has_nonterminal = false;

		//vector<int> subdegree (n+1,0); //vertex degrees in the subgraph
		for (int e=1; e<=m; e++) { //visit all edges in the original graph
			// edge only exists if edge itself has not been fixed to zero
			// and if neither endpoint has been removed
			if (efix[e] == 0) continue; //edge itself removed
			if (efix[e] == 1) fatal ("not supported"); // cannot support fixed edges yet
			int v,w;
			g.GetEndpoints(e,v,w);
			if (vfix[v]==0 || vfix[w]==0) continue; 
			int newid = sub.AddEdge(v,w,g.GetCost(e));
			//subdegree[v] ++;
			//subdegree[w] ++;
			mapper.e2new[e] = newid; //remember edge ids for output
			//fprintf (stderr, ".%d ", newid);
		}

		//fprintf (stderr, "\n");
		/*
		for (int v=1; v<=n; v++) {

			if (subdegree[v] <= 2 && subdegree[v]>0) {
				if (subdegree[v] == 1 || !sub.IsTerminal(v)) {
					fprintf (stderr, "%d:%d%c ", v, subdegree[v], sub.IsTerminal(v) ? '*' : ' ');
				}
			}
		}*/

		//fprintf (stderr, "NT%d ", has_nonterminal);
		//if (!has_nonterminal) fprintf (stderr, "INSTANCE ONLY HAS TERMINALS!\n");

		sub.Commit();
		
		const bool verbose = false;
		
		if (verbose) fprintf (stderr, " sub:%d:%d:%d ", sub.VertexCount(), sub.EdgeCount(), sub.TerminalCount());
	}




	static void OutputSubproblem (Graph &g, vector<int> &vfix, vector<int> &efix, double bbcount, EdgeCost primal) {
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int ecount = 0;

		Graph sub;
		GraphMapper mapper;
		GetSubgraph(g,sub,vfix,efix, mapper);


		fprintf (stderr, "Original graph has %d vertices and %d edges.\n", g.VertexCount(), g.EdgeCount());
		fprintf (stderr, "Should be outputting a graph with %d vertices and %d edges.\n", sub.VertexCount(), sub.EdgeCount());

		// only keep high-degree vertices
		int vcount = 0;
		int zerocount = 0;
		int acount = 0;
		vector<int> old2new(n+1, -1);
		for (int v=1; v<=n; v++) {
			if (sub.GetDegree(v)==0) {zerocount ++;}
			vcount ++;
			old2new[v] = vcount;
		}
		fprintf (stderr, "There are %d vertices with degree %d.\n", zerocount, 0);

		char filename [2048];
		sprintf (filename, "subproblem-%d-%07.0f-%d.stp", primal, bbcount, zerocount); 
		fprintf (stderr, "<%s>", filename);

		FILE *file = fopen (filename, "w");
		if (!file) {
			fprintf (stderr, "Could not open <%s> for writing.\n", filename);
		} else {
			fprintf (file, "33D32945 STP File, STP Format Version 1.0\n\n");
			fprintf (file, "SECTION Graph\n");
			fprintf (file, "Nodes %d\n", n);
			fprintf (file, "Edges %d\n", sub.EdgeCount());

			int sm = sub.EdgeCount();
			for (int e=1; e<=sm; e++) {
				int v, w;
				sub.GetEndpoints(e, v, w);
				EdgeCost cost = sub.GetCost(e);
				//fprintf (stderr, "WARNING: check double output\n");
				fprintf (file, "E %d %d %.0f\n", v, w, (double)cost);
			}
			fprintf (file, "END\n\n");

			fprintf (file, "SECTION Terminals\n");
			fprintf (file, "Terminals %d\n", sub.TerminalCount());
			for (int v=1; v<=n; v++) {
				if (sub.IsTerminal(v)) {
					fprintf (file, "T %d\n", v);
				}
			}
			fprintf (file, "END\n\n");
			fprintf (file, "EOF\n");

			fclose(file);

		}
	}

};

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

#define DEBUG_VERTEX -1
#undef INFINITY

using namespace std;

class GlobalInfo {
	EdgeCost bestfound;
	int solved;

public:
	int bbpruned; //true iff bb pruned at a node that was not yet solved
	EdgeCost fixed; //hack because bb cannot handle fixed costs
	// 
	EdgeCost UpdateBestFound(EdgeCost bf) {
		double answer = bestfound;
		if (bf != answer) {
#pragma omp critical
			{
				if (bf < bestfound) {
					bestfound = bf;
					fprintf (stderr, "[[[ %.2f ]]] ", bestfound);
				}
				answer = bestfound;
			}
		}
		return answer;
	}

	inline bool IsSolved() {
		return (solved!=0);
	}

	void MakeSolved() {
		solved = 1;
	}

	GlobalInfo() {
		bestfound = INFINITE_COST;
		fixed = 0;
		solved = 0;
		bbpruned = 0;
	}
};


struct LocalArc {
public:
	int v;
    int w;
    EdgeCost cost;

	LocalArc(int _v, int _w, EdgeCost _cost) {
		v = _v;
		w = _w;
		cost = _cost;
	}

	void Set(int _v, int _w, EdgeCost _cost) {
		v = _v;
        w = _w;
        cost = _cost;
	}
    
	LocalArc() {v=w=0; cost=0;}
};


class CutRecorder {
public:
	vector<int> cutlist;

	void Reset() {
		cutlist.clear();
	}

	void Reset(int size) {
		cutlist.reserve(size);
		cutlist.clear();
	}

	inline void AddArc(int alabel) {
		cutlist.push_back(alabel);
	}

	inline void CloseCut() {
		cutlist.push_back(-1);
	}

};

class GraphMapper {
private:
	void Init() {
		oldn = oldm = 0;
		v2new = e2new = NULL;
	}

public:
	int *v2new;
	int *e2new;
	int oldn, oldm;

	GraphMapper() {
		Init();
	}

	void Reset(int _oldn, int _oldm) {
		oldn = _oldn;
		oldm = _oldm;
		v2new = new int [oldn+1];
		e2new = new int [oldm+1];
		for (int e=1; e<=oldm; e++) e2new[e] = -1;
		for (int v=1; v<=oldn; v++) v2new[v] = -1;
	}

	~GraphMapper() {
		if (v2new) delete [] v2new;
		if (e2new) delete [] e2new;
	}

	void Destroy() {
		if (v2new) delete [] v2new;
		if (e2new) delete [] e2new;
		Init();
	}

};


class SPGSolver {
private:
	static void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
	}

	static void warning (const string &msg) {
		fprintf (stderr, "%s", msg.c_str());
		fflush(stderr);
	}

	/*
	static void ShowUsage() {
		fprintf (stderr, "Usage: steiner <filename> <upper bound> [-prep 1] [-seed <s>]\n");
	}*/

	static void ExtractFilename (char *filename, char *name) {
		strcpy (name, filename);

		//fprintf (stderr, "<%s> ", name);
		//exit(-1);

		char *lastsep = NULL;
		char *lastdot = NULL;
		
		char *end = &name[strlen(name)];
		char *p = name;

		for (p=name; p<=end; p++) {
			if (*p=='/' || *p=='\\') lastsep = p; 
			if (*p=='.') lastdot = p;
		}

		if (lastsep == NULL) {
			lastsep = name;
		} else {
			lastsep ++;
		}

		if (lastdot == NULL) {
			lastdot = end;
		} 

		// avoid weird cases
		if (lastdot <= lastsep) {
			lastdot = end;
			lastsep = name;
		}

		*lastdot = 0;

		int i = 0;
		while (1) {
			name[i] = *lastsep;
			if (name[i] == 0) break;
			i++;
			lastsep ++;
		}

		//fprintf (file, "graph %s\n", lastsep);

		//fprintf (file, "graph %s\n", name);
		//delete [] name;
	}



	static void OutputGraphStats (FILE *file, Graph &g, char *filename) {
		fprintf (file, "file %s\n", filename);
		fprintf (file, "nvertices %d\n", g.VertexCount());
		fprintf (file, "nedges %d\n", g.EdgeCount());
		fprintf (file, "nterminals %d\n", g.TerminalCount());


	}

public:

	static EdgeCost ReadBound(FILE *logfile, char *graphname) {
		EdgeCost answer = -1;
		
		FILE *file = fopen ("bounds.txt", "r");
		if (!file) {
			fprintf (stderr, "Could not find bounds file.\n");	
			fflush (stderr);
			return answer;
		}
		fprintf (stderr, "Reading bounds file... ");
		fflush (stderr);

		double solution;
		const int BUFSIZE = 65534; //1048574;
		char buffer [BUFSIZE+2];
		char bufseries [BUFSIZE+2];
		char bufclass [BUFSIZE+2];
		char bufinstance [BUFSIZE+2];

		sprintf (bufseries, "undefined");
		sprintf (bufclass, "undefined");

		while (fgets(buffer, BUFSIZE, file)!=0) {
			//fprintf (stderr, "<%s> ", buffer);
			if (sscanf (buffer, "#series %s", bufseries) > 0) {
				continue;
			}
			if (sscanf (buffer, "#class %s", bufclass) > 0) {
				continue;
			}
				

			if (sscanf (buffer, "%s %lg", bufinstance, &solution)>0) {
				if (strcmp (bufinstance, graphname) == 0) {
					answer = (EdgeCost)solution;
					fprintf (logfile, "instance %s\n", graphname);
					fprintf (logfile, "series %s\n", bufseries);
					fprintf (logfile, "class %s\n", bufclass);
					fprintf (logfile, "bestknown %.20f\n", (double)solution);
				}
			}
		}
		fclose(file);
		//fprintf (stderr, "done (solution is %.0f).\n", (double)answer);
		//fflush (stderr);
		return answer;
	}

	typedef enum { MS_PLAIN, MS_COMBINATION, MS_BINARY, MS_MULTILEVEL, MS_TIMEBOUNDEDCOMBINATION, MS_TIMEBOUNDEDMULTILEVEL, MS_TIMEBOUNDEDADAPTIVE, MS_NUMBER } MSType;
	static void ShowUsage () {
		fprintf (stderr, "Usage: steiner <stp_file> [-bb] [-ub] [-prep] [-msit] [-seed] [-mstype]\n");
		fprintf (stderr, "Valid types: plain(%d) combination(%d) binary(%d) multilevel(%d)\n", MS_PLAIN, MS_COMBINATION, MS_BINARY, MS_MULTILEVEL);
		exit (-1);
	}

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
			root = PickRandomTerminal(g,random); 
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
		ReportResults(stdout, "bb", times, primal, bestknown);
	}

	static void ReportResults (FILE *file, const string &prefix, double seconds, EdgeCost solvalue, EdgeCost bestknown) {
		fprintf (file, "%ssolution %.20f\n", prefix.c_str(), (double)solvalue);
		fprintf(file, "%stimeus %.3f\n", prefix.c_str(), 1000000.0 * seconds);
		fprintf(file, "%stimems %.6f\n", prefix.c_str(), 1000.0 * seconds);
		fprintf(file, "%stimes %.9f\n", prefix.c_str(), seconds);
		double ratio = (double)solvalue / (double)bestknown;
		double error = ratio - 1;
		fprintf(file, "%sratio %.20f\n", prefix.c_str(), ratio);
		fprintf(file, "%serror %.20f\n", prefix.c_str(), error);
		fprintf(file, "%spcterror %.20f\n", prefix.c_str(), 100.0 * error);

	}
	
	static void RunMultistart(Graph &g, int mstype, int msit, EdgeCost &gbestfound, EdgeCost &bestknown, SteinerConfig &config, char *name, GlobalInfo *ginfo, ExecutionLog *executionLogPtr, SteinerSolution *outSolution) {
		double mstime = 0;
		EdgeCost mscost = g.GetFixedCost();
		int elite = -1;
		EdgeCost bestfound = gbestfound;

		if (g.EdgeCount()>0 && g.TerminalCount()>1) {
			fprintf(stderr, "MULTISTARTING WITH %d\n", mstype);
			mscost = INFINITE_COST;

			// If the multistart type is the time bounded combination multistart type, do not do the
			// incremental number of iterations.
			if (mstype == MS_TIMEBOUNDEDCOMBINATION || mstype == MS_TIMEBOUNDEDMULTILEVEL || mstype == MS_TIMEBOUNDEDADAPTIVE) {
				SteinerSolution solution(&g);
				RFWTimer mstimer(true);
				TimeBoundedMultistart(solution, mstype, msit, config.OUTPUT_INCUMBENT ? name : NULL, -1, &config, ginfo, executionLogPtr);
				mstime = mstimer.getTime();
				if (outSolution != nullptr) {
					outSolution->CopyFrom(&solution);
					mscost = solution.GetCost();
				}
			}
			else {
				int doneit = 0;
				for (int maxit = (config.TIME_LIMIT <= 0 ? (msit <= 0 ? 2 : msit) : 2);
					 (msit <= 0 || doneit < msit) && (config.TIME_LIMIT <= 0 || executionLogPtr->timerPtr->getTime() < config.TIME_LIMIT);
					 maxit *= 2)
				{
					SteinerSolution solution(&g);
					int curit = maxit;
					if (msit > 0 && doneit + maxit > msit)
						curit = msit - doneit;
					fprintf(stderr, "RUNNING %d ITERATIONS (%d ALREADY DONE IN %.2f SEC).\n", curit, doneit, executionLogPtr->timerPtr->getTime());
					RFWTimer mstimer(true);
					switch (mstype) {
					case MS_PLAIN: PlainMultistart(solution, curit, -1, &config); break;
					case MS_COMBINATION: CombinationMultistart(solution, curit, elite, config.OUTPUT_INCUMBENT ? name : NULL, -1, &config, ginfo, executionLogPtr); break;
					case MS_BINARY:	BinaryMultistart(solution, curit, name, &config); break;
					case MS_MULTILEVEL: MultilevelMultistart(solution, curit, curit, elite, config.OUTPUT_INCUMBENT ? name : NULL, -1, &config); break;
					};
					doneit += curit;
					mstime += mstimer.getTime();

					if (solution.GetCost() < mscost) {
						fprintf(stderr, "FOUND BETTER SOLUTION. IMPROVING COST FROM %.0f TO %.0f.\n", mscost, solution.GetCost());
						mscost = solution.GetCost();
						if (outSolution != nullptr && (outSolution->GetCost() > solution.GetCost() || outSolution->EdgeCount() == 0)) {
							outSolution->CopyFrom(&solution);
						}
					}
					//if (outSolution != nullptr && (outSolution->GetCost() > solution.GetCost() || outSolution->EdgeCount() == 0)) {
					//	outSolution->CopyFrom(&solution);
					//	fprintf(stderr, "FOUND BETTER SOLUTION. IMPROVING COST FROM %.0f TO %.0f.\n", mscost, solution.GetCost());
					//	mscost = solution.GetCost();
					//}
				}
			}
		} 
		if (mscost < bestfound) bestfound = mscost;

		fprintf (stderr, "Solution is %.12f.\n", (double)mscost);

		/*
		double ratio = (double)mscost / (double)bestknown;
		double error = ratio - 1;
		fprintf (stdout, "ratio %.12f\n", ratio);
		fprintf (stdout, "error %.12f\n", error);
		fprintf (stdout, "pcterror %.12f\n", 100.0 * error);
		*/
		ReportResults(stdout, "ms", mstime, mscost, bestknown); 

		fprintf (stdout, "msiterations %d\n", msit);
		fprintf(stdout, "mstype %d\n", mstype);
		//fprintf (stdout, "mssolution %.0f\n", mscost);
		//fprintf (stdout, "mstimeseconds %.12f\n", mstime);
		gbestfound = bestfound;
	}

	static void Solve (int argc, char **argv) {
		const uint64_t version = 201504021540;
		cout << "version " << version << endl
			<< "precision " << fixed << setprecision(20) << EDGE_COST_PRECISION << endl;

		if (argc < 2) ShowUsage();
		Graph g;
		g.ReadSTP(argv[1]);
		char *filename = argv[1];

		char *name = new char[strlen(filename)+1];
		ExtractFilename(filename, name); 
		OutputGraphStats (stdout, g, argv[1]);
		fprintf (stdout, "graph %s\n", name);
		EdgeCost bestknown = ReadBound(stdout, name); 
		EdgeCost bestfound = INFINITE_COST;

		//ReadBound("bounds.txt", );

		EdgeCost solcost = 0;
		bool MSTPRUNE = false;
		bool PREPROCESS = false;
		bool SAFE_PREPROCESS = false;
		bool DUALASCENT = false;
		bool BRANCHBOUND = false;
		bool LOCALSEARCH = false;
		bool MULTISTART = false;
		bool BINARYMULTISTART = true;
		int mstype = MS_COMBINATION;
		int msit = 0;

		int seed = 17;
		EdgeCost primal = INFINITE_COST;

		//fprintf (stderr, "ARGC is %d\n", argc);
		SteinerConfig config;

		for (int i=2; i<argc; i+=2) {
			if (i == argc-1) ShowUsage();
			if (strcmp(argv[i], "-ub")==0) {
				primal = atoi(argv[i+1]);
				fprintf (stderr, "Setting upper bound to %.0f.\n", primal);
				continue;
			}

			if (strcmp(argv[i], "-bb")==0) {
				if (atoi(argv[i+1])==0) {
					BRANCHBOUND = false;
				} else {
					BRANCHBOUND = true;
				}
				continue;
			}

			if (strcmp(argv[i], "-prep")==0) {
				if (atoi(argv[i+1])!=0) {
					fprintf (stderr, "Will do preprocessing.\n");
					PREPROCESS = true;
					if (atoi(argv[i + 1]) > 1)
						SAFE_PREPROCESS = true;
				}
				continue;
			}

			if (strcmp(argv[i], "-ls")==0) {
				if (atoi(argv[i+1])!=0) {
					fprintf (stderr, "Will do localsearch.\n");
					LOCALSEARCH = true;
				}
				continue;
			}

			if (strcmp(argv[i], "-bms")==0) {
				if (atoi(argv[i+1])!=0) {
					fprintf (stderr, "Will do binary.\n");
					BINARYMULTISTART = true;
				}
				continue;
			}

			if (strcmp(argv[i], "-msit")==0) {
				msit = (atoi(argv[i+1]));
				fprintf (stderr, "Will run %d multistart iterations.\n", msit);
				continue;
			}

			if (strcmp(argv[i], "-mstype")==0) {
				mstype = (atoi(argv[i+1]));
				fprintf (stderr, "Will run multistart type %d.\n", mstype);
				assert(mstype != 14); // This is corrupt (BB with TimeBoundedMultistart)
				continue;
			}


			if (strcmp(argv[i], "-seed")==0) {
				seed = atoi(argv[i+1]);
				fprintf (stderr, "Set seed to %d.\n", seed);
				continue;
			}

			//maybe config knows what to do with this parameter
			config.ReadParameter (argv[i], argv[i+1]);
		}

		fflush (stderr);

		if (config.EARLY_STOP_BOUND < 0) {config.EARLY_STOP_BOUND = bestknown;}

		bestfound = primal;
		config.Output(stdout);
		fprintf (stdout, "seed %d\n", seed);

		// if there is a best known solution, output it whenever we find it
		if (config.OUTPUT_THRESHOLD<0 && bestknown>=0) {
			config.OUTPUT_THRESHOLD = bestknown;
		}
		
		RFWTimer timer(true);

		// solution cost log.
		ExecutionLog executionLog(&g, &timer, config.TIME_LIMIT);


		/*
		if (argc>2) {
			primal = atoi(argv[2]);
			fprintf (stderr, "Primal bound is %d.\n", primal);
		}*/

		/*
		if (argc>3) {
			seed = atoi(argv[3]);
			fprintf (stderr, "Seed is %d.\n", seed);
		}*/

		MULTISTART = (msit != 0);

		RFWRandom::randomize(seed);


		bool APPLY_PERTURBATION = false;
		if (APPLY_PERTURBATION) {
			fprintf (stderr, "Applying perturbation... ");
			int m = g.EdgeCount();
			vector<EdgeCost> pertcost (m+1,-1);
			RFWLocalRandom random(seed+17);
			ApplyPerturbation(g, pertcost, random, 1, 1.0001);
			g.ApplyCosts(pertcost);
			for (int i=1; i<std::min(10, m); i++) {
				fprintf (stderr, "%.5f ", (double)g.GetCost(i));
			}
			fprintf (stderr, "done.\n");
		}



		/*
		if (DUALASCENT) {
			vector<int> dummy;
			DualAscent(g, primal, dummy);

			exit(-1);
		}*/

		if (PREPROCESS) {
			fprintf (stderr, "Should be preprocessing.\n");
			RFWTimer preptime(true);
			RunPreprocessing(g, !SAFE_PREPROCESS);
			fprintf (stdout, "preptime %.6f\n", preptime.getTime());
			fprintf (stdout, "prepfixed %.6f\n", g.GetFixedCost());
			//PrepBottleneck(g);
			//exit(-1);
		}

		bool GUARDED_MULTISTART = false;
		if (mstype > MS_NUMBER) {
			GUARDED_MULTISTART = true;
			mstype = mstype % 10;
			fprintf (stderr, "Will run guarded mode %d.\n", mstype); 
		}

		double second_time = 0;
		double first_time;

		SteinerSolution bestSolution(&g);

		if (MULTISTART && GUARDED_MULTISTART) {
			fprintf (stderr, "There are %d threads, %d processes.\n", omp_get_max_threads(), omp_get_num_procs());

			GlobalInfo ginfo;
			ginfo.fixed = g.GetFixedCost();
			config.DEPTH_LIMIT = 128;
			fprintf (stderr, "Setting depth limit to %d.\n", config.DEPTH_LIMIT);


			omp_set_num_threads(2);
			RFWTimer ftimer(true);
#pragma omp parallel
			{
				Graph thread_g = g;
				EdgeCost thread_bestfound = bestfound;
				fprintf (stderr, "<%d> ", omp_get_thread_num());
				if (omp_get_thread_num() == 0) {
					RunMultistart(thread_g, mstype, msit, thread_bestfound, bestknown, config, name, &ginfo, &executionLog, nullptr);
					fprintf (stderr, "DONE RUNNING MULTISTART AND FOUND %.3f\n", thread_bestfound);
					ginfo.MakeSolved();
				} else if (omp_get_thread_num() == 1) {
					RFWTimer stimer(true);
					fprintf (stderr, "Should be running something smarter here.\n");
					int bbseed = seed;
					if (bbseed > 0) bbseed = -bbseed;
					RunBranchAndBound(thread_g, bbseed, thread_bestfound, thread_bestfound, bestknown, config, &ginfo, &executionLog);
					second_time += stimer.getTime();
					fprintf (stderr, "DONE RUNNING BRANCHING-AND-BOUND!\n");
					if (!ginfo.bbpruned) ginfo.MakeSolved();
					else fprintf (stderr, "Did not find the optimal solution, though.\n");
				}
#pragma omp critical 
				{
					if (thread_bestfound < bestfound) {
						bestfound = thread_bestfound;
						fprintf (stderr, "THREAD %d UPDATED TO %.3f\n", omp_get_thread_num(), bestfound);
					}
				}
			}
#pragma omp barrier
			{
				fprintf (stderr, "Done with this.\n");
			}
			first_time = ftimer.getTime();
			fprintf (stdout, "bbpruned %d\n", ginfo.bbpruned);
			MULTISTART = false;
			BRANCHBOUND = false;
		}


		if (MULTISTART) {
			RunMultistart(g, mstype, msit, bestfound, bestknown, config, name, NULL, &executionLog, &bestSolution);
		}

		if (LOCALSEARCH) {
			int maxit = 10;
			for (int m=0; m<5; m++) {
				double besttime = 99999999;
				fprintf (stderr, "Running %d... ", m); fflush(stderr);
				if (m==3) {
					fprintf (stderr, "WARNING! SKIPPING METHOD 3.\n");
					continue;
				}
				//if (m != 2) continue;
				for (int i=0; i<maxit; i++) {
					RFWTimer timer(true);
					SteinerSolution solution(&g);
					switch (m) {
						case 0: MSTPrim (g, solution); break;				
						case 1: MSTKruskal (g, solution); break;
						case 2: SPH (g, solution, NULL, 1); break;
						case 3: FullBoruvka (g, solution); break;
						case 4: TestLocalSearch (g, solution); break;
					}
					if (m>=2 && MSTPRUNE) {
						MSTPrune(g,solution);
					}
					solcost = solution.GetCost();
					double t = timer.getTime();
					if (t < besttime) besttime = t;
				}
				fprintf (stderr, "Method %d found solution of cost %d in %.3f milliseconds (best of %d runs).\n", m, solcost, besttime * 1000.0, maxit);
				fflush(stderr);
			}
		}

		if (BRANCHBOUND) {
			RunBranchAndBound(g, seed, primal, bestfound, bestknown, config, NULL, &executionLog);
		}
		double walltime = timer.getTime();
		fprintf (stdout, "totalwalltimeseconds %.12f\n", walltime);
		//fprintf (stdout, "totaltimeseconds %.12f\n", walltime + second_time);
		//fprintf (stdout, "bestsolution %.0f\n", bestfound);

		ReportResults(stdout, "total", walltime + first_time, bestfound, bestknown);
		ReportResults(stdout, "totalcpu", walltime + second_time, bestfound, bestknown);

		// Dump official output logs.
		if (!config.LOG_FILENAME.empty()) {
			ofstream logFile(config.LOG_FILENAME.c_str());
			if (logFile.is_open()) {
				logFile << "SECTION Comment" << endl
					<< "Name \"" << name << "\"" << endl
					<< "Problem \"SPG\"" << endl
					<< "Program \"puw\"" << endl
					<< "Version \"" << version << "\"" << endl
					<< "End" << endl
					<< endl
					<< "SECTION Solutions" << endl;
				for (size_t i = 0; i < executionLog.solCost.size(); ++i) {
					logFile << "Solution " << fixed << executionLog.solCost[i].second << " " << fixed << executionLog.solCost[i].first << endl;
				}
				logFile << "End" << endl
					<< endl
					<< "SECTION Run" << endl
					<< "Threads 1" << endl
					<< "Time " << fixed << walltime << endl
					<< "Dual 0" << endl;
				if (executionLog.solCost.empty())
					logFile << "Primal inf" << endl;
				else
					logFile << "Primal " << fixed << executionLog.bestSolution.GetCost() << endl;
				logFile << "End" << endl
					<< endl;

				if (!executionLog.solCost.empty()) {
					logFile << "SECTION Finalsolution" << endl;

					size_t numVertices = 0;
					for (int v = 1; v <= g.VertexCount(); ++v) {
						if (executionLog.bestSolution.GetDegree(v) > 0 || g.IsTerminal(v))
							++numVertices;
					}
					logFile << "Vertices " << numVertices << endl;
					for (int v = 1; v <= g.VertexCount(); ++v) {
						if (executionLog.bestSolution.GetDegree(v) > 0 || g.IsTerminal(v))
							logFile << "V " << v << endl;
					}
					logFile << "Edges " << executionLog.bestSolution.EdgeCount() << endl;
					for (size_t e = 1; e <= g.EdgeCount(); ++e) {
						if (executionLog.bestSolution.Contains(e))
							logFile << "E " << g.GetFirstEndpoint(e) << " " << g.GetSecondEndpoint(e) << endl;
					}
					logFile << "End" << endl;
				}
 				logFile.close();
			}
		}

		delete [] name;
	}

	static void CombineSolutions(SteinerSolution &target, SteinerSolution &sa, SteinerSolution &sb, RFWLocalRandom &random, SteinerConfig *config) {
		Graph &g = *sa.g;
		int n = g.VertexCount();
        int m = g.EdgeCount();

		const bool verbose = false;


		/*
            UniverseSet baselist = new UniverseSet(n);
            VoronoiData voronoi = new VoronoiData(n);
            UnionFind uf = new UnionFind(n);
            SteinerSolution solution = new SteinerSolution(g);
            BinaryHeap<ArcCost> heap = new BinaryHeap<ArcCost>(n);
            ArcCost [] pertcost = new ArcCost [m+1];

            SteinerSolution target = new SteinerSolution(g);
            Console.Error.Write("{0} x {1}:", sa.GetCost(), sb.GetCost());
			*/

		vector<EdgeCost> pertcost(m+1,-1);

		//was: 1000, (100,500), 1

        for (int e = 1; e <= m; e++) {
			int tcount = 0;
            if (sa.Contains(e)) tcount++;
            if (sb.Contains(e)) tcount++;

            int mult = 1;
            if (tcount == 0) mult = 1000; //random.GetInteger(200,300); //edge in neither solution: very expensive
            else if (tcount == 1) mult = random.GetInteger(100, 500);  //split edge: intermediate cost
			else { mult = 1; } //random.GetInteger(100,200); } //edge in both: keep it
            pertcost[e] = g.GetCost(e) * mult; // g.GetCost(a);
        }

		int root = PickRandomTerminal(g, random);
		SPH (g, target, &pertcost[0], root);
		MSTPrune(g,target);
		RunLocalSearch(g, target, random, -1, config);

		if (verbose) fprintf (stderr, "%d x %d -> %d\n", sa.GetCost(), sb.GetCost(), target.GetCost());


		/*
            baselist.Reset();
            for (int v = 1; v <= n; v++) {if (g.IsTerminal(v)) baselist.Insert(v);}
            ComputeVoronoi(voronoi, baselist, heap, pertcost);

            uf.Reset();
            target.Reset();
            Boruvka(target, voronoi, uf, pertcost);

            ArcCost borcost = target.GetCost();
            FullLocalSearch(target);

            ArcCost newcost = target.GetCost();
            Console.Error.WriteLine(" {0}", newcost);

            return target;
        }*/
	}



	static void BinaryMultistart (SteinerSolution &solution, int maxit, char *outprefix, SteinerConfig *config) {
		int nlevels = 32;

		vector<SteinerSolution*> levelsol (nlevels, NULL); //solution[i]: solution obtained by combining 2^i solutions
		Graph &g = *solution.g;
		int n = g.VertexCount();
		int m = g.EdgeCount();
		RFWTimer timer(true);
		RFWLocalRandom random (RFWRandom::getInteger(1,999999999));

		EdgeCost bestcost = 999999999;
		
		const bool verbose =false;
		bool USE_PERTURBATION = true;
		bool ADAPTIVE_PERTURBATION = false;
		fprintf (stderr, "Running multistart for %d iterations and %d levels (perturbation=%d).\n", maxit, nlevels, USE_PERTURBATION);
		fflush(stderr);
		vector<EdgeCost> pertcost (m+1,-1);

		//bool USE_PERTURBATION = true;
		//bool ADAPTIVE_PERTURBATION = false;
		SteinerSolution cursol(&g);
		SteinerSolution bestsol(&g);
		SteinerSolution combsol(&g);
		//EdgeCost bestcost = 999999999;
		
		for (int i=0 ; i<maxit; i++) {
			int root = random.GetInteger(1,n); //PickRandomTerminal(g, random);
			if (USE_PERTURBATION) {
				InitPerturbation(g, pertcost, random, config);
			}
			SPH (g, cursol, USE_PERTURBATION ? &pertcost[0] : NULL, root);
			MSTPrune(g,cursol);
			RunLocalSearch(g, cursol, random, -1, config);

			if (cursol.GetCost() < bestcost) {
				bestcost = cursol.GetCost();
				bestsol.CopyFrom(&cursol);
			}

			int j;
			for (j=0; j<nlevels; j++) {
				fprintf (stderr, "%10d : %2d : %.0f : ", (int)i, j, bestcost);
				if (levelsol[j] == NULL) {
					fprintf (stderr, "+%.0f\n", cursol.GetCost());
					levelsol[j] = new SteinerSolution(&cursol);
					break;
				}

				// so there is a solution at level j
				
				//SteinerSolution *refsol = elite.GetReference(random.GetInteger(1, elite.Count()));
				int maxtries = 5;
				int t;
				for (t=0; t<maxtries; t++) {
					CombineSolutions(combsol, cursol, *levelsol[j], random, config);
					//fprintf (stderr, "%.0f x %.0f -> %.0f\n", cursol.GetCost(), levelsol[j]->GetCost(), combsol.GetCost());
					if (combsol.GetCost() < cursol.GetCost() && combsol.GetCost() < levelsol[j]->GetCost()) {
						break; //cursol.CopyFrom(&combsol);
					}
				}
				fprintf (stderr, "<%d>", t);

				if (combsol.GetCost() > cursol.GetCost()) {
					combsol.CopyFrom(&cursol);
					fprintf (stderr, "a");
					//fprintf (stderr, "!");
				}

				if (combsol.GetCost() > levelsol[j]->GetCost()) {
					combsol.CopyFrom(levelsol[j]);
					fprintf (stderr, "b");
				}


				cursol.CopyFrom(&combsol);
				delete levelsol[j]; 
				levelsol[j] = NULL;
				//if (cursol.GetCost() < cursol.

				if (cursol.GetCost() < bestcost) {
					bestcost = cursol.GetCost();
					bestsol.CopyFrom(&cursol);
				}
			}

			if (i%10==0) fflush(stderr);

			if (j==nlevels) {
				fprintf (stderr, "Ran out of levels!\n");
				break;
			}


		}


		for (int i=0; i<nlevels; i++) {
			if (levelsol[i]) delete levelsol[i];
		}

		solution.CopyFrom(&bestsol);
	}



	static void EliteMultistart (SteinerSolution &solution, int maxit, int capacity, char *outprefix, SteinerConfig *config) {
		Graph &g = *solution.g;

		SteinerSolution bestsol(&g);
		SteinerSolution combsol(&g);
		EdgeCost bestcost = INFINITE_COST;

		RFWLocalRandom random (RFWRandom::getInteger(1,999999999));

		SolutionPool elite(maxit);

		for (int i=0; i<999999; i++) {
			CombinationMultistart(solution, maxit, capacity, NULL, 0);
			EdgeCost curcost = solution.GetCost();

			fprintf (stderr, "Should be adding solution %.0f to capacity %d.\n", (double)curcost, capacity);
			fflush(stderr);
			int pos1 = elite.Add(&solution);
			//int pos1 = -1;
			CascadedCombination(solution, combsol, elite, -1, random, config);
			curcost = solution.GetCost();
			int pos2 = elite.Add(&solution);

			fprintf (stderr, "[%d,%d:%.0f] ", pos1, pos2, (double)solution.GetCost());

			if (curcost < bestcost) {
				bestsol.CopyFrom(&solution);
				bestcost = curcost;
			}

			fprintf (stderr, "Iteration %d: %.0f\n", i, (double)bestsol.GetCost());

			elite.Output(stderr, 8);
		}

		solution.CopyFrom(&bestsol);


	}

	//--------------------------
	// Repeatedly combine initial solution with others from the elite set, then add the result to elite set itself.
	// Combinations continue until the algorithm fails to improve 'maxfail' times.
	static void CascadedCombination(SteinerSolution &solution, SteinerSolution &combsol, SolutionPool &elite, int maxfail, RFWLocalRandom &random, SteinerConfig *config) {
		if (maxfail < 0) maxfail = config->MAX_COMB_FAIL;
		int failures_to_go = maxfail;
		const bool verbose = false;

		fprintf (stderr, "f%d ", maxfail);

		
		//bool verbose = true;

		//fprintf (stderr, "Adding original solution...\n");
		//elite.Add(&solution); //add initial solution to the pool
		//global.UpdateSolution(solution.GetCost()); //make sure we remember it's the best 

		//if (i >= COMBINATION_THRESHOLD) {
		//fprintf (stderr, "<%d> ", i);
			//fprintf (stderr, "Trying various combinations...\n");
		if (verbose) fprintf (stderr, "%d->", solution.GetCost());
		while (failures_to_go > 0) {
			SteinerSolution *refsol = elite.GetReference(random.GetInteger(1, elite.Count()));
			CombineSolutions(combsol, solution, *refsol, random, config);
			//fprintf (stderr, "<%.0f> ", solution.GetCost());
			if (!combsol.IsBetter(&solution)) {
				failures_to_go --;
				//if (verbose) fprintf (stderr, "! ");
			} else {
				solution.CopyFrom(&combsol);
				//if (verbose) fprintf(stderr, "  ");
			}
		}
		if (verbose) fprintf (stderr, "%d\n", solution.GetCost());
		//KeyVertexInsertion(solution);

		elite.Add(&solution);
	}

	static void AddSmallPerturbation (Graph &g) {
		fatal ("deprecated function");
		int m = g.EdgeCount();
		vector<EdgeCost> pertcost(m+1);
		for (int e=1; e<=m; e++) {
			EdgeCost c = 10000 * g.GetCost(e) + RFWRandom::getInteger(0,10);
			pertcost[e] = c;
		}
		g.ApplyCosts(pertcost);

	}

	struct DistancePair {
		int source;
		EdgeCost distance;
	};

	class ClosenessData {
	private: 
		int k;
		int maxid;

		void GetBounds (int v, int &first, int &last) {
			first = v*maxid;
			last = (v+1)*maxid;
		}
	public:
		vector<DistancePair> data;

		ClosenessData(int _k, int _maxid) {
			k = _k;
			maxid = _maxid;

			data.resize(k*maxid+1);			
		}

	};

	static void FindKClosest (Graph &g, int k, vector<int> sources, ClosenessData &cdata) {
		int n = g.VertexCount();

		for (int v=1; v<=n; v++) {
			//cdata[
		}

	}



	static void KeyVertexInsertion (SteinerSolution &solution, RFWLocalRandom &random) {
		//return;
		Graph &g = *solution.g;
		int n = g.VertexCount();

		SteinerSolution tempsol(&g);

		RFWStack<int,true> tempterm(n+1);

		//fprintf (stderr, "Should be finding improvements!\n");

		fprintf (stderr, "k");
		const bool MOVE_SIDEWAYS = true;

		vector<int> perm(n+1);
		for (int v=1; v<=n; v++) {perm[v] = v;}
		for (int i=1; i<n; i++) {
			int j = random.GetInteger(i,n);
			std::swap(perm[i], perm[j]);
		}

		int improvements = 1;
		while (improvements > 0) {
			improvements = 0;
			//for (int v=1; v<=n; v++) {
			for (int p=1; p<=n; p++) {
				int v = perm[p]; //RFWRandom::getInteger(1,n); //warning! Should have a real permutation.
				//fprintf (stderr, "%d ", v);
				if (g.IsTerminal(v)) continue;
				//if (!solution.Contains(v)) continue; //MUCH WEAKER VERSION OF THE SEARCH
				/*
				if (solution.GetDegree(v) <= 2) {
					//fprintf (stderr, ".");
					continue; //WARNING: THIS IS FOR TESTING ONLY
				}*/
				//fprintf (stderr, "v%d ", v);
				tempterm.reset();
				//fflush (stderr);

				// push all current terminals
				for (int w=1; w<=n; w++) {
					if (g.IsTerminal(w)) continue;
					if (solution.GetDegree(w) <= 2) continue;
					//fprintf (stderr, "t%d:%d ", w, solution.GetDegree(w));
					tempterm.push(w);
				}
				//fprintf (stderr, "Created %d new terminals.\n", tempterm.getNElements());

				// push v itself (if not already pushed)
				//if (solution.GetDegree(v) <= 2) tempterm.push(v);
				int oldt = g.TerminalCount();
				//fprintf (stderr, "oldt=%d ", g.TerminalCount());

				for (int i=1; i<=tempterm.getNElements(); i++) {
					int w = tempterm.peek(i);
					//fprintf (stderr, "x");
					if (g.IsTerminal(w)) fatal ("something wrong!\n");
					g.MakeTerminal(w);
					if (g.TerminalCount() == oldt) fatal ("insertion did nothing");
				}

				//fprintf (stderr, "newt=%d ", g.TerminalCount());

				tempsol.Reset();
				//fprintf (stderr, "%.0f+%.0f = ", solution.GetCost(), tempsol.GetCost());
				SPH(g, tempsol, NULL, v);
				//fprintf (stderr, "%.0f>>%.0f ", oldcost, newcost);
				for (int i=1; i<=tempterm.getNElements(); i++) {
					int w = tempterm.peek(i);
					g.UnmakeTerminal(w);
				}
				MSTPrune(g, tempsol);
				EdgeCost oldcost = solution.GetCost();
				EdgeCost newcost = tempsol.GetCost();

				//if (newcost - oldcost > 0) fprintf (stderr, "%.0f ", newcost - oldcost);
				//fflush(stderr);

				double improvement = oldcost - newcost;

				// remember the new solution if there is an improvement or a tie (if MOVE_SIDEWAYS is true)
				if (improvement > EDGE_COST_PRECISION || (MOVE_SIDEWAYS && (improvement > -EDGE_COST_PRECISION))) {
					solution.CopyFrom(&tempsol);
					if (improvement > EDGE_COST_PRECISION) {
						fprintf (stderr, ".");
						improvements ++;
						fprintf (stderr, "i%.0f ", newcost);
						fflush(stderr);
						return;
					}
				}
				
			}
		}	
	}

	static void DecayLocalSearch (SteinerSolution &solution, vector<EdgeCost> &original, RFWLocalRandom &random, int decaysteps, double exponent, SteinerConfig *config, ExecutionLog *executionLogPtr = nullptr) {
		Graph &g = *solution.g;
		int m = g.EdgeCount();
		vector<EdgeCost> current(m+1);

		double precision = EDGE_COST_PRECISION;

		int decaytogo = decaysteps;
		for (;;) {
			EdgeCost prevcost = solution.GetCost();

			// run one round of local seach
			MSTPrune(g,solution);
			RunLocalSearch(g, solution, random, 1, config, executionLogPtr);

			EdgeCost newcost = solution.GetCost();
			if (decaytogo == 0) break;
			if (prevcost - newcost <= precision) {
				fprintf (stderr, "<%d> ", decaysteps - decaytogo);
				break;
			}
			decaytogo --;

			g.RetrieveCosts(current);
			for (int e=1; e<=m; e++) {
				//current[e] = original[e] + (current[e] - original[e]) * exponent; //NOT REALLY AN EXPONENT
				current[e] = original[e] + (current[e] - original[e]) * exponent; //NOT REALLY AN EXPONENT
			}
			g.ApplyCosts(current);
			solution.UpdateCost();
		}

		g.ApplyCosts(original);
		solution.UpdateCost();
		EdgeCost before = solution.GetCost();
		//fprintf (stderr, "d%.0f->", solution.GetCost());
		//SPH (g, solution, NULL, root);
		//fprintf (stderr, "[%d->", solution.GetCost());
		MSTPrune(g,solution);
		RunLocalSearch(g, solution, random, -1, config);		
		EdgeCost after = solution.GetCost();
		fprintf (stderr, " %3.0f", before - after);
		//fprintf (stderr, "%d]", solution.GetCost());
	}

	static void GenerateRandomizedSolution(SteinerSolution &solution, int root, vector<EdgeCost> &pertcost, RFWLocalRandom &random, SteinerConfig *config, ExecutionLog *executionLogPtr = nullptr) {
		Graph &g = *solution.g;
		int m = g.EdgeCount();
		const bool VERBOSE_STEP = false;
				
		vector<EdgeCost> original (m+1);
		//remember original costs
		g.RetrieveCosts(original);

		// find local optimum with perturbed costs
		g.ApplyCosts(pertcost);
		SPH (g, solution, NULL, root);
		if (executionLogPtr != nullptr) {
			fprintf(stderr, "ADDING SPH SOLUTION WITH COST %f.\n", solution.GetCost());
			executionLogPtr->AddSolution(solution);
		}
		if (VERBOSE_STEP) {
			fprintf (stderr, "Done after SPH (%d).\n", solution.GetCost());
			fflush(stderr);
		}

		int LS_PERT_ROUNDS= std::max(999,g.VertexCount()); 
		double LS_PERT_EXPONENT = 1.0; //no decay
		if (config) {
			LS_PERT_ROUNDS = config->LS_PERT_ROUNDS;
			LS_PERT_EXPONENT = config->LS_PERT_EXPONENT; 
		}

		if (LS_PERT_EXPONENT <= 0) fatal ("invalid perturbation exponent");

		bool DECAY_LOCAL_SEARCH = (LS_PERT_EXPONENT <= 0.999999); 

		if (DECAY_LOCAL_SEARCH) {
			DecayLocalSearch(solution, original, random, LS_PERT_ROUNDS, LS_PERT_EXPONENT, config, executionLogPtr);
		} else {
			MSTPrune(g,solution);
			if (VERBOSE_STEP) {
				fprintf (stderr, "Done after MSTPrune (%d).\n", solution.GetCost());
				fflush(stderr);
			}
			RunLocalSearch(g, solution, random, LS_PERT_ROUNDS, config, executionLogPtr); //, 2);
			if (VERBOSE_STEP) {
				fprintf (stderr, "Done after LocalSearch (%d).\n", solution.GetCost());
				fflush(stderr);
			}
		
			// find real local optimum
			g.ApplyCosts(original);
			solution.UpdateCost();
			//SPH (g, solution, NULL, root);
			//fprintf (stderr, "[%d->", solution.GetCost());
			MSTPrune(g,solution);
			RunLocalSearch(g, solution, random, -1, config);
			//fprintf (stderr, "%d]", solution.GetCost());
		}
	}


	static int ComputeCapacity(int maxit, SteinerConfig *config) {
		double denominator = 1.0;
		if (config) denominator = config->ELITE_DENOMINATOR;
		int capacity = (int)ceil(sqrt((double)maxit / denominator));
		return capacity;

	}



	static void PlainMultistart (SteinerSolution &solution, int maxit, int capacity, SteinerConfig *config) {
		if (capacity<0) capacity = ComputeCapacity(maxit, config);
		SolutionPool elite (capacity);
		RFWLocalRandom random (RFWRandom::getInteger(1,999999999));
		FlexibleMultistart (solution, maxit, elite, NULL, maxit, config);
		//CombinationMultistart (solution, maxit, capacity, NULL, maxit, config);
		SolutionPool children (capacity);

		SolutionPool *a, *b;

		a = &elite;
		b = &children;
		fprintf (stderr, "There are %d elite solutions.\n", elite.GetCount());

		int maxfail = 100;
		int failures = maxfail;

		EdgeCost curbest = elite.FindBestCost();
		//fprintf (stderr, "Initial best is %.0f\n", curbest);
		//exit(-1);

		while (1) {
			fprintf (stderr, "HERE!\n");
			fflush(stderr);
			RecombineGeneration(*a, *b, random, config);
			//if (b->FindBestCost() >= a->FindBestCost()) break;

			//fprintf (stderr, "Doing stuff now.\n");
			//fflush (stderr);

			//EdgeCost newbest = curbest +1;
			EdgeCost newbest = b->FindBestCost();
			if (newbest >= curbest) {
				failures --;
				if (failures == 0) break;
			} else {
				failures = maxfail;
				curbest = newbest;
			}
			fprintf (stderr, "New best solution is %.0f\n", curbest);
			//break;

			//fprintf (stderr, "Swapping (%d,%d)...\n", a->GetCount(), b->GetCount());
			//fflush(stderr);
			std::swap(a,b);
			//fprintf (stderr, "Resetting...\n");
			//fflush(stderr);
			b->HardReset();
 		}

		fprintf (stderr, "Ended with solution with cost %.0f\n", curbest);
	}

	static void CombinationMultistart(SteinerSolution &solution, int maxit, int capacity, char *outprefix, int COMBINATION_THRESHOLD = -1, SteinerConfig *config = NULL, GlobalInfo *ginfo = NULL, ExecutionLog *executionLogPtr = nullptr) {
		if (capacity<0) capacity = ComputeCapacity(maxit, config);
		SolutionPool elite (capacity);
		if (!config->AGGRESSIVE_COMBINATION) COMBINATION_THRESHOLD = capacity;
		fprintf (stderr, "<<<<< %p >>>>>>\n", ginfo);
		FlexibleMultistart (solution, maxit, elite, outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr);
	}


	static void TimeBoundedMultistart(SteinerSolution &solution, int mstype, int maxitupper, char *outprefix, int COMBINATION_THRESHOLD = -1, SteinerConfig *config = NULL, GlobalInfo *ginfo = NULL, ExecutionLog *executionLogPtr = nullptr) {
		SolutionPool tentativeElite(2);
		RFWTimer tentativeTimer(true);
		FlexibleMultistart(solution, 1, tentativeElite, outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr, false, false);
		double tentativeTime = tentativeTimer.getTime();

		int maxit = maxitupper;
		const double blowupFactor = 2.5; // this fell from the sky after careful consideration
		if (tentativeTime > 0 && config->TIME_LIMIT > 0)
			maxit = static_cast<int>(ceil(config->TIME_LIMIT / tentativeTime / blowupFactor));
		if (maxit > maxitupper)
			maxit = maxitupper;
		fprintf(stderr, "TENTATIVE TIME: %.3f SEC; WILL COMPUTE %d ITERATIONS (UPPER BOUND = %d).\n", tentativeTime, maxit, maxitupper);

		int capacity = ComputeCapacity(maxit, config);
		SolutionPool elite(capacity);
		if (!config->AGGRESSIVE_COMBINATION) COMBINATION_THRESHOLD = capacity;

		// Determine which multistart to call depending on mstype.
		int localtype = mstype;
		if (mstype == MS_TIMEBOUNDEDADAPTIVE) {
			localtype = maxit > 2048 ? MS_TIMEBOUNDEDMULTILEVEL : MS_TIMEBOUNDEDCOMBINATION;
		}
		if (localtype == MS_TIMEBOUNDEDCOMBINATION) {
			FlexibleMultistart(solution, maxitupper, elite, outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr);
		}
		else if (localtype == MS_TIMEBOUNDEDMULTILEVEL) {
			MultilevelMultistart(solution, maxit, maxitupper, capacity, outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr);
		}
		else {
			fatal("Not supported.");
		}

		if (tentativeElite.FindBestCost() < solution.GetCost()) {
			solution.CopyFrom(tentativeElite.GetReference(tentativeElite.FindBestPosition()));
			fprintf(stderr, "USING SOLUTION FROM FIRST TENTATIVE ITERATION.\n");
		}

		fprintf(stdout, "actualmstype %d\n", localtype);
	}


	template<class T> static void Permute(vector<T> &array, RFWLocalRandom &random) {
		int s = (int)array.size();
		for (int i=0; i<s-1; i++) {
			int j = random.GetInteger(i,s-1);
			std::swap(array[i], array[j]);
		}
	}

	static void MultilevelMultistart(SteinerSolution &solution, int maxit, int maxitupper, int capacity, char *outprefix, int COMBINATION_THRESHOLD = -1, SteinerConfig *config = NULL, GlobalInfo *ginfo = NULL, ExecutionLog *executionLogPtr = nullptr) {
		fprintf (stderr, "SHOULD BE RUNNING MMS FROM SOLUTION %.0f\n", solution.GetCost());
		fflush(stderr);
		if (capacity<0) capacity = ComputeCapacity(maxit, config);
		const int SUBCOUNT = 4;

		SolutionPool *subelite[SUBCOUNT];
		if (!config->AGGRESSIVE_COMBINATION) COMBINATION_THRESHOLD = capacity;

		int localit = maxit / (2*SUBCOUNT);
		int subcapacity = ComputeCapacity(localit, config);
		int restit = maxit - SUBCOUNT*localit;
		int restcap = ComputeCapacity(restit, config);
		SolutionPool elite(restcap);

		// Only do phase one of the algorithm if there are enough iterations.
		if (localit > 0) {
			for (int i = 0; i < SUBCOUNT; i++) {
				subelite[i] = new SolutionPool(subcapacity);
			}

			// Runs subcount independent multistarts using half the total number of iterations.
			for (int i = 0; i < SUBCOUNT; i++) {
				fprintf(stderr, "SUBPROBLEM %d\n", i);
				FlexibleMultistart(solution, localit, *subelite[i], outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr);
				fprintf(stderr, "\n\n");
			}

			vector<SteinerSolution*> solpointers;

			for (int i = 0; i < SUBCOUNT; i++) {
				subelite[i]->Output(stderr, 8);
				fprintf(stderr, "\n");
				int count = subelite[i]->GetCount();
				for (int j = 1; j <= count; j++) solpointers.push_back(subelite[i]->GetReference(j));
			}

			RFWLocalRandom random(RFWRandom::getInteger(1, 1000000));
			Permute(solpointers, random);
			//std::sort(solpointers.begin(), solpointers.end(), [&](SteinerSolution *x, SteinerSolution *y) {return x->GetCost() >= y->GetCost();});
			//sort (&perm[0], &perm[perm.size()], [&](int x, int y) {return totalbound[x]/count[x] > totalbound[y]/count[y];});

			for (int i = 0; i < (int)solpointers.size(); i++) {
				SteinerSolution *sol = solpointers[i];
				if (sol->GetCost() < solution.GetCost()) { solution.CopyFrom(sol); }
				elite.Add(sol);
			}

			for (int i = 0; i<SUBCOUNT; i++) {
				delete subelite[i];
			}

			/*
			for (int i=0; i<SUBCOUNT; i++) {
			subelite[i]->Output(stderr, 8);
			fprintf (stderr, "\n");
			int count = subelite[i]->GetCount();
			for (int j=1; j<=count; j++) {
			SteinerSolution *sol = subelite[i]->GetReference(j);
			if (sol->GetCost() < solution.GetCost()) {solution.CopyFrom(sol);}
			elite.Add(sol);
			}
			}*/
			elite.Output(stderr, 8);

		}


		fprintf (stderr, "SUPERPROBLEM (from solution %.0f):\n", solution.GetCost());
		FlexibleMultistart(solution, maxitupper - SUBCOUNT*localit, elite, outprefix, COMBINATION_THRESHOLD, config, ginfo, executionLogPtr);
	}

	static void RecombineGeneration (SolutionPool &parents, SolutionPool &children, RFWLocalRandom &random, SteinerConfig *config) {
		int pcount = parents.GetCount();
		if (pcount == 0) return;

		fprintf (stderr, "Should be combining.\n");

		Graph &g = *(parents.GetReference(1)->g);

		SteinerSolution newsol(&g);
		SteinerSolution bestsol(&g);
		bestsol.CopyFrom(parents.GetReference(parents.FindBestPosition()));
		//fprintf (stderr, "WARNING: THIS IS VERY WRONG.\n");

		for (int i=1; i<pcount; i++) {
			SteinerSolution *a = parents.GetReference(i);
			fprintf (stderr, " %.0f", a->GetCost());
			for (int j=i+1; j<=pcount; j++) {
				SteinerSolution *b = parents.GetReference(j);
				CombineSolutions(newsol, *a, *b, random, config);
				//fprintf (stderr, "%.0f x %.0f: %.0f\t", a->GetCost(), b->GetCost(), newsol.GetCost());
				if (newsol.GetCost() < bestsol.GetCost()) {
					bestsol.CopyFrom(&newsol);
				}
				children.Add(&newsol);
			}
		}

		// make sure the best solution (even if from a previous iteration) is preserved
		children.Add(&bestsol);

		fprintf (stderr, "Best solution is %.0f\n", bestsol.GetCost());
	}

	static void FlexibleMultistart(SteinerSolution &solution, int maxit, SolutionPool &elite, char *outprefix, int COMBINATION_THRESHOLD = -1, SteinerConfig *config = NULL, GlobalInfo *ginfo = NULL, ExecutionLog *executionLogPtr = nullptr, bool USE_PERTURBATION = true, bool outputStats = true) {
        Graph &g = *solution.g;

		//AddSmallPerturbation(g);
		int itbits = 6;
		
		//int maxit = 1 << itbits;
        //int capacity = 1 << (itbits / 2);
		/*
		if (capacity < 0) {
			double denominator = 1.0;
			if (config) denominator = config->ELITE_DENOMINATOR;
			capacity = (int)ceil(sqrt((double)maxit / denominator));
		}*/
		int n = g.VertexCount();

		SteinerSolution bestsol(&g); //best solution found so far
		SteinerSolution combsol(&g); //combined solution
		bestsol.CopyFrom(&solution);
		//SolutionPool elite (capacity);
		int capacity = elite.GetCapacity();

		//Graph &g = *solution.g;

		int m = g.EdgeCount();
		RFWTimer timer(true);
		RFWLocalRandom random (RFWRandom::getInteger(1,999999999));

		EdgeCost bestcost = INFINITE_COST;
		if (elite.GetCount() > 0) {
			int p = elite.FindBestPosition();
			bestsol.CopyFrom(elite.GetReference(p));
			bestcost = bestsol.GetCost();
		}
		
		const bool verbose = false;

		//bool  = perturbation;
		bool ADAPTIVE_PERTURBATION = false;
		bool RESILIENT_PERTURBATION = USE_PERTURBATION;
		//int COMBINATION_THRESHOLD = -1; //capacity; 

		fprintf (stderr, "Running multistart for %d iterations and %d elite solutions (perturbation=%d).\n", maxit, capacity, USE_PERTURBATION);
		//fflush (stderr);
		vector<EdgeCost> pertcost (m+1,-1);

		const bool VERBOSE_STEP = false;

		int i;
		for (i=0; i<maxit; i++) {
			//if (ginfo) fprintf (stderr, "THERE IS A GINFO.\n");
			if (ginfo && ginfo->IsSolved()) {
				fprintf (stderr, "Stopping at iteration %d.\n", i);	
				break;
			}
			if (config->TIME_LIMIT > 0 && executionLogPtr->timerPtr->getTime() >= config->TIME_LIMIT) {
				fprintf(stderr, "Stopping at iteration %d because of time limit of %2.f sec (%.2f sec passed).\n", i, config->TIME_LIMIT, executionLogPtr->timerPtr->getTime());
				break;
			}

			int root = random.GetInteger(1,n); //PickRandomTerminal(g, random);


			if (USE_PERTURBATION) {
				ADAPTIVE_PERTURBATION = false; //((i % 2) == 1); //random.GetInteger(0,1)==0;
				//ADAPTIVE_PERTURBATION = true;
				if (ADAPTIVE_PERTURBATION) {
					AdaptivePerturbation(g, pertcost, elite, random);
				} else {
					//fprintf (stderr, "Here!\n");
					InitPerturbation(g, pertcost, random, config);
					//fflush (stderr);
				}
			}

			if (RESILIENT_PERTURBATION) {
				GenerateRandomizedSolution (solution, root, pertcost, random, config);
			} else {
				SPH (g, solution, USE_PERTURBATION ? &pertcost[0] : NULL, root);
				if (executionLogPtr != nullptr)
					executionLogPtr->AddSolution(solution);
				MSTPrune(g,solution);
				RunLocalSearch(g, solution, random, -1, config, executionLogPtr);
			}

			if (executionLogPtr != nullptr) {
				executionLogPtr->AddSolution(solution);
			}

			//fprintf (stderr, "Adding original solution...\n");
			elite.Add(&solution); //add initial solution to the pool
			//global.UpdateSolution(solution.GetCost()); //make sure we remember it's the best 

			if (i >= COMBINATION_THRESHOLD) {
				CascadedCombination(solution, combsol, elite, -1, random, config);
			} else {
				fprintf (stderr, "! ");
			}


			EdgeCost solcost = solution.GetCost();
			//fprintf (stderr, "Found solution costing %d.\n", solcost);
			if (solcost < bestcost) {
				if (executionLogPtr != nullptr) {
					executionLogPtr->AddSolution(solution);
				}
				bestcost = solcost;
				bestsol.CopyFrom(&solution);
				//fprintf (stderr, "HERE (%.2f %p %p %.2f).\n", bestcost, outprefix, config, config->OUTPUT_THRESHOLD);
				if (outprefix) {
					if (config && bestcost < config->OUTPUT_THRESHOLD) {
						fprintf (stderr, "\n<outputting solution with value %0f>\n", bestcost);
						config->OUTPUT_THRESHOLD = bestcost;
						bestsol.Output(outprefix);
					}
				}
				if (ginfo) ginfo->UpdateBestFound(bestcost);
			}
			bool localverbose = ((i % 10) == 0);
			if (localverbose) {
				fprintf (stderr, "%6d : %6d : %10.0f : %10.5f\n", i, root, (double)solcost, (double)bestcost);
				fflush(stderr);
			}

			if (config && bestcost <= config->EARLY_STOP_BOUND) {
				fprintf (stderr, "Early stop!\n");
				break;
			}


			//elite.Output(stderr);
			//if (i % 10 == 0) fflush (stderr);
		}

		if (outputStats) {
			//fprintf (stdout, "msiterations %d\n", maxit);
			fprintf(stdout, "actualiterations %d\n", i);
			fprintf(stdout, "earlystop %d\n", maxit != i);
			fprintf(stdout, "mselite %d\n", capacity);
			//fprintf (stdout, "totaltimeseconds %.8f\n", timer.getTime());
			//fprintf (stdout, "mssolution %.0f\n", (double)bestcost);
		}
		fflush (stderr);

		solution.CopyFrom(&bestsol);

    }

	class Perturbator {
		int PERTURBATION_MODE;
		double base;
		double range;
		double threshold;
		SteinerConfig *config;
	public:
		void SetParameters (int n, SteinerConfig *_config) {
			config = _config;
			PERTURBATION_MODE = config->PERTURBATION_MODE;
			base = config->PERT_FLOOR;
			range = config->PERT_RANGE;
			threshold = 0;
			if (PERTURBATION_MODE == 2) {
				threshold = sqrt((double)n) / (double)n;
			} else if (PERTURBATION_MODE == 4) {
				threshold = sqrt((double)n) / (double)(2*n);
			} else {
				threshold = (log(n)/log(2)) / (double)n;
			}
		}

		double GetPerturbation(double r) {
			switch (PERTURBATION_MODE) {
				case 0: return base + r * range; break; //expected value is 1.5
				case 1: return base + r + 1/(10.0*r); break; //NOT TRUE EXPECTATION
				case 2: return (r >= threshold) ? 1 + range * r : range/r; break;
				case 3: return base + range * r*r; break;
				case 4: return (r >= threshold) ? 1 + range * r : range/r; break;
				case 5: return base + range * (1.0 - r*r); break;
				case 6: return base + 1.5 * r * r; break; //expected value is 1.5
				case 7: return (r >= threshold) ? 1 + range * r : max(r/threshold, 0.0000001); break;
				default: fatal ("invalid perturbation mode");
			}
			return 0;
		}

		void ResetRange (double r) {
			if (config->PERT_EXTRA > 0) {
				range = config->PERT_RANGE + config->PERT_EXTRA * r;
			}
			/*
			if (VARIABLE_RANGE) {
				if (3*r < 1) range = .25;
				else if (3*r > 2) range = 1;
				else range = 4;
			}*/
		}
	};

	static void VertexPerturbation (Graph &g, vector<EdgeCost> &pertcost, RFWLocalRandom &random, SteinerConfig *config) {
		int n = g.VertexCount();
		int m = g.EdgeCount();

		Perturbator p;
		p.SetParameters(n, config);
		p.ResetRange(random.GetDouble());

		static bool first = true;
		if (first) {
			fprintf (stderr, "USING PERTURBATION MODE %d\n", config->PERTURBATION_MODE);
			first = false;
		}
		vector<double> vpert(n+1);

		for (int v=0; v<=n; v++) {
			vpert[v] = p.GetPerturbation(random.GetDouble());
		}

		// each edge is perturbed by the average of its endpoints
		for (int e = 1; e <= m; e++) {
			int v, w;
            g.GetEndpoints(e, v, w);
            double p = (vpert[v] + vpert[w]) / 2.0;
			EdgeCost cost =g.GetCost(e) * (EdgeCost)p;
            pertcost[e] = cost;
		}
	}



	static void OldVertexPerturbation (Graph &g, vector<EdgeCost> &pertcost, RFWLocalRandom &random) {
		const bool debug = false;
		//fprintf (stderr, "vp");

		int n = g.VertexCount();
        int m = g.EdgeCount();
        int divisor = 1;

        EdgeCost mine = g.GetMinCost();
		while (mine > 1000) {
			divisor *= 10;
            mine /= 10;
		}

		bool UNIFORM = false;

		int t = g.TerminalCount();

		vector<EdgeCost> vpert(n+1);
		for (int v=0; v<=n; v++) {
			if (UNIFORM) {
				// it was 100,200
				vpert[v] = random.GetInteger(100,200); //it's 100,120 for local stuff
				//vpert[v] = random.GetInteger(100, 200); 100 + 20 * g.GetDegree(v));
			} else {
				double r = random.GetDouble();
				//r = (1 - r*r);
				//vpert[v] = 100 + 100 * (r + 1/r);
				vpert[v] = 100*r + 10/r;
				//vpert[v] = 100 + 20*r + 1/r;

				//fprintf (stderr, "%.1f ", vpert[v]);
			}
		}

		for (int e = 1; e <= m; e++) {
			int v, w;
            g.GetEndpoints(e, v, w);
            bool average = true;

            EdgeCost p;
			if (average) {
				p = (vpert[v] + vpert[w]) / 2;
			} else {
				/*
                    ArcCost minp = (vpert[v] < vpert[w]) ? vpert[v] : vpert[w];
                    ArcCost maxp = (vpert[v] > vpert[w]) ? vpert[v] : vpert[w];
                    p = (3*minp + maxp) / 4;
					*/
			}

                //ArcCost p = (vpert[v] < vpert[w]) ? vpert[v] : vpert[w]; 

                //Console.Error.Write("<{0}>", p);

			EdgeCost cost = (g.GetCost(e) / divisor) * p;
            pertcost[e] = cost;
		}
	}

	static void AdaptivePerturbation (Graph &g, vector<EdgeCost> &pertcost, SolutionPool &elite, RFWLocalRandom &random) {

		int m = g.EdgeCount();
		int solcount = elite.Count();

		//fprintf (stderr, "There should be %d elite solutions.\n", solcount);
		//fflush (stderr);

		for (int e=1; e<=m; e++) pertcost[e] = 0;

		// count the number of times each edge is used
		for (int i=1; i<=solcount; i++) {
			SteinerSolution *s = elite.GetReference(i);
			if (!s) fatal ("Something really bad happened.\n");
			for (int e=1; e<=m; e++) {
				//fprintf (stderr, "%d ", e);
				//fflush (stderr);
				if (s->Contains(e)) pertcost[e] ++;
			}
		}

		//fprintf (stderr, "Done computing multiplicities.\n");

		double factor = (90.0 / (double)solcount);

		int maxc = 0;
		for (int e=1; e<=m; e++) {
			int c = (int)pertcost[e];
			if (c > maxc) maxc = c;
			pertcost[e] = g.GetCost(e) * random.GetInteger(100, 110 + (int)ceil((double)c * factor));
		}

		fprintf (stderr, "%d ", maxc, solcount);
	}

    static void InitPerturbation(Graph &g, vector<EdgeCost> &pertcost, RFWLocalRandom &random, SteinerConfig *config) {
		//double p = config->PERT_VERTEX;
		bool USE_VERTEX_PERTURBATION = (random.GetDouble() < config->PERT_VERTEX);

		//bool USE_VERTEX_PERTURBATION = config->PERT_VERTEX; //random.GetInteger(1,2)==1; //
        if (USE_VERTEX_PERTURBATION) {
			VertexPerturbation(g, pertcost, random, config);
			return;
        }
		fprintf (stderr, "e");
		Perturbator p;
		p.SetParameters(g.VertexCount(), config); 
		p.ResetRange(random.GetDouble());

		int m = g.EdgeCount();
		for (int e=1; e<=m; e++) {
			//pertcost[e] = (base + range * random.GetDouble()) * g.GetCost(e);
			pertcost[e] = p.GetPerturbation(random.GetDouble()) * g.GetCost(e);
		}
	}

	static void ApplyPerturbation (Graph &g, vector<EdgeCost> &pertcost, RFWLocalRandom &random, double minfactor, double maxfactor) {
		double range = (maxfactor - minfactor);

		int m = g.EdgeCount();
		for (int e=1; e<=m; e++) {
			double factor = minfactor + range * random.GetDoubleOpen();
			pertcost[e] = (EdgeCost)factor * g.GetCost(e);
		}
	}


	static void Multistart (SteinerSolution &solution, SteinerConfig *config) {
		int maxit = 9999;
		Graph &g = *solution.g;

		int m = g.EdgeCount();
		RFWTimer timer(true);
		RFWLocalRandom random (RFWRandom::getInteger(1,999999999));

		EdgeCost bestcost = INFINITE_COST;
		
		bool USE_PERTURBATION = true;

		fprintf (stderr, "Running multistart for %d iterations with perturbation=%d.\n", maxit, USE_PERTURBATION);

		vector<EdgeCost> pertcost (m+1,-1);

		for (int i=0; i<maxit; i++) {
			int root = PickRandomTerminal(g, random);
			if (USE_PERTURBATION) InitPerturbation(g, pertcost, random, NULL);
			SPH (g, solution, USE_PERTURBATION ? &pertcost[0] : NULL, root);
			MSTPrune(g,solution);
			RunLocalSearch(g, solution, random, -1, config);

			EdgeCost solcost = solution.GetCost();
			if (solcost < bestcost) {
				bestcost = solcost;
			}
			if (i % 10 == 0) {
				fprintf (stderr, "%6d : %6d : %6d : %6d\n", i, root, solcost, bestcost);
				fflush(stderr);
			}
			//if (i % 10 == 0) fflush (stderr);
		}
		fflush (stderr);

		fprintf (stderr, "totaltimeseconds %.8f\n", timer.getTime());
		fprintf (stderr, "solution %d\n", bestcost);

		

		/*
		if (LOCALSEARCH) {
			int maxit = 10;
			for (int m=0; m<5; m++) {
				double besttime = 99999999;
				fprintf (stderr, "Running %d... ", m); fflush(stderr);
				if (m==3) {
					fprintf (stderr, "WARNING! SKIPPING METHOD 3.\n");
					continue;
				}
				//if (m != 2) continue;
				for (int i=0; i<maxit; i++) {
					RFWTimer timer(true);
					SteinerSolution solution(&g);
					switch (m) {
						case 0: MSTPrim (g, solution); break;				
						case 1: MSTKruskal (g, solution); break;
						case 2: SPH (g, solution, NULL, 1); break;
						case 3: FullBoruvka (g, solution); break;
						case 4: TestLocalSearch (g, solution); break;
					}
					if (m>=2 && MSTPRUNE) {
						MSTPrune(g,solution);
					}
					solcost = solution.GetCost();
					double t = timer.getTime();
					if (t < besttime) besttime = t;
				}
				fprintf (stderr, "Method %d found solution of cost %d in %.3f milliseconds (best of %d runs).\n", m, solcost, besttime * 1000.0, maxit);
				fflush(stderr);
			}
		}*/
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



	static void RunLocalSearch (Graph &g, SteinerSolution &solution, RFWLocalRandom &random, int maxrounds, SteinerConfig *config, ExecutionLog *executionLogPtr = nullptr) {
		bool RUN_Q = false;
		bool RUN_V = false;
		bool RUN_U = false;
		bool RUN_K = false;
		bool RESTRICT_K = false; //
		bool verbose = false;
		static bool first = true;
		
		char *lstype = config->LSTYPE;

		bool wait = false;
		for (int i=0; ;i++) {
			char c = lstype[i];
			if (c==0) break;
			if (c=='w') {
				wait = true;
				continue;
			}
			switch (c) {
				case 'k': RUN_K = true; RESTRICT_K = wait; break;
				case 'v': RUN_V = true; break;
				case 'u': RUN_U = true; break;
				case 'q': RUN_Q = true; break;
				default: fprintf (stderr, "WARNING: invalid local search parameter (%c).\n", c);
			}
			wait = false;
		}

		//fprintf (stderr, "<%d%d>", RUN_K, RESTRICT_K);



		//return;

		//int maxrounds = 999;
		if (maxrounds < 0) maxrounds = 999; //999; //large number



		if (first) {
			fprintf (stderr, "Running with maxrounds %d.\n", maxrounds);
		}
		
		//SPH (g, solution, NULL, PickRandomTerminal(g,random));
		//if (verbose) fprintf (stderr, "SPH found solution %d.\n", solution.GetCost());
		int n = g.VertexCount();
		EdgeCost oldcost = solution.GetCost();
		RFWTimer timer(true);
		int rounds = 0;
		int i=0;
		for (i=0; i<maxrounds; i++) {
			rounds ++;
			if (RUN_V) {
				VertexInsertion(g, solution, n, random);
				MSTPrune(g, solution);
				if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
			}
			//fprintf (stderr, "Starting key vertex elimination!\n");
			//fflush (stderr);
			
			if (RUN_Q) {
				KeyVertexElimination(g, solution, random);
				//fprintf (stderr, "Ending key vertex elimination!\n");
				if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
			} 
			
			if (RUN_U) {
				VertexElimination(g, solution, random);
				if (verbose) fprintf (stderr, " u%d ", solution.GetCost());
			}

			if (RUN_K) {
				if (!RESTRICT_K || solution.GetCost() == oldcost) {
					KeyVertexInsertion(solution, random);
				} 
			}

			EdgeCost newcost = solution.GetCost();
			if (newcost > oldcost + EDGE_COST_PRECISION) fatal ("invalid result");
			if (newcost > oldcost - EDGE_COST_PRECISION) break; //stop if result did not improve
			if (executionLogPtr != nullptr) executionLogPtr->AddSolution(solution);
			oldcost = newcost;
		}
		const bool VERBOSE_ROUNDS = false;
		if (VERBOSE_ROUNDS) fprintf (stderr, "%d ", i);
		//fprintf (stderr, "%d ", rounds);
		if (verbose) fprintf (stderr, "\n\nDone with local search: %d (%.2f ms, %.2f ms average)\n", solution.GetCost(), 1000 * timer.getTime(), 1000 * timer.getTime() / (double)rounds); 
		//exit(-1);
		first = false;
	}




	static void TestLocalSearch (Graph &g, SteinerSolution &solution) {
		RFWLocalRandom random(RFWRandom::getInteger(1,999999999));
		bool RUN_Q = true;
		bool RUN_V = true;
		bool RUN_U = false;
		bool verbose = false;
		
		SPH (g, solution, NULL, PickRandomTerminal(g,random));
		if (verbose) fprintf (stderr, "SPH found solution %d.\n", solution.GetCost());
		int n = g.VertexCount();
		EdgeCost oldcost = solution.GetCost();
		RFWTimer timer(true);
		int rounds = 0;
		for (int i=0; i<10; i++) {
			rounds ++;
			if (RUN_V) {
				VertexInsertion(g, solution, n, random);
				MSTPrune(g, solution);
				if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
			}
			//fprintf (stderr, "Starting key vertex elimination!\n");
			//fflush (stderr);
			
			if (RUN_Q) {
				KeyVertexElimination(g, solution, random);
				//fprintf (stderr, "Ending key vertex elimination!\n");
				if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
			} 
			
			if (RUN_U) {
				VertexElimination(g, solution, random);
				if (verbose) fprintf (stderr, " u%d ", solution.GetCost());
			}
			EdgeCost newcost = solution.GetCost();
			if (newcost > oldcost) fatal ("invalid result");
			if (newcost == oldcost) break;
			oldcost = newcost;
		}
		if (verbose) fprintf (stderr, "\n\nDone with local search: %d (%.2f ms, %.2f ms average)\n", solution.GetCost(), 1000 * timer.getTime(), 1000 * timer.getTime() / (double)rounds); 
		//exit(-1);
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

				//fprintf (stderr, " %d:%d->%d ", w, arc->label, newcost); 
				//fflush (stderr);
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
			//fprintf (stderr, "\n%d : %d :   ", v, g.GetDegree(v));
			//fprintf (stderr, "  ", v, g.GetDegree(v));
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
				//fprintf (stderr, " %d:%d", zeroside ? 0 : 1, dual);
				//fflush(stderr);
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
				SPH (g, solution, NULL, PickRandomTerminal(g,random));
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
							VertexInsertion (g, solution, 5*g.VertexCount(), random);
							//MSTPrune(g,solution);
							VertexElimination(g, solution, random);
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
						int moves = VertexElimination(g, solution, random);
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
					if (debug) CheckSolution(g, solution);
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

				//fprintf (stderr, " %d:%d->%d ", w, arc->label, newcost); 
				//fflush (stderr);
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


	// Comput minimum spanning tree of the graph g using Prim's algorithm (ignores terminals)
	// 'solution' will contain the MST edges

	static void MSTPrim (Graph &g, SteinerSolution &solution) {
		bool verbose = false;
		int n = g.VertexCount();
		BinaryHeap<EdgeCost> heap(n); // = new BinaryHeap<ArcCost>(n);
		vector<int> parc (n+1); //not need to initialize
		unsigned int r = PickRandomTerminal(g);
		parc[r] = 0;
		int nscanned = 0;
		solution.Reset();
		//int inscount = 0;

		heap.Insert(r, 0);
		while (!heap.IsEmpty()) {
			unsigned int v;
			EdgeCost acost;
			heap.RemoveFirst(v,acost); //v, out acost);
            if (v!=r) solution.Insert(parc[v]); //add parent edge to the solution

			//scan outgoing arcs
            nscanned ++;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //neighbor
				if (solution.GetDegree(w) > 0) continue; //ignore if already in solution
				if (heap.Insert(w, a->cost)) {
					parc[w] = a->label;
					//inscount ++;
				}
            }
        }
		//fprintf (stderr, "%.3f ", (double)inscount / (double)n);
		//if (nscanned != n) fprintf (stderr, "Warning: graph is not connected");
	}

	//--------------------------------------------------------------
	// Kruskal's algorithm to compute the MST of the full graph.
	// (Could be made faster for some instances with partial sort.)
	//--------------------------------------------------------------
	static void MSTKruskal (Graph &g, SteinerSolution &solution) {
		// create list of all edges sorted in increasing order of weight
		const bool verbose = false;
		int i, m = g.EdgeCount();
		vector<int> elist(m+1);
		for (i=0; i<m; i++) {elist[i] = i+1;}
		sort(&elist[0], &elist[m],  [&](int x, int y) {return g.GetCost(x)<g.GetCost(y);});

		// create empty solution and union find with singletons
		solution.Reset();
		int n = g.VertexCount();
		UnionFind uf(n);
		int togo = n-1; //number of unions left

		// process all edges in order
		for (i=0; i<m; i++) {
			int e = elist[i];
			int v, w;
			g.GetEndpoints(e,v,w);	
			if (uf.Union(v,w)) { //if successfully joined...
				solution.Insert(e); //...we have a new edge in the tree
				if (--togo==0) break;
			}
		}

		if (verbose) fprintf (stderr, "%.3f ", (double)i/(double)m);
	}


    /// Compute the MST of the distance network of the subgraph induced
    /// by the vertices in bases.
	/// <param name="solution">Final solution (output)</param>
    /// <param name="bases">Set of bases (key vertices)</param>
	/*
	static void DNH(Graph &g, SteinerSolution &solution, UniverseSet &baselist) {
            int n = g.VertexCount();
            int m = g.EdgeCount();
            VoronoiData voronoi = new VoronoiData(n);
            UnionFind uf = new UnionFind(n);
            BinaryHeap<ArcCost> heap = new BinaryHeap<ArcCost>(n);

            
            ArcCost [] pertcost = null;
            

            ComputeVoronoi(voronoi, baselist, heap, pertcost);
            solution.Reset();
            Boruvka(solution, voronoi, uf, pertcost);
            //Console.Error.WriteLine("DNH");
        }
		*/

	/// Boruvka-based implementation of DNH (given a Voronoi diagram and the associated union-find data structure).
	/// Traverses a list of edge IDs in each pass, eliminating those that are no longer boundary.
	/// Seems to be worse than the version that actually traverses graphs.

	static void Boruvka(Graph &g, SteinerSolution &solution, VoronoiData &voronoi, UnionFind &uf, EdgeCost *pertcost) {
		int v, n = g.VertexCount();
		int m = g.EdgeCount();
		EdgeCost solvalue = 0;
        const bool verbose = false;

		//count boundary regions in the current diagram
        int nregions = 0;
        for (v=1; v<=n; v++) {
			// it's not clear find is needed, unless some merges happened before
            if (uf.Find(voronoi.GetBase(v))==v) nregions ++;
		}

		//Boruvka's algorithm
        int *minarc = new int[n+1]; //minimum outgoing edge from the region based in v
        EdgeCost *minvalue = new EdgeCost[n+1]; //value associated with the neighbor
		int *elist = new int [m]; //list of all potential boundary edges
		for (int i=0; i<m; i++) elist[i] = (i+1);

		int ecount = m; //number of edges in edge list
        int rounds = 0;
        bool changes = true;

		while (changes && nregions > 1) {
			rounds ++;
			//if (rounds > 3) break;
			changes = false;

			//initially, we don't know what are the arcs out of each component
			for (v=1; v<=n; v++) {
				minarc[v] = -1;
				minvalue[v] = 0;
			}

			int nextpos = 0;
			//fprintf (stderr, "%d ", ecount);
			for (int i=0; i<ecount; i++) {
				int e = elist[i]; //get edge in the current position
				int v, w;
				g.GetEndpoints(e,v,w);
				int bv = voronoi.GetBase(v);
				int bw = voronoi.GetBase(w);
				if (bv == bw) continue; //same base

				bv = uf.Find(bv); 
				bw = uf.Find(bw);
				if (bv == bw) continue; //same component

				//found a boundary edge: move it forward in the list
				if (i!=nextpos) elist[nextpos++] = e;
				
				// get length of actual edge
				EdgeCost cost = (pertcost!=NULL) ? pertcost[e] : g.GetCost(e);
                cost += voronoi.GetDistance(v) + voronoi.GetDistance(w);

				//update bv and bw if better
				if (minarc[bv]==-1 || (cost<minvalue[bv])) {minarc[bv]=e; minvalue[bv]=cost;}
				if (minarc[bw]==-1 || (cost<minvalue[bw])) {minarc[bw]=e; minvalue[bw]=cost;}
			}
			ecount = nextpos; //fewer edges for next round

			//join each region to its best neighbor
			int bvisited = 0;
			int b;
			for (b=1; b<=n; b++) {
				if (minarc[b] >= 0) { //best neighbor defined...
					bvisited ++;
					int w = 0;
					g.GetEndpoints(minarc[b], v, w);
					int bv = voronoi.GetBase(v);
					int bw = voronoi.GetBase(w);
					if (uf.Union(bv,bw)) { // if they were not already joined...
						changes = true;
						nregions--;
						solution.Insert(minarc[b]);

						while (v != bv) {
							int f = voronoi.GetParentArc(v);
							if (!solution.Insert(f)) break;
							v = g.GetOther(f, v);
						}
						while (w != bw) {
							int f = voronoi.GetParentArc(w);
							if (!solution.Insert(f)) break;
							w = g.GetOther(f, w);
						}
					}
				}
			}
		}
		if (verbose) fprintf(stderr, "\n"); 
		delete [] minvalue;
		delete [] minarc;
		delete [] elist;
	}


	static void DNHPrim (Graph &g, SteinerSolution &solution, VoronoiData &voronoi, UnionFind &uf) {
		int n = g.VertexCount();
		int m = g.EdgeCount();

		vector <int> new2old (m+1,-1);

		Graph dg;
		dg.SetVertices(n);
		dg.SetEdges(m); //maximum tentative number of edges

		// build appropriate subgraph of the distance network
		int ecount = 0;
		for (int v=1; v<=n; v++) {
			int bv = uf.Find(voronoi.GetBase(v));
			EdgeCost vdist = -1;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (v>=w) continue;
				int bw = uf.Find(voronoi.GetBase(w));
				if (bv == bw) continue;
				dg.MakeTerminal(bv);
				dg.MakeTerminal(bw);
				if (vdist<0) vdist = voronoi.GetDistance(v);
				EdgeCost cost = a->cost + vdist + voronoi.GetDistance(w);
				dg.AddEdge(bv,bw,cost);
				new2old[++ecount] = a->label;
			}
		}
		dg.Commit(); //create actual graph

		// find a solution in the new graph		
		SteinerSolution ds(&dg);
		MSTPrim(dg,ds);

		// transform into solution in the new graph
		solution.Reset();
		int newm = dg.EdgeCount();
		for (int e=1; e<=newm; e++) {
			if (!ds.Contains(e)) continue;
			int orige = new2old[e]; //original boundary edge
			if (orige<1 || orige>m) {fatal ("Edge out of range.\n");}
			int v,w;
			g.GetEndpoints(orige,v,w);
			solution.Insert(orige);

			for (;;) {
				int f = voronoi.GetParentArc(v);
				if (f<1 || !solution.Insert(f)) break;
				v = g.GetOther(f,v);
			}
			for (;;) {
				int f = voronoi.GetParentArc(w);
				if (f<1 || !solution.Insert(f)) break;
				w = g.GetOther(f,w);
			}
		}
	}


        /// <summary>
        /// Find MST of the distance network associated with a given set of bases.
        /// We can use uf to force the algorithm to think of two bases as a single one;
        /// uf will be altered during the algorithm, as different regions are joined.
        /// (Note: the algorithm could be implemented without uf---by coloring---but
        /// it's simpler to just do everything with uf.)
        /// </summary>
        /// <param name="solution">Output: the solution.</param>
        /// <param name="voronoi">Current voronoi diagram.</param>
        /// <param name="uf">List of regions; a voronoi region must be entirely within the same uf region.</param>
        /// <param name="pertcost">Perturbed edge costs (may be null).</param>

	static void BoruvkaGraph(Graph &g, SteinerSolution &solution, VoronoiData &voronoi, UnionFind &uf, EdgeCost *pertcost) {
		int v, n = g.VertexCount();
		int m = g.EdgeCount();
		EdgeCost solvalue = 0; //needed
        const bool verbose = false;

		//remember vertices with no outgoing boundary edges
		bool *hasneighbors = new bool [n+1];

		//count boundary regions in the current diagram
        int nregions = 0;
        for (v=1; v<=n; v++) {
			hasneighbors[v] = true; //as far as we know, everybody has a neighbor
            if (uf.Find(voronoi.GetBase(v)) == v) nregions ++;
		}

		//Boruvka's algorithm
        int *minarc = new int[n+1]; //minimum outgoing edge from the region based in v
        EdgeCost *minvalue = new EdgeCost[n+1]; //value of this neighbor

		int fakenumber = 0;
        int rounds = 0;
        bool changes = true;

		bool *useless = new bool [m+1];
		for (int e=1; e<=m; e++) {useless[e] = false;}

		const bool debug = false;
		while (changes && nregions > 1) {
			rounds ++;
			//if (rounds > 0) break;
			if (debug && rounds>4) break;
			if (verbose) fprintf (stderr, "%d:%d ", rounds, nregions);
			changes = false; //will be true if there is any progress in this round

			//initially, we don't know what are the arcs out of each component
			for (v=1; v<=n; v++) {
				minarc[v] = -1;
				minvalue[v] = 0;
			}

			int skipped = 0;

			//look for cheapest edge out of each component
			for (v=1; v<=n; v++) {
				if (!hasneighbors[v]) {
					skipped ++;
					if (verbose) fprintf (stderr, ".");
					continue; //skip internal vertices
				} else {
					if (verbose) fprintf (stderr, "|");
				}
				int bv = uf.Find(voronoi.GetBase(v)); //component containing v
				EdgeCost vdist = voronoi.GetDistance(v); //distance from v to its base
				int neighcount = 0;

				hasneighbors[v] = false; //unless we find an actual edge, this is false
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;					
					if (w>v) continue; //look at each edge once (COULD SKIP IF GRAPH GUARANTEED TO BE SORTED)
					if (useless[a->label]) {continue;}

					int bw = uf.Find(voronoi.GetBase(w)); //component containing w
					if (bv==bw) {
						useless[a->label] = true;
						continue; //no longer a boundary edge!
					}

					if (debug) {
						fakenumber += bw;
						changes = true;
						hasneighbors[v] = true;
						continue;
					}


					neighcount ++;
					hasneighbors[v] = true; //found a neighbor! (ACTUALLY, THIS SHOULD ONLY BE TRUE IF WE FIND AT LEAST 2, no?)
					//changes = true;
					//continue;

					//compute cost of corresponding boundary path
					int label = a->label;
					EdgeCost cost = (pertcost!=NULL) ? pertcost[label] : a->cost;
                    cost += vdist + voronoi.GetDistance(w);
					if (cost < 0) fatal ("Something very wrong.");

					//update if better for bv
					if (minarc[bv]==-1 || (cost<minvalue[bv])) {minarc[bv]=label; minvalue[bv]=cost;}

					//update if better for bw
					if (minarc[bw]==-1 || (cost<minvalue[bw])) {minarc[bw]=label; minvalue[bw]=cost;}
				}

				// WARNING: THIS LOOKS VERY WRONG; I GUESS BOTH VERTICES COULD HAVE A SINGLE NEIGHBOR, AND NOBODY WILL BE VISITED. 
				// SOME TIEBREAKING IS NEEDED
				const bool AGGRESSIVE_OPTIMIZATION = false; //true;
				if (AGGRESSIVE_OPTIMIZATION) {
					if (neighcount <= 1) hasneighbors[v] = false; //if this vertex has a single neighbor, it will be merged with it in this iteration
				}
			}

			if (debug && fakenumber==314123411) {fprintf (stderr, "Miracle!\n");}

			//join each region to its best neighbor
			int bvisited = 0;
			int b;
			for (b=1; b<=n; b++) {
				if (minarc[b] >= 0) { //best neighbor defined...
					bvisited ++;
					int w = 0;
					g.GetEndpoints(minarc[b], v, w);
					int bv = voronoi.GetBase(v);
					int bw = voronoi.GetBase(w);
					if (uf.Find(bv) != uf.Find(bw)) { //join regions if not already joined [COULD UNITE DIRECTLY AND CHECK BOOL]
						changes = true;
						uf.Union(bv, bw);
						nregions--;
						solution.Insert(minarc[b]);

						// WARNING: THIS IS NOT STOPPING AS SOON AS IT SHOULD
						while (v != bv) {
							int f = voronoi.GetParentArc(v);
							if (!solution.Insert(f)) break;//; // fprintf (stderr, ".");
							v = g.GetOther(f, v);
						}
						while (w != bw) {
							int f = voronoi.GetParentArc(w);
							if (!solution.Insert(f)) break;
							w = g.GetOther(f, w);
						}
					}
				}
			}
			//fprintf (stderr, "s%.3f b%.3f r%d   ", (double)skipped / (double)n, (double)bvisited /(double)n, nregions);
		}
		//fprintf (stderr, "\n");


		if (verbose) fprintf(stderr, "\n"); //Done computing shortest paths; solution value so far is %d.\n", solvalue);
            //Console.Error.WriteLine("rounds:{0}", rounds);
		delete [] minvalue;
		delete [] minarc;
		delete [] hasneighbors;
		delete [] useless;
	}

        




	/// Run DNH (Boruvka implementation) to create a solution from scratch. 
	/// If 'r' is null, uses the original edge weights; otherwise, perturbs
	/// edge weights before running the algorithm.
	/// <param name="solution">The output solution (will be deleted).</param>
	/// <param name="r">Random number generator.</param>
    static void FullBoruvka(Graph &g, SteinerSolution &solution) { //, OptRandom r) {
		int n = g.VertexCount();
        int m = g.EdgeCount();
		UniverseSet baselist(n); // = new UniverseSet(n);
		VoronoiData voronoi(n); // = new VoronoiData(n);
		UnionFind uf(n); // = new UnionFind(n);
		BinaryHeap<EdgeCost> heap(n); // = new BinaryHeap<ArcCost>(n);

		/*
            ArcCost [] pertcost = null;
            if (r != null) {
                pertcost = new ArcCost[m + 1];
                InitPerturbation(pertcost, r);
            }
			*/

		baselist.Reset();
		for (int v = 1; v <= n; v++) {
			if (g.IsTerminal(v)) baselist.Insert(v);
		}

		//fprintf (stderr, "Should be computing Voronoi.\n");
        ComputeVoronoi(g, voronoi, baselist, heap, NULL); //pertcost);


		solution.Reset();
//        return; 
		//fprintf (stderr, "Missing boruvka!\n");
        //Boruvka(g, solution, voronoi, uf, NULL);
        BoruvkaGraph(g, solution, voronoi, uf, NULL);
		//DNHPrim(g,solution,voronoi,uf);
            //Console.Error.WriteLine("t3:{0} ", timer.GetTime());
	}




        /// <summary>
        /// Compute the Voronoi diagram of the current graph, given a set of bases
        /// and maybe the perturbed cost of the edges.
        /// </summary>
        /// <param name="voronoi">Output: description of the Voronoi diagram</param>
        /// <param name="baselist">List of bases.</param>
        /// <param name="heap">Preallocated heap to be used in the computation (will be reset)</param>
        /// <param name="pertcost">Edge costs (use original costs if null).</param>
	static void ComputeVoronoi(Graph &g, VoronoiData &voronoi, UniverseSet &baselist, BinaryHeap<EdgeCost> &heap, EdgeCost *pertcost) {
		const bool GLOBAL_USE_VORONOI_TIE_BREAKER = false; // NOT CLEAR WHERE THIS THING IS SUPPOSED TO BE DEFINED
			
		const bool verbose = false;
		voronoi.Reset();
        int nbases = 0;
        heap.Reset();

		// initialize with all bases
		int p, pend;
		for (baselist.GetBounds(p,pend); p<pend; p++) {
			int b = baselist.PickPos(p);
			nbases ++;
            voronoi.MakeBase(b);
            heap.Insert(b, 0);
		}

        if (verbose) fprintf (stderr, "%d vertices marked as bases.\n", nbases);

        int count = 0;
        //WARNING: RANDOMIZING THE CHOICE SEEMS TO BE A GOOD IDEA
        bool randomize = false;
        bool PREFER_TERMINALS = true;
        bool USE_TIEBREAKERS = GLOBAL_USE_VORONOI_TIE_BREAKER && (randomize || PREFER_TERMINALS);

        //perform multisource Dijkstra
        while (!heap.IsEmpty()) {
			unsigned int v;
            EdgeCost dist;
            heap.RemoveFirst(v, dist);
			count++;
			if (verbose) fprintf (stderr, "%d ", dist);

            //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //g.GetArcHead(pa);
				EdgeCost newdist = dist;

				if (pertcost == NULL) newdist += a->cost; //g.GetArcCost(pa);
				else newdist += pertcost[a->label];
                    
				bool improve = false;

				if (voronoi.GetBase(w) == 0) improve = true;
				else if (newdist <= voronoi.GetDistance(w)) improve = true; //using leq here to prefer shorter edges...
				else if (USE_TIEBREAKERS && newdist == voronoi.GetDistance(w)) {
					if (randomize) {
						fatal ("randomization not implemented!\n"); //(stderr, "NOT IMPLEMENTED!\n
                        //improve = (random.GetInteger(0, 1) == 0); //(arc.cost < g.GetCost(voronoi.GetParentArc(w)));
					} else if (PREFER_TERMINALS) {
						improve  = g.IsTerminal(voronoi.GetBase(v));
					}
                    //improve = (arc.cost < g.GetCost(voronoi.GetParentArc(w)));
				}

				if (improve) { //make w a tentative child of v
					heap.Insert(w, newdist);
					voronoi.Update(w, voronoi.GetBase(v), a->label, newdist);
				}
			}
		}
	}

	// add all vertices in the current solution to solnodes
	// (vertices with incident edges)
	static void MarkSolutionNodes(Graph &g, SteinerSolution &solution, UniverseSet &solnodes) {
		int n = g.VertexCount();
		for (int v=1; v<=n; v++) {
			if (solution.GetDegree(v)>0) solnodes.Insert(v); 
		}
	}

	/// <summary>
	/// Add to 'keynodes' all key vertices in the solution; if 'terminals'=true, adds
	/// the terminals as well. Takes O(n) time.
	/// </summary>
	/// <param name="solution">Input: the solution.</param>
	/// <param name="keynodes">Output with every keynode marked.</param>
    static void MarkKeyNodes(Graph &g, SteinerSolution &solution, UniverseSet &keynodes, bool terminals) {
		int n = g.VertexCount();

        for (int v=1; v<=n; v++) {
			int d = solution.GetDegree(v); //degree in the solution
            if (d==0) continue; //not there: fine
            if (g.IsTerminal(v)) { //add terminals if 'terminals' is true
				if (terminals) keynodes.Insert(v);
			} else { //add nonterminals of degree at least three
				if (d>2) keynodes.Insert(v);
                else if (d==1) {
					fprintf (stderr, "Solution has degree-one nonterminal!");   
				}
			}
		}
	}





	/// Modifies a solution by computing the MST of the subgraph induced by its vertices 
    /// and removing all vertices of degree one. Changes the solution.
	/// <param name="svertices">Vertices in the solution (doesn't change).</param>
    /// <param name="solution">Edges in the solution (may change).</param>
	static bool MSTPrune(Graph &g, SteinerSolution &solution) {
            //return 0;
            EdgeCost original = solution.GetCost();
            int n = g.VertexCount();
            UniverseSet svertices(n); 
            MarkSolutionNodes(g, solution, svertices);
            MST(g,solution,svertices);
            Prune(g,solution);
			//fprintf (stderr, "Solution costs %d, gain is %d.\n", solution.GetCost(), original - solution.GetCost());
            //Prune(g,solution);
			//fprintf (stderr, "g%d ", original - solution.GetCost());
            return (solution.GetCost() < original) ? 1 : 0;
        }

	/// Iteratively removes all degree-one vertices in the solution.
	/// Takes O(1) if there are no such vertices. If there are, takes
	/// O(n) + degree of all vertices removed.
	/// <param name="solution">Original solution (will be modified).</param>
	static void Prune(Graph &g, SteinerSolution &solution) {
		if (solution.LeafCount()==0) return; //WARNING: THIS SHOULD BE THERE!

        int n = g.VertexCount();
        for (int v=1; v<=n; v++) {
			int t = v;
            while (solution.GetDegree(t)==1 && !g.IsTerminal(t)) {
				//find the unique incident solution edge
				SPGArc *a, *end;
				for (g.GetBounds(t,a,end); a<end; a++) {
					int alabel = a->label;
					if (!solution.Contains(alabel)) continue; 
					solution.Remove(alabel);
                    t = g.GetOther(alabel,t); //process the other endpoint next
                    break;
				}
			}
		}
	}

	static int VeryOldPickTerminal(Graph &g, RFWLocalRandom &random) {
		int t = 0;
		int count = 0;
		int n = g.VertexCount();
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) {
				if (t==0 || g.GetDegree(v) < g.GetDegree(t)) {
					t = v;
				}
			}
		}
		return t;
	}

	// this is old; the terminal is not random
	static int WrongPickRandomTerminal(Graph &g) {
		//fprintf (stderr, "r");
		int n = g.VertexCount();
		for (int v=1; v<=n; v++) if (g.IsTerminal(v)) return v;
		fatal ("could not find terminal");
		return 0;
	}

	static int PickRandomTerminal(Graph &g) {
		//fprintf (stderr, "r");
		int n = g.VertexCount();
		int count = 0;
		int target = RFWRandom::getInteger(1,g.TerminalCount());
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) {
				if (++count == target) return v;
			}
		}
		fatal ("could not find terminal");
		return 0;
	}

	static int PickRandomTerminal(Graph &g, RFWLocalRandom &random) {
		static bool first = true;
		if (first) {
			fprintf (stderr, "PICKRANDOMTERMINAL IS NOT PROPERLY SET.\n");	
			first = false;
		}
		int n = g.VertexCount();
		int count = 0;
		int target = random.GetInteger(1,g.TerminalCount());
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) {
				if (++count == target) return v;
			}
		}
		fatal ("could not find terminal");
		return 0;
	}





	/// <summary>
        /// Compute the MST of the subgraph induced by svertices
        /// (assumed to contain all terminals---DO I NEED THIS?)
        /// </summary>
        /// <param name="svertices">list of vertices to be spanned</param>
        /// <param name="solution">solution in which the MST will be stored</param>
        /// <returns>Number of vertices scanned.</returns>
	static int MST (Graph &g, SteinerSolution &solution, UniverseSet &svertices) {
            bool verbose = false;
            int n = g.VertexCount();
            BinaryHeap<EdgeCost> heap(n);
			vector<int> parc (n+1);
            int r = PickRandomTerminal(g);
            if (!svertices.Contains(r)) fatal ("Terminal does not appear to belong to the solution.");
            parc[r] = 0;
            int nscanned = 0;

            //run Prim's algorithm
			solution.Reset();
            heap.Insert(r, 0);
            while (!heap.IsEmpty()) {
                unsigned int v;
                EdgeCost acost;
                heap.RemoveFirst(v,acost);
                if (v!=r) solution.Insert(parc[v]); //add edge (p(v),v) to solution

                //scan vertices
                nscanned ++;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
                    int w = a->head; 
                    if (!svertices.Contains(w)) continue; //we only care about svertex
                    if (solution.GetDegree(w) > 0) continue; //vertex already in the new tree
                    if (heap.Insert(w, a->cost)) parc[w] = a->label;
                }
            }
			//if (solution.Count() < g.TerminalCount() - 1) {fatal ("solution does not have enough vertices");}
            return nscanned;
        }



	template <class GRAPH> static void PruneNeighbors(GRAPH &g, int v, vector<bool> &arcs, vector<int> &degree, RFWStack<int> &candidates) {
		SPGArc *a, *end;
		for (g.GetBounds(v,a,end); a<end; a++) {
			if (!arcs[a->label]) continue;
			arcs[a->label] = false;
			//g.HideEdge(a->label); // this only works for dynamic graph
			int w = a->head;
			if (--degree[w] <= 2) {
				if (!g.IsTerminal(w)) candidates.push(w);
			}
		}
		degree[v] = 0;
	}

	static void DynamicToStatic (DynamicGraph &oldg, Graph &newg) {
		int newn = 0;
		int newm = 0;
		int oldn = oldg.VertexCount();
		const bool verbose = false;

		vector<int> old2new(oldn+1, -1); 
		if (verbose) fprintf (stderr, "Creating static graph (from %d)...\n", oldn); fflush(stderr);
		for (int v=1; v<=oldn; v++) {
			if (!oldg.IsTerminal(v) && (oldg.GetDegree(v) == 0)) continue;
			old2new[v] = ++newn;
		}
		if (verbose) fprintf (stderr, "New graph should have %d vertices.\n", newn);

		newg.SetVertices(newn);
		newg.coord.SetMaxId(newn);
		for (int v=1; v<=oldn; v++) {
			if (old2new[v]<0) continue;
			double x, y;
			oldg.coord.GetCoordinates(v,x,y);
			newg.coord.SetCoordinates(old2new[v],x,y);
			if (oldg.IsTerminal(v)) newg.MakeTerminal(old2new[v]);
			SPGArc *a, *end;
			for (oldg.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (w <= v) continue;
				newg.AddEdge(old2new[v],old2new[w],a->cost);
			}
		}

		newg.SetFixedCost(oldg.GetFixedCost());
		newg.Commit();

		if (verbose) fprintf (stderr, "Created graph with %d vertices, %d edges, and %d terminals (and fixed cost %.3f).\n", newg.VertexCount(), newg.EdgeCount(), newg.TerminalCount(), newg.GetFixedCost());
		
	}

	static void CopyFromStaticGraph (DynamicGraph &dg, Graph &sg, vector<bool> &arcs) {
		const bool verbose = false;
		if (verbose) fprintf (stderr, "Copying from static graph.\n");

		dg.SetVertices(sg.VertexCount());
		dg.SetEdges(sg.EdgeCount());
		for (int e=1; e<=sg.EdgeCount(); e++) {
			int v, w;
			sg.GetEndpoints(e,v,w);
			dg.BatchAddEdge(v,w,sg.GetCost(e));
		}
		for (int v=1; v<=sg.VertexCount(); v++) {
			if (sg.IsTerminal(v)) dg.MakeTerminal(v);
		}

		dg.Commit();
		dg.SetFixedCost(sg.GetFixedCost());
		dg.coord = sg.coord;

		// WARNING: THIS SHOULD BE BATCHED
		int hidden = 0;
		if (verbose) fprintf (stderr, "Hiding edges that have already been deleted!\n");
		for (int e=1; e<=sg.EdgeCount(); e++) {
			if (!arcs[e]) dg.HideEdge(e);
		}
		if (verbose) fprintf (stderr, "done hiding.\n");
		if (hidden>0) fprintf (stderr, "h%d ", hidden);
	}

	static void OneTwoPrep (Graph &sg, vector<bool> &arcs) {
		DynamicGraph g;
		CopyFromStaticGraph (g, sg, arcs);
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int v, w;

		const int VERBOSE_LEVEL = 1;
		bool DO_NOTHING = false;
		const bool BOTTLE_TEST = true;
		const bool verbose = false;

		int MAX_DEGREE = 10;


		if (VERBOSE_LEVEL>=1) fprintf (stderr, "OneTwoPrep(%d,%d) ", n, m);
		//fprintf (stderr, "Vertex %d has degree %d\n", DEBUG_VERTEX, g.GetDegree(DEBUG_VERTEX));
		EdgeCost infinity = 10e32;
		vector<EdgeCost> distance (n+1, infinity);

		//vector<int> edgelist; // stack of edges to process
		vector<int> onelist;
		vector<int> twolist;
		vector<int> bottlelist; 
		// first copy all edges incident vertices of degree 1 and 2
		for (int e=1; e<=m; e++) {
			if (!arcs[e]) continue;
			g.GetEndpoints(e,v,w);
			if (g.GetDegree(v)==1 || g.GetDegree(w)==1) {onelist.push_back(e); continue;}
			if (g.GetDegree(v)<=2 || g.GetDegree(w)<=2) {twolist.push_back(e); continue;}
			if (BOTTLE_TEST) bottlelist.push_back(e);
		}

		while ((onelist.size()>0 || twolist.size()>0 || bottlelist.size()>0) && !DO_NOTHING) {
			int e;
			if (onelist.size() > 0) {
				e = onelist.back();
				onelist.pop_back();
			} else if (twolist.size() > 0) {
				e = twolist.back();
				twolist.pop_back();
			} else if (bottlelist.size() > 0) {
				if (!BOTTLE_TEST) fatal ("bottleneck processing should not be active");
				//fprintf (stderr, "B");
				e = bottlelist.back();
				bottlelist.pop_back();
			}
			if (!arcs[e]) continue;

			g.GetEndpoints(e,v,w);
			if (v==w) {
				if (VERBOSE_LEVEL>=3) fprintf (stderr, "FOUND POTENTIAL SELF-LOOP!");
				if (g.IsTerminal(v)) {
					fprintf (stderr, "Apparently %d is a selfloop.\n", v);
					fatal ("And it's a terminal.\n");
				}
				arcs[e] = false;
				g.HideEdge(e);
				continue;
			}

			if ((g.IsTerminal(v) && g.GetDegree(v)>1) || g.GetDegree(v) > g.GetDegree(w)) std::swap(v,w);

			// v has the lowest degree / is the non-terminal


			/*
			if (g.IsTerminal(v)) {
				if (g.GetDegree(v) == 1) {
					arcs[e] = false;
					g.HideEdge(e);
					g.GetFixedCost(g.GetCost(e));
					if (!g.IsTerminal(
				}
				continue;
			}*/
			//if (g.IsTerminal(v)) continue;

			//fprintf (stderr, "[%d,%d] ", g.GetDegree(v), g.GetDegree(w));
			// v is the vertex we will eliminate
			if (g.GetDegree(v)==1) {
				if (verbose) fprintf (stderr, "Removing (%d,%d) with degrees (%d,%d) and length %.0f.\n", v, w, g.GetDegree(v), g.GetDegree(w), g.GetCost(e));
				arcs[e] = false;
				g.HideEdge(e);

				
				if (g.IsTerminal(v)) {
					g.IncFixedCost(g.GetCost(e));
					g.UnmakeTerminal(v);
					if (!g.IsTerminal(w)) g.MakeTerminal(w);
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "Removed %d, made %d terminal (%d).\n", v, w, g.GetDegree(w));
				}


				if (g.GetDegree(w)==2) {
					twolist.push_back(g.GetFirstEdgeLabel(w));
				} else if (g.GetDegree(w)==1) {
					onelist.push_back(g.GetFirstEdgeLabel(w));
				}

				if (g.GetDegree(w)==0 && g.IsTerminal(w)) {
					int vleft = 0;
					for (int y=1; y<=g.VertexCount(); y++) {
						if (g.GetDegree(y)>0) vleft ++;
					}
					fprintf (stderr, "neighbor (%d) of degree-one vertex (%d) was a degree-one terminal; vertices left=%d.\n", w, v, vleft);
					if (vleft > 0) fatal ("that shouldn't happen");
					else break;
				}

				// should already be in bottlelist, I guess
				continue;
			}

			if (g.IsTerminal(v)) continue;

			if (g.GetDegree(v)==2) {
				int neighbors[2];
				int edges[2];
				int ncount = 0;
				EdgeCost length = 0;

				// get list of neighbors
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (!arcs[a->label]) fatal ("should not see deleted vertices");
					if (ncount == 2) {
						fprintf (stderr, "Vertex %d already has %d %d as neighbors.\n", v, neighbors[0], neighbors[1]);
						fatal ("too many neighbors");
					}
					neighbors[ncount] = a->head;
					edges[ncount] = a->label;
					length += a->cost;
					arcs[a->label] = false;
					ncount ++;
				}
				if (ncount != 2) fatal ("not enough neighbors");


				g.HideEdge(edges[0]);
				g.HideEdge(edges[1]);


				if (neighbors[0] == neighbors[1]) {
					if (verbose) fprintf (stderr, "Removing potential self-loop!\n");
					int u = neighbors[0];
					if (g.GetDegree(u)==1) onelist.push_back(g.GetFirstEdgeLabel(u));
					else if (g.GetDegree(u)==2) twolist.push_back(g.GetFirstEdgeLabel(u));
					else if (BOTTLE_TEST) bottlelist.push_back(g.GetFirstEdgeLabel(u));
				} else {
					if (verbose) fprintf (stderr, "Removed <%d-%d-%d> with length %.8f + %.8f.\n", neighbors[0], v, neighbors[1], g.GetCost(edges[0]), g.GetCost(edges[1]));
					if (g.IsTerminal(v)) fatal ("About to remove a terminal!\n");

					int eid = g.InsertEdge(neighbors[0], neighbors[1], length);

					if (eid != (int)arcs.size()) {
						fprintf (stderr, "New id is %d, current size is %d.\n", eid, arcs.size());
						fatal ("invalid ID for new edges");
					}
					arcs.push_back(true);

					int x = neighbors[0];
					int y = neighbors[1];

					//fprintf (stderr, "
					if (g.GetDegree(x)==1 && g.GetDegree(y)==1) onelist.push_back(eid);
					else if (g.GetDegree(x)==2 && g.GetDegree(y)==2) twolist.push_back(eid);
					else if (BOTTLE_TEST) bottlelist.push_back(eid);
				}
				continue;

			}

			// OK, we'll try to remove the edge by bottleneck distance
			if (!(g.IsTerminal(w) && g.GetDegree(w)==1)) {
				//fprintf (stderr, "B");

				if (g.GetDegree(v)+g.GetDegree(w) >= 2*MAX_DEGREE) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "!");
					continue;

				}

				EdgeCost length = g.GetCost(e);

				SPGArc *a, *end;
				bool prune = false;

				// find distances from v
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label] && a->label!=e) {distance[a->head] = a->cost;}
				}

				if (distance[w] <= length) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "Prunning at equality!\n");
					prune = true;
				} else {
					// try to find matches
					for (g.GetBounds(w,a,end); a<end; a++) {
						if (arcs[a->label] && a->label!=e) {
							if (distance[a->head] + a->cost <= length) {
								prune = true; break;
							}
						}
					}
				}

				// reset distances
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label]) {
						if (a->head == w) {
							if (a->cost > length) {
								if (BOTTLE_TEST) bottlelist.push_back(a->label); //CAREFUL: MAY ADD SEVERAL TIMES
							}
						}
						distance[a->head] = infinity;
					}
				}

				if (prune) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "REMOVING EDGE (%d,%d)\n", v, w);
					arcs[e] = false;
					g.HideEdge(e);
					if (g.GetDegree(v) == 2) twolist.push_back(g.GetFirstEdgeLabel(v));
					else if (g.GetDegree(v) == 1) onelist.push_back(g.GetFirstEdgeLabel(v));
					if (g.GetDegree(w) == 2) twolist.push_back(g.GetFirstEdgeLabel(w));
					else if (g.GetDegree(w) == 1) onelist.push_back(g.GetFirstEdgeLabel(w));
				}

			}

			//fprintf (stderr, "Should not have gotten here.");
		}
		//fprintf (stderr, "Done with OneTwo.\n");
		bool OUTPUT_STUFF = false;
		if (OUTPUT_STUFF) {
			int kept = 0;
			vector<bool> vkeep(n+1, true);
			for (int v=1; v<=n; v++) {
				if (g.GetDegree(v) == 0) vkeep[v] = false;
				else kept ++;
			}
			GraphDrawer::DrawSubgraph("test.gv", g, arcs, vkeep);
			//OutputSTP("dtest.stp", g, arcs);		
			OutputCompactSTP("dtest.stp", g, arcs);		
			fprintf (stderr, "Kept %d/%d vertices.\n", kept, n);
		}

		Graph newg;
		DynamicToStatic(g, newg);
		sg = newg;
	}

	// recursively remove corners
	static void PrepCorners(Graph &sg, vector<bool> &arcs, vector<int> &degree) {
		DynamicGraph g;
		CopyFromStaticGraph (g, sg, arcs);
		
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int v, w;



		fprintf (stderr, "Preprocessing corners (%d,%d)!\n", n, m);

		EdgeCost infinity = 10e32;
		vector<EdgeCost> distance (n+1, infinity);


		//reset degrees
		for (int v=1; v<=n; v++) {
			degree[v] = 0;
		}

		//compute actual degrees
		for (int e=1; e<=m; e++) {
			if (arcs[e]) {
				g.GetEndpoints(e,v,w);
				fprintf (stderr, "(%d,%d)", v, w);
				degree[v]++;
				degree[w]++;
			}
		}

		RFWStack<int> candidates(n);
		for (int v=1; v<=n; v++) {
			if (degree[v] != g.GetDegree(v)) fatal ("invalid degree in the graph");
			if (degree[v] <= 2 && !g.IsTerminal(v)) candidates.push(v);
		}

		while (!candidates.isEmpty()) {
			v = candidates.pop();
			fprintf (stderr, "Checking vertex %d with degree %d.\n", v, degree[v]);

			if (degree[v] == 1) {
				PruneNeighbors(g, v, arcs, degree, candidates);
				continue;
			}

			if (degree[v] == 2) {
				fprintf (stderr, "Vertex %d has degree %d.\n", v, degree[v]);
				int neighbor[2];
				EdgeCost cost[2];
				int ncount =0;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label]) {
						cost[ncount] = a->cost;
						neighbor[ncount] = a->head;
						if (++ncount == 2) break;
					}
				}
				fprintf (stderr, "Vertex %d has neighbors %d and %d.\n", v, neighbor[0], neighbor[1]);
				EdgeCost length = cost[0] + cost[1];
				bool prune = false;
				if (neighbor[0] == neighbor[1]) {
					prune = true;
				} else {
					for (g.GetBounds(neighbor[0],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {distance[a->head] = a->cost;}
					}

					for (g.GetBounds(neighbor[1],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {if (distance[a->head] + a->cost <= length) prune = true;}
					}

					for (g.GetBounds(neighbor[0],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {distance[a->head] = infinity;}
					}
				}

				if (prune) {
					PruneNeighbors(g, v, arcs, degree, candidates);
				}
			}
		}
	}

	//-----------------------------------------------
	// Recursively remove all vertices of degree one 
	//-----------------------------------------------
	static void PrepDegreeOne(Graph &g, vector<bool> &arcs, vector<int> &degree) {
		int n = g.VertexCount();

		// warning: could loop over edges instead
		int v;
		for (v = 1; v <= n; v++) {
			degree[v] = 0;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				if (arcs[a->label]) degree[v]++;
			}
		}

		int removed = 0;
		int terminals = 0;
		int i;
		for (i=1; i<=n; i++) {
			v = i;
			while (!g.IsTerminal(v) && degree[v]==1) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int e = a->label;
					if (arcs[e]) {
						removed ++;
						arcs[e] = false;
						degree[v] --;
						int w = a->head;
						degree[w] --;
                        if (degree[w] == 1) {v = w;}
						break;
					}
				}
				if (a>=end) fatal ("did not find parent edge");
			}
            
			if (degree[v] == 1) {
				terminals++;
			}
		}
		fprintf (stderr, "Removed %d vertices, and there are still %d terminals of degree 1.", removed, terminals);
	}


	/*-------------------------------------------------------------
	 | Create a new graph containing only the edges marked as kept
	 *------------------------------------------------------------*/
	static void CommitSubgraph (Graph &newg, Graph &g, vector<bool> &edgekeep) {
		int oldn = g.VertexCount();
		int oldm = g.EdgeCount();
		int newm = 0;
		int v, w, e;
		const bool verbose = false;

		if (verbose) fprintf (stderr, "Creating subgraph after preprocessing.\n");

		vector<int> degree(oldn+1, 0);
		
		// figure out the degree
		for (e=1; e<=oldm; e++) {
			if (edgekeep[e]) {
				g.GetEndpoints(e, v, w);
				degree[v] ++;
				degree[w] ++;
				newm ++;
			}
		}

		int newcount = 0;
		vector<int> old2new(oldn+1, -1); 
		int onecount = 0;
		int twocount = 0;
		for (v=1; v<=oldn; v++) {
			if (degree[v] == 0) continue;
			if (degree[v] == 1 && !g.IsTerminal(v)) {onecount ++;} // continue;}
			if (degree[v] == 2 && !g.IsTerminal(v)) {twocount ++;}
			old2new[v] = ++newcount;
		}

		if (verbose) fprintf (stderr, "Skipped %d vertices of degree 1 and counted %d of degree 2; graph has %d vertices.\n", onecount, twocount, newcount);
		//fprintf (stderr, "Graph has %d vertices.\n", newcount);

		// FINALLY CREATE NEW GRAPH
		newg.SetVertices(newcount);
		newg.coord.SetMaxId(newcount);
		newg.SetEdges(newm);
		newg.SetFixedCost(g.GetFixedCost());
		for (e=1; e<=oldm; e++) {
			if (edgekeep[e]) {
				g.GetEndpoints(e, v, w);
				newg.AddEdge(old2new[v],old2new[w],g.GetCost(e));
			}
		}
		for (v=1; v<=oldn; v++) {
			if (old2new[v] < 0) continue;
			if (g.IsTerminal(v)) newg.MakeTerminal(old2new[v]);
			double x, y;
			g.coord.GetCoordinates(v,x,y);
			newg.coord.SetCoordinates(old2new[v], x, y);
		}
		newg.Commit();

		if (verbose) fprintf (stderr, "Created actual graph with %d vertices and %d edges.\n", newg.VertexCount(), newg.EdgeCount());

	}

	static void RunPreprocessing (Graph &sg, const bool runOneTwo) {
		
		int origm = sg.EdgeCount();
		double threshold = 0.01;
		int rounds;
		for (rounds=0; rounds<999; rounds++) {
			int oldedges = sg.EdgeCount();
			fprintf (stderr, "Running iteration %d (fixed:%.3f).", rounds, sg.GetFixedCost());
			PrepBottleneck(sg);
			//fprintf (stderr, "(fixed after prepbottle:%.3f).", sg.GetFixedCost());
			//GraphDrawer::DrawGraph("graph1.gv", sg);
			int m = sg.EdgeCount();		

			if (runOneTwo) {
				vector<bool> keep(m + 1, true);
				OneTwoPrep(sg, keep);
			}
			//fprintf (stderr, "<<<< F%.3f >>>>\n", sg.GetFixedCost());

			int newedges = sg.EdgeCount();
			if (newedges == 0) {
				fprintf (stderr, "Instance solved to optimality!\n");
				break;
			}
			double reduction = (double)(oldedges - newedges) / (double)oldedges; 
			fprintf (stderr, "R%d:%.3f ", rounds, reduction);
			if (reduction < threshold) {
				fprintf (stderr, "!\n"); //"Reduction was only %.3f; stopping.\n", reduction);				
				break;
			}
		}
		//GraphDrawer::DrawGraph("graph2.gv", sg);

		fprintf (stderr, "Done preprocessing.\n");
		fflush(stderr);
		fprintf (stdout, "preprounds %d\n", rounds+1);
		fprintf (stdout, "prepthreshold %.3f\n", threshold);
		
		double reduction = (origm - sg.EdgeCount()) / (double)origm;
		fprintf (stdout, "prepreduction %.6f\n", reduction); 
		fprintf (stdout, "prepreductionpct %.6f\n", 100.0 * reduction); 
		fprintf (stdout, "prepsolved %d\n", sg.EdgeCount()==0);
		fprintf (stdout, "prepvertices %d\n", sg.VertexCount());
		fprintf (stdout, "prepedges %d\n", sg.EdgeCount());
		fprintf (stdout, "prepterminals %d\n", sg.TerminalCount());

		/*
		DynamicGraph g;
		vector<bool> arcs (sg.EdgeCount()+1, true);
		CopyFromStaticGraph (g, sg, arcs);
		PrepBottleneck(sg);*/
	}



    // Preprocess the current instance using bottleneck information: remove edge (v,w) if
	// the bottleneck distance between the endpoints is less than or equal to its length
	static void PrepBottleneck(Graph &g) {
		int n = g.VertexCount();
		int m = g.EdgeCount();
        int v;
        fprintf (stderr, "Preprocessing instance... ");

		const bool verbose = true;
		SteinerSolution solution(&g);
		UniverseSet baselist(n);
		//vector<int> degree(n+1);

		for (v=1; v<=n; v++) {
			if (g.IsTerminal(v)) baselist.Insert(v);
			//degree[v] = g.GetDegree(v);
		}

		BinaryHeap <EdgeCost> binheap(n); // = new BinaryHeap<ArcCost> (n);
		VoronoiData vordata(n);// = new VoronoiData(n);
		UnionFind uf(n); // = new UnionFind(n);
		ComputeVoronoi(g,vordata,baselist,binheap,NULL); 
		BoruvkaGraph(g,solution,vordata,uf,NULL);

		if (verbose) fprintf (stderr, "Found solution worth %.3f.\n", solution.GetCost());

		// find maximum bottleneck length within the graph
		EdgeCost maxbottle = 0;
		for (int e=1; e<=m; e++) {
			if (!solution.Contains(e)) continue;
			int v, w;
			g.GetEndpoints(e,v,w);
			if (vordata.GetBase(v)!=vordata.GetBase(w)) { //found bridge
				EdgeCost pathcost = vordata.GetDistance(v) + vordata.GetDistance(w) + g.GetCost(e);
				if (pathcost > maxbottle) maxbottle = pathcost;
			}
		}
		
		// remove stupid edges
		int neliminated = 0;
		vector<int> elist;
		vector<bool> keep (m+1,true);
		vector<EdgeCost> priority (m+1,-1); 
		for (int e=1; e<=m; e++) {
			if (solution.Contains(e)) {
				int v,w;
				g.GetEndpoints(e,v,w);
				int bv = vordata.GetBase(v);
				int bw = vordata.GetBase(w);
				if (bv!=bw) {
					priority[e] = g.GetCost(e) + vordata.GetDistance(v) + vordata.GetDistance(w); 
				} 
			} else {
				EdgeCost cost = g.GetCost(e);
				priority[e] = cost; //min(min(dv,dw),cost); //to make sure it appears after all edges
				if (cost < maxbottle) continue;
				EdgeCost dv = vordata.GetDistance(g.GetFirstEndpoint(e));
				if (cost < dv) continue;
				EdgeCost dw = vordata.GetDistance(g.GetSecondEndpoint(e));
				if (cost < dw) continue;
				neliminated ++;
				keep[e] = false;
			}
		}
		for (int e=1; e<=m; e++) {
			if (keep[e] && priority[e]>=0) {
				elist.push_back(e);
				//fprintf (stderr, " %d", g.GetCost(e));
			}
		}
		if (verbose) fprintf (stderr, "\n");
		int kept = (int)elist.size();

		//fprintf (stderr, "Maximum bottleneck length is %.3f.\n", maxbottle);
		if (verbose) fprintf (stderr, "Maxbotttle=%.3f eliminated %d/%d edges, keeping %d(%d).\n", maxbottle, neliminated, m, m - neliminated, kept);
		//fflush(stderr);
		//const bool verbose = false;

		bool FANCY_BOTTLENECK = true;
		if (FANCY_BOTTLENECK) {
			if (verbose) fprintf(stderr, "elist.size=%d, kept=%d, priority.size=%d, edge.size=%d\n", elist.size(), kept, priority.size(), solution.EdgeCapacity());
			
			// Sort all relevant edges by their bottleneck relevance. In case of ties, solution edges are inserted first.
			sort(elist.begin(), elist.begin() + kept, [&](int x, int y) {
				if (verbose) { 
					assert(x < priority.size()); 
					assert(y < priority.size()); 
					assert(x < solution.EdgeCapacity()); 
				} 
				return (priority[x] < priority[y] || (priority[x] == priority[y] && solution.Contains(x) && !solution.Contains(y)));
			});
			
			if (verbose) fprintf(stderr, "Sorted.\n");

			EdgeCost curbottle = 0;
			UnionFind newuf(n);
			for (int i=0; i<kept; i++) {
				//fprintf (stderr, "-");
				int e = elist[i];
				int v, w;
				g.GetEndpoints(e,v,w);
				bool localdebug = false;

				/*
				if (v==107 && w==115) {
					fprintf (stderr, "Processing!\n");
					localdebug = true;
				}*/

				if (solution.Contains(e)) {
					if (localdebug) fprintf (stderr, "Edge %d belongs to the solution.\n", e);
					int bv = vordata.GetBase(v);
					int bw = vordata.GetBase(w);
					newuf.Union(bv,bw);
					if (verbose && (curbottle!=priority[e])) fprintf (stderr, "   B%d ", priority[e]);
					curbottle = priority[e];
					//if (curbottle==2) exit(-1);
				} else {
					//this is an edge we have to check
					EdgeCost cost = priority[e];
					if (cost!=g.GetCost(e)) fprintf (stderr, ".");
					int v, w;
					g.GetEndpoints(e,v,w);
					if (vordata.GetParentArc(v)==e || vordata.GetParentArc(w)==e) {
						if (localdebug) fprintf (stderr, "Edge %d is a parent edge!\n", e);
						if (verbose) fprintf (stderr, ".%.0f", g.GetCost(e));
						keep[e] = true;
					} else {
						int bv = vordata.GetBase(v);
						int bw = vordata.GetBase(w);
						if (localdebug) fprintf (stderr, "Edge %d is a candidate!\n", e);
						if (newuf.Find(bv) == newuf.Find(bw)) {
							if (localdebug) fprintf (stderr, "Edge %d has both endpoints in the same component!\n", e);
							EdgeCost dv = vordata.GetDistance(v);
							EdgeCost dw = vordata.GetDistance(w);

							if (cost < curbottle) fatal ("invalid order");
							EdgeCost localbottle = max(dv,dw);
							if (cost >= localbottle) {
								
								if (DEBUG_VERTEX>=0) {
									if (v==DEBUG_VERTEX || w==DEBUG_VERTEX) {
										fprintf (stderr, "Removing edge (%d,%d).\n", v, w);
										fprintf (stderr, "dv=%.3f dw=%.3f cost=%.3f (v=%d w=%d bv=%d bw=%d) (e=%d pw=%d pv=%d) lb=%.3f cb=%.3f\n", dv, dw, cost, v, w, bv, bw, e, vordata.GetParentArc(v), vordata.GetParentArc(w), localbottle, curbottle);
									}
								}
								
								keep[e] = false;
							}
						}
					}
					//fprintf (stderr, "%d", (int)keep[e]);
				}
			}
		}

		int kcount = 0;
		for (int e=1; e<=m; e++) {
			if (keep[e]) {
				kcount++;
			}
		}
		if (verbose) fprintf (stderr, "\nKept %d edges eventually.\n", kcount); fflush(stderr);

		const bool REMOVE_NON_SP_EDGES = true;
		if (REMOVE_NON_SP_EDGES) {
			// remove edges with 
			int removed = 0;
			vector<int> last(n+1,-1);
			vector<EdgeCost> value(n+1,0);
			for (int e=1; e<=m; e++) {
				if (!keep[e]) continue;
				EdgeCost ecost = g.GetCost(e);
				int v,w;
				g.GetEndpoints(e,v,w);
				for (int s=0; s<=1; s++) {
					int x = (s==0) ? v : w; 
					SPGArc *a, *end;
					for (g.GetBounds(x,a,end); a<end; a++) {
						if (!keep[a->label]) continue;
						if (e==a->label) continue;
						int y = a->head;
						EdgeCost cost = a->cost;
						if (s==0) {
							if (last[y]!=x || value[y]>cost) {
								last[y] = x;
								value[y] = cost;
							}
						} else {
							if ((last[y]==v) && (cost+value[y] <= ecost)) {
								if (verbose) fprintf (stderr, "Removing %d-%d-%d\n", v, y, w);
								keep[e] = false;
							}
						}
					}
				}
				if (!keep[e]) removed++;
				//fprintf (stderr, "%d ", (int)keep[e]);
			}
			fprintf (stderr, "Removed %d additional edges, left with %d.\n", removed, kcount - removed);
		}

		fprintf (stderr, "Removed non-sp edges.\n");
		fflush(stderr);

		bool SHOW_KEPT = false;
		if (SHOW_KEPT) {
			fprintf (stderr, "\nKept edges: ");
			for (int e=1; e<=m; e++) {
				if (keep[e]) fprintf (stderr, "%d ", g.GetCost(e));
			}
			fprintf (stderr, "\n");
		}

		vector<int> degree(n+1,0);
		int maxdegree = 0;
		for (int e=1; e<=m; e++) {
			if (keep[e]) {
				int dv = ++degree[g.GetFirstEndpoint(e)];
				int dw = ++degree[g.GetSecondEndpoint(e)];
				maxdegree = max(maxdegree,max(dv,dw));
			}
		}

		const bool REPORT_DEGREES = false;
		if (REPORT_DEGREES) {
			vector<int> degcount(maxdegree+1,0);
			for (int v=1; v<=n; v++) {
				degcount[degree[v]] ++;
			}
			fprintf (stderr, "Degrees:\n");
			for (int d=0; d<=maxdegree; d++) {
				if (degcount[d]>0) fprintf (stderr, "%3d: %5d\n", d, degcount[d]);
			}
			fprintf (stderr, "done reporting degrees\n");
			fflush(stderr);
		}

		
		
		
		//for (int e=1; e<=m; e++) keep[e] = true;
		//OneTwoPrep(g,keep);
		//PrepDegreeOne(g,keep,degree);
		
		//exit(-1);

		//PrepCorners(g,keep,degree);
		
		/*
		vector<int> x;
		int *y = &x[0];*/

		//bool *bkeep = &(keep[0]);
		bool OUTPUT_FILES = false;

		if (OUTPUT_FILES) {

			vector<bool> vkeep(n+1, true);
			for (int v=1; v<=n; v++) {
				if (degree[v] == 0) vkeep[v] = false;
			}

			fprintf (stderr, "About to draw...\n"); fflush(stderr);
			GraphDrawer::DrawSubgraph("test.gv", g, keep, vkeep);
			fprintf (stderr, "Drawn.\n"); fflush(stderr);
			OutputSTP("test.stp", g, keep);
		}

		Graph newg;
		CommitSubgraph (newg, g, keep);
		g = newg;



		/*
		Boruvka(solution, vordata, uf, null);
            Console.Error.WriteLine(" Boruvka found solution worth (value = {0})", solution.GetCost());

            STEdgeLinear dyntree = new STEdgeLinear(n);
            STEdgeLinear keytree = new STEdgeLinear(n);

            int added = 0;
            foreach (int e in solution.ElementEnumerator()) {
                int bv = vordata.GetBase(g.GetFirstEndpoint(e));
                int bw = vordata.GetBase(g.GetSecondEndpoint(e));
                if (bv != bw) {
                    added++;
                    //Console.Error.WriteLine("{0}/{1} ", bv, bw);
                    ArcCost cost = vordata.GetDistance(g.GetFirstEndpoint(e)) + vordata.GetDistance(g.GetSecondEndpoint(e)) + g.GetCost(e);
                    dyntree.Evert(bw);
                    dyntree.Link(bw, bv, -cost);
                }
            }

            Console.Error.WriteLine("Done building tree ({0} edges inserted, {1} in solution, {2} added).", dyntree.GetEdgeCount(), baselist.Count(), added);

            int eliminated = 0;
            bool[] keep = new bool[m + 1];
            int f;
            for (f = 1; f <= m; f++) {
                keep[f] = true;
                if (solution.Contains(f)) continue;
                int fv = g.GetFirstEndpoint(f);
                int fw = g.GetSecondEndpoint(f);
                ArcCost cost = g.GetCost(f);

                int bv = vordata.GetBase(fv);
                int bw = vordata.GetBase(fw);

                //Console.Error.WriteLine("Checking ({0},{1}) with bases ({2},{3}) {4} {5}.", fv, fw, bv, bw, g.IsTerminal(bv), g.IsTerminal(bw));

                //distance from v to closest terminal
                ArcCost tv = vordata.GetDistance(fv);
                if (tv > cost) continue;
                if (tv == cost && vordata.GetParentArc(fv) == f) continue;

                //distance from w to closest terminal
                ArcCost tw = vordata.GetDistance(fw);
                if (tw > cost) continue;
                if (tw == cost && vordata.GetParentArc(fw) == f) continue;

                //bottleneck distance inside the tree
                ArcCost tt = 0;
                if (bv != bw) {
                    dyntree.Evert(bw);
                    int x = dyntree.GetMinCost(bv);
                    //Console.Error.WriteLine("Querying {0}.", x);
                    tt = -dyntree.GetCost(x);
                    if (tt > cost) continue; //DO I NEED EQUALITY HERE?
                }

                eliminated ++;
                degree[fv]--;
                degree[fw]--;
                keep[f] = false;
                //Console.Error.WriteLine("Should be eliminating {0} ({1} {2} {3})", cost, tv, tw, tt);
            }

            PrepDegreeOne(keep, degree);

            int[] degcount = new int[n];
            int i;
            for (i=0; i<n; i++) {
                degcount[i] = 0;
            }

            for (v = 1; v <= n; v++) {
                if (!g.IsTerminal(v))
                degcount[degree[v]]++;
            }

            for (i=0; i<n; i++) {
                if (degcount[i] > 0) {
                    Console.Error.WriteLine("{0}: {1}", i, degcount[i]);
                }
            }

            Console.Error.WriteLine(" done preprocessing (value = {0}), {1}/{2} eliminated, {3} remaining.", 
                solution.GetCost(), eliminated, m, m - eliminated);

            bool simplify = true;
            if (simplify) {
                //remaining vertices
                int [] vlist = new int [n+1];
                int [] elist = new int [m+1];
                for (v=1; v<=n; v++) {
                    if (g.IsTerminal(v)) vlist[v] = 1;
                    else if (degree[v]>0) {vlist[v] = 2;} 
                    else {vlist[v] = 0;}
                }
                for (int e=1; e<=m; e++) {
                    if (keep[e]) {elist[e] = 2;} 
                    else {elist[e] = 0;}
                }
                WeightedGraph subgraph = new WeightedGraph(g,vlist,elist);
                Console.Error.WriteLine("Created subgraph, but did nothing with it!");
            }

            OutputSubgraph(keep);*/
	}

	// WARNING: IS THIS WORKING FOR A GRAPH WITH A SINGLE (TERMINAL) NODE?

	template <class GRAPH> static void OutputCompactSTP (const string &filename, GRAPH &g, vector<bool> &keep) {
		FILE *file = fopen (filename.c_str(), "w");
		int oldm = g.EdgeCount();
		int oldn = g.VertexCount();
		int newm = 0;

		int n = 0;
		vector<int> old2new (oldn+1,-1);
		for (int v=1; v<=oldn; v++) {
			if (g.GetDegree(v) > 0) {
				n++;
				old2new[v] = n;
			}
		}

		for (int e=1; e<=oldm; e++) {
			if (keep[e]) newm++;
		}

		fprintf (stderr, "Compact graph should have %d vertices and %d edges.\n", n, newm);
		fflush(stderr);

		
		fprintf (file, "33d32945 STP File, STP Format Version 1.00\n");
		fprintf (file, "Section Comment\n");
		fprintf (file, "End\n\n");

		fprintf (file, "Section Graph\n");
		fprintf (file, "Nodes %d\n", n);
		fprintf (file, "Edges %d\n", newm);
		for (int e=1; e<=oldm; e++) {
			if (!keep[e]) continue;
			int v,w;
			g.GetEndpoints(e,v,w);
			EdgeCost cost = g.GetCost(e);
			fprintf (file, "E %d %d %.0f\n", old2new[v], old2new[w], (double)cost);
		}
		fprintf (file, "End\n\n");

		fprintf (file, "Section Terminals\n");
		fprintf (file, "Terminals %d\n", g.TerminalCount());
		int counted = 0;
		for (int v=1; v<=oldn; v++) {
			if (g.IsTerminal(v)) {
				if (g.GetDegree(v) == 0) {
					fprintf (stderr, "ERROR! Terminal %d has degree 0.\n", v);
					fatal ("that's bad");
				}
				fprintf (file, "T %d\n", old2new[v]);
				counted ++;
			}
		}
		if (counted != g.TerminalCount()) {
			fprintf (stderr, "Found %d terminals, counted %d.\n", g.TerminalCount(), counted);
			fatal ("invalid number of terminals");
		}

		fprintf (file, "End\n\n");

		if (g.coord.CoordinatesAvailable()) {
			fprintf (file, "Section Coordinates\n");
			for (int v=1; v<=oldn; v++) {
				if (g.GetDegree(v) == 0) continue;
				double x, y;
				g.coord.GetCoordinates(v,x,y);
				fprintf (file, "DD %d %.0f %.0f\n", old2new[v], x, y);
			}
			fprintf (file, "End\n\n");
		}
		fclose(file);
	}


	template <class GRAPH> static void OutputSTP (const string &filename, GRAPH &g, vector<bool> &keep) {
		FILE *file = fopen (filename.c_str(), "w");
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int newm = 0;

		for (int e=1; e<=m; e++) {
			if (keep[e]) newm++;
		}
		
		fprintf (file, "33d32945 STP File, STP Format Version 1.00\n");
		fprintf (file, "Section Comment\n");
		fprintf (file, "End\n\n");

		fprintf (file, "Section Graph\n");
		fprintf (file, "Nodes %d\n", n);
		fprintf (file, "Edges %d\n", newm);
		for (int e=1; e<=m; e++) {
			if (!keep[e]) continue;
			int v,w;
			g.GetEndpoints(e,v,w);
			EdgeCost cost = g.GetCost(e);
			fprintf (file, "E %d %d %.0f\n", v, w, (double)cost);
		}
		fprintf (file, "End\n\n");

		fprintf (file, "Section Terminals\n");
		fprintf (file, "Terminals %d\n", g.TerminalCount());
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) fprintf (file, "T %d\n", v);
		}
		fprintf (file, "End\n\n");

		if (g.coord.CoordinatesAvailable()) {
			fprintf (file, "Section Coordinates\n");
			for (int v=1; v<=n; v++) {
				double x, y;
				g.coord.GetCoordinates(v,x,y);
				fprintf (file, "DD %d %.0f %.0f\n", v, x, y);
			}
			fprintf (file, "End\n\n");
		}
		fclose(file);
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
			root = PickRandomTerminal(*g,random);
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
			root = PickRandomTerminal(*g, random);
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



	static void MarkSolutionVertices (SteinerSolution &solution, UniverseSet &svertices) {
		int n = solution.g->VertexCount();
		for (int v=1; v<=n; v++) {
			if (solution.GetDegree(v)>0) svertices.Insert(v);
		}
	}

        /// <summary>
        /// Steiner-vertex insertion local search.
        /// </summary>
        /// <param name="solution">Initial solution (may change).</param>
        /// <param name="maxv">Maximum number of vertices checked.</param>
        /// <returns></returns>
        static int VertexInsertion(Graph &g, SteinerSolution &solution, int maxv, RFWLocalRandom &random) {
            int n = g.VertexCount();
            int v = 0;

			bool changed = false;
            bool CHANGE_ON_TIE = false;
            const bool debug_mode = false;
            const bool randomize = true; //randomize order in which vertices are tested
            const bool high_verbose = false;
            const bool verbose = false;
            if (verbose) {
                fprintf (stderr, "Running vertex-insertion local search.");
                fprintf (stderr, "Solution costs %.0f.", solution.GetCost());
            }

			RFWStack<LocalArc> removed (n); //edges temporarily removed from the tree
			RFWStack<int> inserted (n); //edges temporarily inserted into the tree

            //list of all vertices in the solution
            UniverseSet svertices(n); 
            MarkSolutionVertices(solution, svertices);
			if (verbose) fprintf (stderr, "Solution has %d vertices.\n", svertices.Count());

            //compute number of solution neighbors for each nonsolution vertex
            int *solneighbors = new int[n + 1]; //this is zero for terminals!
            for (v=1; v<=n; v++) {solneighbors[v] = 0;}
			int p, pend;
			for (svertices.GetBounds(p,pend); p<pend; p++) {
				int v = svertices.PickPos(p);
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
					if (!svertices.Contains(w)) solneighbors[w]++;
				}
			}


			if (verbose) fprintf (stderr, "Inserting initial edges to the tree.\n");
            EdgeCost cost0 = solution.GetCost();
            STEdgeLinear<EdgeCost> dyntree (n);
            //foreach (int alabel in solution.ElementEnumerator()) {
			int m = g.EdgeCount();
			for (int e=1; e<=m; e++) { // WARNING: LOOPING OVER ALL EDGES
				if (!solution.Contains(e)) continue;
                int a, b;
                g.GetEndpoints(e, a, b);
                EdgeCost cost = g.GetCost(e);
                if (debug_mode && (!svertices.Contains(a) || !svertices.Contains(b))) fatal ("Fishy!");
                dyntree.Evert(a);
                dyntree.Link(a, b, -cost);
            }



            int impmoves = 0; //number of improving moves
            EdgeCost curcost = cost0; //current solution value
            int failures = 0; //number of consecutive failures to improve

			// create a random permutation of all vertices
			// (there should be a separate function for this)
            int *perm = new int[n + 1];
            for (v = 1; v <= n; v++) {perm[v] = v;}
            if (randomize) {
                for (v = 1; v < n; v++) {
                    int w = random.GetInteger(v, n);
                    int t = perm[v];
                    perm[v] = perm[w];
                    perm[w] = t;
                }
            }

			// prepare for the main loop
            int pos = 0;
            v = 1;
            int count = 0;

			const bool TIME_EVERTS = false;
			RFWTimer evertime;
			if (TIME_EVERTS) {
				evertime.start();
				evertime.pause();
			}


			while (failures < n) {
                if (count >= maxv) break; //enough tries---just stop
                count ++;
                if (++pos > n) pos = 1;
                v = perm[pos];                

                //consider only nonterminals with at least two neighbors in the solution
                if (solneighbors[v] <= 1) {failures++; continue;}

                if (debug_mode && svertices.Contains(v)) fprintf (stderr, "Should not have been here (%d, solneighbors %d).", v, solneighbors[v]);
                if (verbose) fprintf (stderr, "Processing %d: ", v);

                EdgeCost increase = 0; //increase in cost
                int inscount = 0; //number of insertions
                int tries = 0;
                bool first = true;

                //invariant: the new vertex (v) will be the root during these operations
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
                    int w = a->head; 
                    if (!svertices.Contains(w)) continue; //only care about neighbors in the solution
                    tries++;
                    EdgeCost cost = a->cost; 

					if (first) { //first edge: just insert it
						if (TIME_EVERTS) evertime.resume();
						//fprintf (stderr, ">");
                        dyntree.Evert(w); //make w the root of its tree
                        if (TIME_EVERTS) evertime.pause();
						if (debug_mode && (dyntree.GetRoot(v) == dyntree.GetRoot(w))) fatal ("v cannot be on the same tree as w");
                        dyntree.Link(w, v, -cost); //v is the new root
                        if (verbose) fprintf (stderr, "+%.2f I(%d,%d) (r:%d)", (double)cost, v, w, dyntree.GetRoot(v));
                        increase += cost;
                        inscount++;
                        inserted.push(a->label); // remember edge (for later removal)
                        first = false;
                    } else {
                        //other edges: insert only if they improve something
                        if (debug_mode && (dyntree.GetRoot(w) != v)) {
                            fprintf (stderr, "In solution: %d.", svertices.Contains(w));
                            fprintf (stderr, "v:%d e:%d", svertices.Count(), -1); //solution.Count());
                            fprintf (stderr, " %d:%d:%d:%d", v, dyntree.GetRoot(v), w, dyntree.GetRoot(w));
                            fatal ("tree does not have the right root");
                        }
                        int z = dyntree.GetMinCost(w); //get mincost edge (z,pz) on the w-v path
                        if (debug_mode) {
                            if (z == 0) fatal ("path should not be empty");
                            if (dyntree.GetRoot(z) == z) fprintf (stderr, "%d is its own root, from v=%d and w=%d and root(w)=%d!", z, v, w, dyntree.GetRoot(w));
                        }
                        EdgeCost zcost = -dyntree.GetCost(z);

						//we might want to remove edge (z,pz) with cost zcost
                        if (zcost > cost) {
                            if (verbose) fprintf (stderr, " %.2f->%.2f", (double)zcost, (double)cost);
                            int pz = dyntree.GetParent(z);

                            //if we are removing an original edge, remember it
                            //(it may need to be reinserted later)
                            if (pz != v) {removed.push(LocalArc(z,pz,zcost));}
                            if (verbose) fprintf(stderr, " R(%d,%d)", z, pz);

                            //remove arc (z,pz)
                            EdgeCost t = -dyntree.cut(z); //remove the edge
                            if (debug_mode && (t != zcost)) fatal ("Internal costs are not consistent");
							
							if (TIME_EVERTS) evertime.resume();
                            dyntree.Evert(w); //make w the root of its component
                            if (TIME_EVERTS) evertime.pause();
							dyntree.Link(w, v, -cost); //add arc from w to v; v remains root

                            if (debug_mode && (dyntree.GetRoot(v) != v)) fatal("tree does not have the right root");

                            increase += (cost - zcost);
                            inscount++;
                            inserted.push(a->label);  //remove the new edge
                            if (verbose) fprintf (stderr, " I(%d,%d):%d ", w, v, a->label); 
                        }
                    }
                }

				//fprintf (stderr, "%d ", tries);

                if (verbose) fprintf (stderr, " = %.2f after %d insertions (tries=%d degree=%d).", increase, inscount, tries, solneighbors[v]);

                // if the new solution is worse, we must restore the original one
                bool restore;
                if (increase <= -EDGE_COST_PRECISION) { //solution improved!
                    impmoves++;
                    failures = 1; //not zero because we don't need v itself
                    restore = false;
                } else { //same or worse
                    restore = (increase > EDGE_COST_PRECISION || !CHANGE_ON_TIE);
                    failures++;
                }

                //if the solution did not improve, restore original tree
                if (restore) {
                    //remove the newly inserted edges; they are all neighbors of v
                    //(note that some may have been removed already)
                    while (!inserted.isEmpty()) {
                        int alabel = inserted.pop();
                        int w = g.GetOther(alabel, v);

                        if (verbose) {
                            fprintf (stderr, "ToRem%d ", alabel);
                            int tx, ty;
                            g.GetEndpoints(alabel, tx, ty);
                            fprintf (stderr, " [%d,%d]", tx, ty);
                        }

                        dyntree.Evert(w); //WARNING: ISN'T IT ENOUGH TO KEEP V AS THE ROOT ALL THE WAY? IT SHOULD BE THE IMMEDIATE PARENT OF EVERYBODY.

                        //if the edge is still there, remove it
                        if (dyntree.GetParent(v) == w) {
                            dyntree.cut(v);
                            if (verbose) fprintf (stderr, "deleting (%d,%d) ", v, w);
                        } else {
                            if (verbose) fprintf(stderr, " @@@ (v:%d pv:%d w:%d) ", v, dyntree.GetParent(v), w);
                        }
                    }

                    //reinsert temporarily removed edges
                    while (!removed.isEmpty()) {
                        LocalArc larc = removed.pop(); //why can't we remember just integers?
                        if (TIME_EVERTS) evertime.resume();
						dyntree.Evert(larc.v); //make one of the endpoints the root
                        if (TIME_EVERTS) evertime.pause();
						dyntree.Link(larc.v, larc.w, -larc.cost); //link one to the other
                        if (verbose) fprintf (stderr, "R");
                    }

                    if (verbose) {
                        int dyncount = dyntree.GetEdgeCount();
                        int solcount = -1; //solution.Count();
                        fprintf (stderr, "Dynamic tree now has {0} edges. Solution has {1} edges.", dyntree.GetEdgeCount(), -1); //solution.Count());
                    }
                } else {
                    //we're going to the neighbor!
                    curcost += increase;
                    removed.reset(); //no need to remember just-removed edges
                    inserted.reset(); //no need to remember just-inserted edges
                    svertices.Insert(v); //solution has one new vertex
                    solneighbors[v] = 0; //v now belongs to the solution
					changed = true;

                    // we will NOT change the list of edges in the solution until we 
					// we reach a local optimal; we only update auxiliary data structures as we go
					SPGArc *a, *end;
					for (g.GetBounds(v,a,end); a<end; a++) {
						int w = a->head;
                        if (!svertices.Contains(w)) {solneighbors[w]++;} //fix degree
                    }
                }

                if (high_verbose && increase < -EDGE_COST_PRECISION) {
					fprintf (stderr, "Vertex %d improves the solution by %.2f.", v, (double)-increase);
                }
            }
            
            if (high_verbose) {
                fprintf (stderr, "End of local search; found %d improving moves.", impmoves);
                fprintf (stderr, "Solution went from %.1f to %.1f.\n", cost0, curcost);
            }

            EdgeCost improvement = cost0 - curcost;

			if (TIME_EVERTS) fprintf (stderr, "%.2f ", 1000*evertime.getTime());

            if (improvement > EDGE_COST_PRECISION) {
				if (!changed) fatal ("inconsistent assessment");
                RebuildSolution(dyntree, solution);
            } else {
				if (changed) fprintf (stderr, "Not using changed version.\n");
			}

            if (verbose) fprintf (stderr, "Found %.2f, mst is %.2f.", (double)curcost, (double)solution.GetCost()); 

			if (solneighbors) delete [] solneighbors;
			if (perm) delete [] perm;
            return impmoves;
        }

        /// <summary>
        /// Rebuild a solution based on the dynamic tree.
        /// </summary>
        static void RebuildSolution (STEdgeLinear<EdgeCost> &dyntree, SteinerSolution &solution) {
            bool verbose = false;
            if (verbose) fprintf (stderr, "Rebuilding solution from dynamic tree.");
			Graph *g = solution.g;
            int n = g->VertexCount();
            solution.Reset();

            for (int v=1; v<=n; v++) {
                int w = dyntree.GetParent(v);
                if (w == 0) continue;
                EdgeCost cost = -dyntree.cut(v);

				// find a matching edge
				SPGArc *a, *end;
				for (g->GetBounds(v,a,end); a<end; a++) {
					if (a->head == w) {
						int e = a->label;
						if (a->cost == cost) { //just in case there are parallel edges
                            solution.Insert(e);
                            break;
                        }
                    }
                }
            }

            if (verbose) fprintf (stderr, "Solution costs %d.", solution.GetCost());
        }

		struct DFSData {
            int *id2dfs;  //maps vertex ids to post-order numbers
            int *dfs2id;  //maps post-order numbers to vertex ids
            int *id2parc; //maps vertex ids to parent arcs in the dfs tree

            DFSData(int n) {
                id2dfs  = new int [n+1];
                dfs2id  = new int [n+1];
                id2parc = new int [n+1];
            }

			~DFSData() {
				delete [] id2parc;
				delete [] dfs2id;
				delete [] id2dfs;
				//fprintf (stderr, "Deleted DFS data.\n");
				//fflush(stderr);
			}
        };

        /// <summary>
        /// Compute nearest common ancestors of all edges with both ends in the solution. 
        /// Puts edge (v,w) in bucket NCA(v,w).
        /// </summary>
        /// <param name="hzt"></param>
        /// <param name="svertices"></param>
        /// <param name="dfsdata"></param>
        /// <param name="root"></param>
        static void ComputeNCA(Graph &g, StaticBuckets &hzt, UniverseSet &svertices, DFSData &dfsdata, int root) {
            bool verbose = false;
            int *id2dfs = dfsdata.id2dfs;
            int *id2parc = dfsdata.id2parc;
            int *dfs2id = dfsdata.dfs2id;

            int v, n = g.VertexCount();
            int *ancestor = new int [n+1]; //auxiliary variable used for nca computations...
            UnionFind ncauf (n); //auxilary data structure for nca computation

            for (v=1; v<=n; v++) ancestor[v] = v;    //and everybody is its own ancestor
            hzt.Reset(); //reset list of visited nodes

            int nhzt = 0; //number of horizontal edges

            //traverse nodes in bottom-up order, ending up in the root
            int count = 0;
            for (int vdfs=1; vdfs<=n; vdfs++) {
                v = dfs2id[vdfs];
				if (verbose) fprintf (stderr, "%d:%d\n", vdfs, v);
                if (v<=0) break;

                //for (int pa = g.GetStart(v); pa < g.GetEnd(v); pa++) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
                    int w = a->head; //g.GetArcHead(pa);

                    //look at (v,w) if w in solution and already processed
                    if (svertices.Contains(w) && (id2dfs[w] < vdfs)) {
                        int nca = ancestor[ncauf.Find(w)]; 
						if (nca<1 || nca>n) fatal ("nca out of range");
                        if (v==nca || w==nca) continue;
                        nhzt++;
                        hzt.Insert(a->label, nca); //g.GetArcLabel(pa), nca);
						if (verbose) fprintf (stderr, "Inserting (%d,%d) into bucket %d.\n", v,w, nca);
                    }
                }

                //unite v with its parent
                if (v!=root) {
                    int parent = g.GetOther(id2parc[v],v);
                    ncauf.Union(parent,v);
                    ancestor[ncauf.Find(v)] = parent; //parent is the ancestor
                }
            }

			delete [] ancestor;
            //edges per horizontal edge
            if (verbose) fprintf (stderr, "%d edges are horizontal, total count is %d.", nhzt, count);
        }




        /// <summary>
        /// Perform DFS on the solution, numbering vertices in reverse post-order.
        /// Returns the number of vertices visited.
        /// </summary>
        /// <param name="r">Root of DFS.</param>
        /// <param name="solution">Current solution.</param>
        /// <param name="dfs2id">Output: map from dfs number to id (-1 if not visited)</param>
        /// <param name="id2dfs">Output: map from id to dfs number (-1 if not visited)</param>
        /// <param name="id2parc">Output: map from it to parent arc (0 if not visited)</param>
        static int DFS (Graph &g, int r, SteinerSolution &solution, DFSData &dfsdata, RFWStack<int> &stack) {
            // this is a funny implementation of dfs: when we first scan a vertex, we simply add to the stack
            // every nonscanned neighbor---even those that are already in the stack. This requires a stack of size m.       
			// WARNING! IF WE ARE ONLY SCANNING EDGES OF THE SOLUTION, THE SIZE IS N
            int *id2dfs = dfsdata.id2dfs;
            int *dfs2id = dfsdata.dfs2id;
            int *id2parc = dfsdata.id2parc;

            int n = g.VertexCount();
            int m = g.EdgeCount();

            stack.reset();

            //id2dfs: -1:unreached  0:scanned  >0:processed
            for (int v=0; v<=n; v++) {
                id2dfs[v] = -1; //everybody unreached, initially
                dfs2id[v] = -1;
                id2parc[v] = 0;
            }

            stack.push(r);
            int nextdfs = 1;

            while (!stack.isEmpty()) {
                int v = stack.pop();
                int vdfs = id2dfs[v];
                if (vdfs > 0) {continue;} //vertex already processed: nothing else to do

                //vertex already scanned, but with no label; we assign it a label now
                if (vdfs == 0) {
                    id2dfs[v] = nextdfs;
                    dfs2id[nextdfs] = v;
                    nextdfs++;
                    continue;    
                }

                //vertex not yet scanned: scan it, put it back on the stack (a label will be assigned later)
                stack.push(v);
                id2dfs[v] = 0;

                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
				SPGArc *a, *end;
                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
				for (g.GetBounds(v,a,end); a<end; a++) {
                    int alabel = a->label; //g.GetArcLabel(pa);
                    if (!solution.Contains(alabel)) continue;
                    int w = a->head; //g.GetArcHead(pa);//arc.head;
                    if (id2dfs[w] >= 0) continue; //w already scanned: no need to go there again
                    id2parc[w] = alabel;
                    stack.push(w);
                }
            }
            return nextdfs - 1;
        }


        //mark the components containing the elements of the stack
        static void MarkComponent(Graph &g, RFWStack<int> &stack, int *id2parc, UniverseSet &marked) {
            //Console.Error.Write("+");
            //Invariant: a vertex becomes marked when it is inserted into the stack
            //A marked vertex is or was on the stack.

            bool verbose = false;
            if (verbose) fprintf (stderr, "Marking component from %d vertices; marked has %d.", stack.getNElements(), marked.Count());

            //make invariants true for original vertices
            for (int i = stack.getNElements(); i >= 1; i--) {
                int v = stack.peek(i);
                if (marked.Contains(v)) fprintf (stderr, "BAD");
                marked.Insert(v);
            }

            int mcount = marked.Count();

            //add all relevant tree children to the stack
            while (!stack.isEmpty()) {
                int v = stack.pop();
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
                    int w = a->head; //g.GetArcHead(pa);
                    if (id2parc[w]==a->label) { //g.GetArcLabel(pa)) {
                        if (marked.Insert(w)) {
                            //mcount ++;
                            stack.push(w);
                        }
                    }
                }
                /*
                foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                    int w = arc.head;
                    if (id2parc[w]==arc.label) {
                        if (marked.Insert(w)) stack.Push(w);
                    }
                }*/
            }

            mcount = marked.Count() - mcount;
            if (verbose) fprintf(stderr, "%d elements marked", mcount);

            /* if (verbose) {Console.Error.WriteLine(" {0} elements marked.", marked.Count());
                foreach (int e in marked.ElementEnumerator()) {Console.Error.Write("+{0} ", e);}}*/
        }



        /*------------------------------------------
         * Local search based on vertex elimination 
         *-----------------------------------------*/
        static int VertexElimination(Graph &g, SteinerSolution &solution, RFWLocalRandom &random) {
            bool verbose = false;
            if (verbose) {
				fprintf (stderr, "\n\nRunning local search based on vertex eliminations.\n");
				fflush(stderr);
			}
            //bool JUST_COUNT = false;
            int improvements = 0;
            int n = g.VertexCount();
            int m = g.EdgeCount();
            int v;

            DFSData dfsdata(n);
            int *id2dfs = dfsdata.id2dfs;
            int *dfs2id = dfsdata.dfs2id;
            int *id2parc = dfsdata.id2parc;

            //bool [] forbidden = new bool [n+1];
            UniverseSet forbidden(n);// = new UniverseSet(n);
            bool *pinned = new bool [n+1];
            for (v=1; v<=n; v++) {pinned[v] = false;}


            //StaticStack<int> stack = new StaticStack<int>(m);
			RFWStack<int> stack(m); 

            //UniverseSet svertices = new UniverseSet(n); //vertices in the original solution
			UniverseSet svertices(n); //vertices in the original solution

            StaticBuckets hzt (m, n); //horizontal edges
            
			UnionFind uf (n); //union find data structure
			UnionFind kuf (n); // union find for Kruskal computation
			BinaryHeap<EdgeCost> kheap(m); //binary heap for Kruskal computation
			
			//UnionFind uf = new UnionFind(n); //union find data structure
            //UnionFind kuf = new UnionFind(n); //union find for kruskal computation
            //BinaryHeap<ArcCost> kheap = new BinaryHeap<ArcCost>(m); //binary heap for kruskal computation

            UniverseSet cedges (m);

            int root = PickRandomTerminal(g,random);
			//fprintf (stderr, "r%d ", root);

			//root = 5;
			if (verbose) fprintf (stderr, "Root is %d.\n", root);

            //profit from each possible move
            EdgeCost *profit = new EdgeCost[n + 1];
            for (v=1; v<=n; v++) profit[v] = 0;

            //CheckSolution(solution);
            int visited = DFS(g, root, solution, dfsdata, stack);
            if (verbose) fprintf (stderr, "DFS visited %d vertices from root %d.", visited, root);

            if (visited < g.TerminalCount()) {
                fprintf (stderr, "Visited %d vertices, but there are %d terminals.", visited, g.TerminalCount());
                fatal ("Something wrong.");
            }

            //find list of all original solution vertices
            MarkSolutionNodes(g, solution, svertices);


            //find set of horizontal edges
            if (verbose) {
				fprintf (stderr, "Solution has %d vertices.\n", svertices.Count());
				fflush(stderr);
			}

            ComputeNCA(g, hzt, svertices, dfsdata, root); //id2dfs, id2parc);


			if (verbose) {
				fprintf (stderr, "Done computing NCAs.\n");
				fflush(stderr);
			}

            //data about subtrees
            PairingHeap<EdgeCost> heap(2 * m, n);
            UniverseSet subtrees (n);

            bool high_verbose = false; //verbose;

            int solsize = svertices.Count();
			if (verbose) {
				fprintf (stderr, "Initiating main loop.\n");
				fflush(stderr);
			}

            for (int i = 1; i < solsize; i++) {
                v = dfs2id[i];

                if (verbose) {
					fprintf (stderr, "\n\nProcessing vertex %d (dfs=%d, pinned=%d): ", v, id2dfs[v], pinned[v]);
					fflush(stderr);
				}

                //------------------------------------
                // find original subtrees of vertex v
                //------------------------------------
                EdgeCost remcost = 0; //total cost of edges removed
                remcost += g.GetCost(id2parc[v]); //arc to the parent
                subtrees.Reset();

				if (high_verbose) {
					fprintf (stderr, "Scanning."); fflush(stderr);
				}

                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
                    int alabel = a->label; //g.GetArcLabel(pa);
                    int w = a->head; //g.GetArcHead(pa);

					if (high_verbose) {
						fprintf (stderr, "Looking at (%d,%d)... ", v, w);
						fflush(stderr);
					}

                    //look at all *original* children (some children may have been removed)
                    if (id2parc[w] == alabel) {
                        remcost += a->cost; //g.GetArcCost(pa); //arc.cost;
                        subtrees.Insert(uf.Find(w));
                        if (verbose) {
                            int f = uf.Find(w);
							if (verbose) {fprintf (stderr, "Vertex %d corresponds to set %d with size %d ", w, f, uf.Size(f)); fflush(stderr);}
                        }
                    }
					if (high_verbose) {
						fprintf (stderr, "End of %d.\n", w);
						fflush (stderr);
					}
                }

                if (high_verbose) {
					fprintf (stderr, "%d children, cost cost is %.2f, terminal %d.\n", subtrees.Count(), remcost, g.IsTerminal(v));
					fflush(stderr);
				}
                    
                //if (!g.IsTerminal(v)) {
                //    if (pinned[v]) Console.Error.Write("PINNED({0})", v);
                //    else Console.Error.Write(".");
                //}

                if (g.IsTerminal(v) || pinned[v]) {
                    //if (pinned[v]) Console.Error.Write(".");
                    if (verbose) {
                        if (pinned[v]) fprintf (stderr, "PINNED(%d)", v);
                        fprintf (stderr, "T! ");
						fflush (stderr);
                    }
                } else {
                    //we may want to remove this vertex (not pinned or a terminal)
                    if (verbose) {
						fprintf (stderr, "R! ");
						fflush(stderr);
					}

                    //penalize edges incident to forbidden nodes?
                    //bool USE_FORBIDDEN_PENALTIES = false;
                    //double badpenalty = 0;

                    //list of edges for Kruskal's algorithm
                    kheap.Reset();

                    //process each subtree: clean up heap to find cheapest edge that is actually
                    //vertical, then insert it into kheap (the list of candidate edges).
                    //foreach (int s in subtrees.ElementEnumerator()) {
					int p, pend;


					for (subtrees.GetBounds(p,pend); p<pend; p++) {
						int s = subtrees.PickPos(p);

                        if (verbose) {
							fprintf (stderr, "Clean up heap %d... ", s);
							fflush(stderr);
						}
                        int x=0, y=0;
                        EdgeCost value = 0;
                        int a = 0;

                        //look for an arc whose head does not belong to any of the subtrees---or to v
                        while (!heap.IsEmpty(s)) {
                            a = heap.FindMin(s, value);
                            g.GetDirectedEndpoints(a, x, y);
                            if (verbose) {
								fprintf (stderr, "a=%d x=%d y=%d m=%d\n", a, x, y, m); fflush(stderr);
							}
                            if (uf.Find(x) != uf.Find(s)) {
								fprintf (stderr, "x:%d y:%d s:%d\n", uf.Find(x), uf.Find(y), uf.Find(s));
								fatal ("Arc's tail should belong to subtree.");
							}

                            //an edge to v or to a subtree (the same or another) is no longer vertical
                            if (y == v || subtrees.Contains(uf.Find(y))) {
                                if (verbose) {
									fprintf (stderr, ".");
									fflush(stderr);
								}
                                heap.DeleteMin(s, value); //get rid of the arc
                                a = 0;
                                continue;
                            } else break; //done cleaning
                        }

						if (verbose) {fprintf (stderr, "Done with heap.\n"); fflush(stderr);}



                        if (a == 0) {
                            if (verbose) {fprintf (stderr, "Empty."); fflush(stderr);}
                        } else {
                            /*
                            if (USE_FORBIDDEN_PENALTIES) {
                                if (forbidden.Contains(x) || forbidden.Contains(y)) badpenalty = 1/(double)(n+1);
                            }*/

							if (verbose) {
								int av, aw;
								g.GetDirectedEndpoints(a, av, aw);
								fprintf (stderr, "Inserting arc %d=(%d,%d) with cost %d into kheap.", a, av, aw, value); 
								fflush(stderr);
							}
                            if (a > m) a -= m;
                            //kheap.Insert(a, value);
                            //kheap.PushBack(a,value + badpenalty);
                            kheap.PushBack(a,value);
							if (verbose) {
								fprintf (stderr, "Done pushing back %d.", a); 
								fflush(stderr);
							}
                        }
						if (verbose) {fprintf (stderr, "%d-done ", s); fflush(stderr);} 
                    }

                    //insert horizontal edges into heap
                    //foreach (int ha in hzt.ElementEnumerator(v)) {

					// WARNING: CHECK THAT THIS IS DOING THE RIGHT THING
					for (int ha=hzt.GetFirst(v); ha>=0; ha=hzt.GetNext(ha)) {

                        /*if (USE_FORBIDDEN_PENALTIES) {int x, y; g.GetDirectedEndpoints(ha, out x, out y);
                            if (forbidden.Contains(x) || forbidden.Contains(y)) badpenalty = 1/(double)(n+1);
                        }*/
                        //kheap.PushBack(ha, g.GetCost(ha) + badpenalty);
						if (verbose) {
							int x,y;
							g.GetEndpoints(ha, x, y);
							fprintf (stderr, "Adding edge %d=(%d,%d) to kheap.\n", ha, x, y);
						}
                        kheap.PushBack(ha, g.GetCost(ha));
                    }
                    kheap.Heapify();

                    if (verbose) {
						fprintf (stderr, "There are %d edges in total in kheap.", kheap.GetSize());
						fflush(stderr);
					}



                    //run kruskal on modified graph
                    EdgeCost altcost = 0;
                    int picked = 0;
                    cedges.Reset();
                    int goal = subtrees.Count();
                    while (!kheap.IsEmpty()) {
                        EdgeCost cost;
                        unsigned int a; //WARNING: MADE THIS UNSIGNED
                        kheap.RemoveFirst(a, cost);
                        int x, y;
                        g.GetEndpoints(a, x, y);

						if (verbose) fprintf (stderr, "h(%d,%d) ", x, y);

                        if (forbidden.Contains(x) || forbidden.Contains(y)) {
                            continue;
                            //break;
                            //continue; //WARNING: BY DOING THIS, WE ARE NO LONGER GUARANTEED TO FIND AN MST
                        }

                        int gx = uf.Find(x);
                        if (!subtrees.Contains(gx)) gx = v;
                        int gy = uf.Find(y);
                        if (!subtrees.Contains(gy)) gy = v;

                        //x and y are now the representatives of their original components
                        if (kuf.Find(gx) != kuf.Find(gy)) {
                            cedges.Insert(a);
                            kuf.Union(gx, gy);
                            altcost += cost;
                            picked++;
                            if (picked == goal) break; //done building the tree
                        }
                    }
                    if (high_verbose) {
						fprintf (stderr, "Alternative has %d edges and cost %.2f (subtrees:%d).\n", picked, altcost, subtrees.Count());
						fflush(stderr);
					}
                    
                    if (verbose && picked!=goal) {fprintf(stderr, "VEDISCONNECTED");}
                    
                    if (picked == goal) {
                        //Console.Error.WriteLine("GOT HERE! NOW TRYING!");
                        if (high_verbose) fprintf (stderr, "Cost: %.2f -> %.2f.", (double)remcost, (double)altcost);
						double impcost = remcost - altcost; //improvement
                        if (impcost > EDGE_COST_PRECISION) { 
                            //found an improving move! change the solution
                            improvements++;

                            bool local_verbose = verbose;
                            if (local_verbose) fprintf(stderr, "I%.2f->%.2f %.2f->", (double)remcost, (double)altcost, (double)solution.GetCost());

                            //remove all neighboring edges
							SPGArc *a, *end;

							if (local_verbose) solution.Output(stderr);

							int removals = 0;
                            //for (int pa = g.GetStart(v); pa < g.GetEnd(v); pa++) {
							for (g.GetBounds(v,a,end); a<end; a++) {
								if (solution.Remove(a->label)) { //will do nothing if the solution is not there
									removals ++;
									if (local_verbose) {
										fprintf (stderr, "Removing %d=(%d,%d) with cost %d.\n", a->label, v, a->head, a->cost);
										if (solution.Contains(a->label)) fatal ("Something bad!\n");
									}
								}
                            }
							
							if (local_verbose) fprintf (stderr, "Removed vertex %d with %d incident edges (remcost=%d altcost=%d).\n", v, removals, remcost, altcost);

                            //add new edges
                            //foreach (int newedge in cedges.ElementEnumerator()) {
							int insertions = 0;
							int p, pend;
							for (cedges.GetBounds(p,pend); p<pend; p++) {
								int newedge = cedges.PickPos(p);
                                int x, y;
                                g.GetEndpoints(newedge, x, y);
                                pinned[x] = pinned[y] = true;
                                solution.Insert(newedge);
								insertions ++;
								if (local_verbose) fprintf (stderr, "Inserting %d=(%d,%d).\n", newedge, x, y);
                            }
							if (removals != insertions+1) fatal ("inconsistent operation");
                            stack.reset();
                            stack.push(v);
							if (local_verbose) {fprintf (stderr, "Marking component.\n"); fflush(stderr);}
                            MarkComponent(g, stack, id2parc, forbidden);

							if (local_verbose) {fprintf (stderr, "%.2f  ", (double)solution.GetCost()); fflush(stderr);} 

							if (verbose) {
								fprintf (stderr, "Checking after move...\n");
								CheckSolution(g, solution);
							}
                            //CheckSolution(solution);
                            profit[v] = remcost - altcost;

                                //Console.Error.WriteLine ("Profit: {0}", profit[v]);
                            //}
                        } 
                    }
                    else {
                        if (high_verbose) {
							fprintf (stderr, "TREE IS NOT EVEN CONNECTED.");
							fflush(stderr);
						}
                    }

					if (verbose) {fprintf (stderr, "Resetting union-find.\n"); fflush(stderr);}

                    //reset the union-find data structure used for this; will explicitly change only the affected elements
                    kuf.Reset(v); 
                    //foreach (int ks in subtrees.ElementEnumerator()) kuf.Reset(ks);
					{
						int p, pend;
						for (subtrees.GetBounds(p,pend); p<pend; p++) {
							int ks = subtrees.PickPos(p);
							kuf.Reset(ks);
						}
					}
                }


                //---------------------------
                // PREPARE FOR NEXT ITERATION
                //---------------------------

                /*---
                 | create heap containing all vertical edges incident into v
                 *--*/
                int insertions = 0;
                //foreach (WeightedGraph.Arc varc in g.ArcEnumerator(v)) {
                //for (int pa = g.GetStart(v); pa < g.GetEnd(v); pa++) {
				{
					SPGArc *a, *end;
					for (g.GetBounds(v,a,end); a<end; a++) {

						int w = a->head; //g.GetArcHead(pa); //varc.head;
						if (!svertices.Contains(w)) continue; //edge must link to other original solution vertex
						if (subtrees.Contains(uf.Find(w))) continue; //edge must not go to subtree
						int alabel = a->label; //g.GetArcLabel(pa); //varc.label;
						if (id2parc[v] == alabel) continue; //edge must not go to original parent

						//get unique label for the directed version of the arc
						int dlabel = (v < w) ? alabel : alabel + m;
						//if (verbose) {fprintf (stderr, "i(%d,%d) with label %d (alabel:%d).", v, w, dlabel, alabel); fflush(stderr);} 
						heap.Insert(v, dlabel, a->cost); //g.GetArcCost(pa)); //varc.cost);
						if (verbose) fprintf (stderr, "Inserted (%d,%d) into heap %d.\n", v, w, v);
						insertions++;
					}
				}

                if (verbose) fprintf(stderr, "Made %d insertions to heap %d.\n", insertions, v);

                //unite heap for v and its subtrees
                int rep = v; //new representative of the heap
                //foreach (int subrep in subtrees.ElementEnumerator()) {
				int p, pend;
				for (subtrees.GetBounds(p,pend); p<pend; p++) {
					int subrep = subtrees.PickPos(p);
                    uf.Union(rep, subrep);
					int t = uf.Find(rep);
					if (verbose) fprintf (stderr, "United %d and %d into %d (size %d).\n", rep, subrep, t, uf.Size(t));
                    if (t == rep) {
                        heap.Merge(rep, subrep);
                    } else {
                        heap.Merge(subrep, rep);
                        rep = t;
                    }
                }

				if (verbose) {fprintf (stderr, "heap(%d)=%d ", v, rep); fflush(stderr);}
            }

            if (verbose) {
                fprintf (stderr, "LS found %d improvements.", improvements);
				fflush(stderr);
                //fprintf (stderr, "Found %d possible candidates for removal.", improvements);
            }

            //Console.Error.WriteLine("FOUND {0} IMPROVEMENTS.", improvements);

            //throw new Exception ("dying");

			delete [] profit;
			delete [] pinned;
			//fprintf (stderr, "Done with local search with cost %d.\n", solution.GetCost());
			//fflush(stderr);
            return improvements;
        }


        static void CheckSolution(Graph &g, SteinerSolution &solution) {
            if (!InnerCheck(g, solution, false)) {
                fatal ("Invalid solution");
            }
        }

        static bool InnerCheck(Graph &g, SteinerSolution &solution, bool verbose) {
            int n = g.VertexCount();
            UniverseSet svertices(n);
            UnionFind uf(n);
            MarkSolutionNodes(g, solution, svertices);
            int ncomp = svertices.Count();

            UniverseSet terminals(n);
            for (int t=1; t<=g.VertexCount(); t++) {
                if (g.IsTerminal(t)) terminals.Insert(t);
            }

            //Console.Error.WriteLine("Term

            //foreach (int e in solution.ElementEnumerator())

			int m = g.EdgeCount();
			int ecount = 0;
			for (int e=1; e<=m; e++) {
				if (!solution.Contains(e)) continue;	
				ecount ++;
                int v, w;
                g.GetEndpoints(e, v, w);
                if (!svertices.Contains(v) || !svertices.Contains(w)) {
					fprintf (stderr, "Edge %d=(%d,%d), membership %d %d.\n", 
						e, v, w, svertices.Contains(v), svertices.Contains(w));
                    fatal ("Inconsistent vertex membership in solution");
                }

                terminals.Remove(v);
                terminals.Remove(w);

                if (uf.Find(v) == uf.Find(w)) {
                    fprintf (stderr, "Vertices %d, %d already in the same component.", v, w);
                    fatal ("Solution has a cycle.");
                } else {
                    uf.Union(v, w);
                    ncomp --;
                }
            }

            if (terminals.Count() > 0) {
                fprintf (stderr, "Terminals not in solution: %d.", terminals.Count());
                fatal ("Missing terminals in solution.");
            }


            fprintf (stderr, "[CHECKING:n=%d:m=%d:%d components]", svertices.Count(), ecount, ncomp);
            if (verbose) {
                for (int v=1; v<=n; v++) {
                    if (svertices.Contains(v)) {
                        fprintf (stderr, "%d:%d ", v, uf.Find(v));
                    }
                }
                fprintf (stderr, "\n");
            }
            
            
            if (ecount != svertices.Count() - 1) return false;
            else return true;
        }


// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF
// NEW STUFF

        /// <summary>
        /// Data structures used for the local search based on key-path exchange.
        /// </summary>
        /// 
        class KeyPathData {
		public:
            VoronoiData voronoi; //main voronoi diagram
            VoronoiData *vorbkp; //backup of the voronoi diagram (to undo changes)
            PairingHeap<EdgeCost> *heap; //heap for original boundary arcs
            //public UniverseSet oldkeynodes; //key nodes in the solution
            UniverseSet *solnodes; //all nodes in the solution
            UniverseSet crucialnodes; //key nodes + terminals
            UniverseSet *removed_arcs; //list of arcs temporarily removed from the solution
            UniverseSet hole; //vertices whose bases belong to the current path
            BinaryHeap<EdgeCost> *binheap; //binary heap for Voronoi computations
            UnionFind uf; //set of already-processed subtrees
            RFWStack<int> *stack; //stack used when marking the hole
            int *bottom; //bottom[v]: bottom vertex on key path containing v

            //data for dfs execution
            DFSData *dfsdata;
            
            //forbidden vertices (those needed for previous paths)
            UniverseSet *pinned;  //vertices that cannot be removed from the solution because they are used in previous moves
            //public bool [] pinned;
            
            UniverseSet *forbdfs; //vertices marked during DFS
            SteinerSolution *newsol;

            KeyPathData(Graph &g) : 
				crucialnodes(g.VertexCount()),
				uf(g.VertexCount()),
				voronoi(g.VertexCount()),
				hole(g.VertexCount())
			{

				const bool verbose = false;
                int n = g.VertexCount();
                int m = g.EdgeCount();

				
                
                //voronoi = new VoronoiData(n);
                vorbkp = new VoronoiData(n);

                //2m elements (arcs) partitioned into n heaps
                heap = new PairingHeap<EdgeCost>(2 * m, n);
                //oldkeynodes = new UniverseSet(n); //key nodes in the solution
                solnodes = new UniverseSet(n); //all nodes in the solution
                //crucialnodes = new UniverseSet(n); //key nodes + terminals
                removed_arcs = new UniverseSet(m); //list of arcs temporarily removed from the solution
                // hole = new UniverseSet(n); //vertices whose bases belong to the current path
                binheap = new BinaryHeap<EdgeCost>(n); //binary heap for Voronoi computations
                stack = new RFWStack<int>(2*m); //stack used for marking the hole and for DFS
                dfsdata = new DFSData(n);

                pinned = new UniverseSet(n);
                forbdfs = new UniverseSet(n);
                newsol = new SteinerSolution(&g);
                bottom = new int[n + 1];
                //uf = new UnionFind(n); //set of already-processed subtrees

				if (verbose) {
					//fprintf (stderr, "Allocated pairing heap %p.\n", heap);
					//fprintf (stderr, "Allocated uf object as %p.\n", uf);
					uf.Find(1);
				}
            }

			~KeyPathData() {
				delete [] bottom;
				delete newsol;
				delete forbdfs;
				delete pinned;

				delete dfsdata;
				delete stack;
				//delete uf;
				delete binheap;
				//delete hole;
				delete removed_arcs;
				//delete crucialnodes;
				delete solnodes;
				delete heap;

				delete vorbkp;
				//delete voronoi;
			}

            void Reset() {
				voronoi.Reset();
                vorbkp->Reset();
                heap->Reset();
                //oldkeynodes.Reset();
                solnodes->Reset();
                //crucialnodes->Reset();
				crucialnodes.Reset();
                removed_arcs->Reset();
                //hole->Reset();
                binheap->Reset();
                uf.Reset();
                stack->reset();
                pinned->Reset();
                forbdfs->Reset();
                newsol->Reset();

                //no need to reset the dfs info---dfs will do it
                //bottom is also reset later...
            }

        };

        //--------------------------------------------
		// Data structure used for key-vertex removal
		//--------------------------------------------
        class KeyVertexData {
		public:
            UniverseSet *cedges; //list of edges (first between components, then on new mst)
            UniverseSet *subgroups; //subgroups that should be avoided
            StaticBuckets *mst; //buckets used for MST computation
            StaticBuckets *hzt; //horizontal arcs
            bool *scanned; //list of scanned vertices
            int *mstparc; //parent in the mst computation

            KeyVertexData(Graph &g) {
                int m = g.EdgeCount();
                int n = g.VertexCount();
                cedges = new UniverseSet(m);
                subgroups = new UniverseSet(n);
                mst = new StaticBuckets(2 * m, n);
                hzt = new StaticBuckets(m, n);
                scanned = new bool [n+1];
                mstparc = new int [n+1];
            }

			~KeyVertexData() {
				delete [] mstparc;
				delete [] scanned;
				delete hzt;
				delete mst;
				delete subgroups;
				delete cedges;
			}
        };


        /// <summary>
        /// Initialize boundary heaps arcs: heap v will contain all the boundary arcs
        /// leaving the voronoi region containing v. Returns a random terminal, as the root.
        /// </summary>
        /// <param name="voronoi">Current voronoi diagram (does not change).</param>
        /// <param name="heap">Pairing heap (will be modified)</param>
        /// <returns></returns>
        static void InitBoundaryHeaps(Graph &g, VoronoiData &voronoi, PairingHeap<EdgeCost> &heap) {
            int bcount = 0;
            bool verbose = false;
            int n = g.VertexCount();
            int m = g.EdgeCount();

			// examine all boundary vertices
            for (int v = n; v > 0; v--) {
                int bv = voronoi.GetBase(v);
                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					//int w = g.GetArcHead(pa); //arc.head;
					int w = a->head;
                    int bw = voronoi.GetBase(w); //neighbor's region
                    if (bv == bw) continue; //inside a boundary region

                    //found a boundary arc
                    bcount++;
                    //int alabel = (v < w) ? arc.label : arc.label + m;
                    //int alabel = (v < w) ? g.GetArcLabel(pa) : g.GetArcLabel(pa) + m;
					int alabel = (v<w) ? a->label : a->label+m;
				
                    EdgeCost cost = voronoi.GetDistance(v) + voronoi.GetDistance(w) + a->cost; //  g.GetArcCost(pa); //arc.cost;
                    //  Console.Error.WriteLine("<({0},{1}):({2},{3}):{4}> ", v, w, bv, bw, cost);
                    heap.Insert(bv, alabel, cost);
                }
            }
            if (verbose) fprintf (stderr, "InitBoundaryHeaps found %d boundary arcs.", bcount);
            //return root;
        }


        static void MergeChildren (Graph &g, int v, KeyPathData &pdata) {
            // fix invariants for next iteration: merge the child paths to k (in uf and heap)
            int *id2parc = pdata.dfsdata->id2parc;
            int *bottom = pdata.bottom;
            UnionFind &uf = pdata.uf;
            PairingHeap<EdgeCost> *heap = pdata.heap;

            //foreach (WeightedGraph.Arc a in g.ArcEnumerator(v)) {
            //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
                //int w = g.GetArcHead(pa); //a.head; 
				int w = a->head;
                if (id2parc[w] == a->label) { //g.GetArcLabel(pa)) { //a.label) { //is w a child of v?
                    int b = bottom[w]; //b is the bottom crucial node of the keypath containing w
                    uf.Union(v,w); 
                    uf.Union(v,b);
                    heap->Merge(v,w);
                    heap->Merge(v,b);
                }
            }
        }






        /// <summary>
        /// Initialize bottom[v] for each vertex v. This is defined as:
        /// - if v is not in the solution, bottom[v] = 0.
        /// - if v is a crucial node, bottom[v]=v.
        /// - if v is noncrucial, bottom[v] is the closest crucial descendent.
        /// The function also groups all internal vertices of a key path into
        /// a single group in both data.uf and data.heap. For the heaps, the
        /// representative will be the top element in the corresponding key
        /// path.
        /// 
        /// The function takes O(n\alpha n) time (with some \alpha).
        /// 
        /// </summary>
        /// <param name="bottom">The array to be filled.</param>
        /// <param name="data">Structure containing the variables kept by the local search.</param>
        static void PrepareKeyPaths (Graph &g, int *bottom, KeyPathData &data) {
            bool verbose = false;
            int n = g.VertexCount();
            int count = 0;
            int ucount = 0;

            UnionFind &uf = data.uf;
            UniverseSet &crucial = data.crucialnodes;
            int *id2parc = data.dfsdata->id2parc;
            PairingHeap<EdgeCost> *heap = data.heap;

            if (verbose) fprintf (stderr, "Marking bottom of each solution vertex now from %d crucial vertices...", crucial.Count());

            //initialization
            uf.Reset();
            for (int v = 1; v <= n; v++) bottom[v] = 0;
            
            //process each key path at once
            //foreach (int b in crucial.ElementEnumerator()) {

			int cp, cpend;
			for (crucial.GetBounds(cp, cpend); cp<cpend; cp++) {
				int b = crucial.PickPos(cp);
                if (verbose) fprintf (stderr, "Processing vertex %d.", b);
                int w = b;
                //move up the key path until reaching another crucial vertex
                while (true) {
                    count ++;
                    bottom[w] = b;
                    int e = id2parc[w];
                    if (e==0) break; //reached the root
                    int p = g.GetOther(e,w);
                    if (crucial.Contains(p)) break; //reached another key path
                    if (w!=b) {
                        ucount++; 
                        uf.Union(w,p);
                        heap->Merge(p,w); //make sure the top element is the heap                        
                    }
                    w = p;
                }
            }
            if (verbose) fprintf (stderr, "enough vertices marked (%d), %d unions performed, %d crucial.\n", count, ucount, crucial.Count());
        }



        /// <summary>
        /// Tarjan's offline nca computation adapted to handle Voronoi regions. Only works if the voronoi
        /// regions are defined by the solution nodes.
        /// </summary>
        /// <param name="hzt">buckets of horizontal edges---will be filled</param>
        /// <param name="voronoi"></param>
        /// <param name="dfsdata"></param>
        /// <param name="root"></param>
        static void ComputeVoronoiNCA(Graph &g, StaticBuckets &hzt, VoronoiData &voronoi, DFSData &dfsdata, int root) {
            bool verbose = false;
            bool debug = false;
            int *id2dfs = dfsdata.id2dfs;
            int *id2parc = dfsdata.id2parc;
            int *dfs2id = dfsdata.dfs2id;

            int v, n = g.VertexCount();
            RFWStack<int> region(n); //used to keep track of the regions visited
            bool *visited = new bool[n + 1]; //nodes already scanned or on the stack (region)
            int *ancestor = new int [n+1];   //auxiliary variable used for nca computations...
            UnionFind ncauf (n); //auxilary data structure for nca computation

            for (v=1; v<=n; v++) {
                visited[v] = false; //nobody visited initially
                ancestor[v] = v;    //and everybody is its own ancestor
            }
            hzt.Reset(); //reset list of visited nodes

            int nhzt = 0; //number of horizontal edges
            int processed = 0; //total number of nodes visited

            //traverse nodes in bottom-up order, ending up in the root
            int count = 0;
            for (int vdfs=1; vdfs<=n; vdfs++) {
                int bv = dfs2id[vdfs];
                if (bv<=0) break;
                if (debug && (bv!=voronoi.GetBase(bv) || visited[bv])) {
					fatal ("invalid violated"); 
				} //throw new Exception ("invariant violated");

                visited[bv] = true;
                region.push(bv);

                //Now we must traverse the adjacency lists of all nodes in the same region as bv.
                while (!region.isEmpty()) {
                    v = region.pop();
                    processed ++;
                    //for (int pa = g.GetStart(v); pa < g.GetEnd(v); pa++) {
					SPGArc *a, *end;
					for (g.GetBounds(v,a,end); a<end; a++) {
						//int w = g.GetArcHead(pa);
						int w = a->head;
                        int bw = voronoi.GetBase(w);

                        //if the other endpoint belongs to the same region, we must make 
                        //sure it is visited; the edge is not on the boundary, though
                        if (bw == bv) {
                            if (!visited[w]) {region.push(w); visited[w] = true;}
                            continue;
                        }

                        //bw in different component; if already processed, it knows the nca
                        if (id2dfs[bw] < vdfs) { 
                            int nca = ancestor[ncauf.Find(bw)]; 
                            if (bv==nca || bw==nca) continue;
                            nhzt++;
                            //hzt.Insert(g.GetArcLabel(pa), nca);
                            hzt.Insert(a->label, nca);
                        }
                    }
                }

                //unite v with its parent
                if (bv!=root) {
                    int parent = g.GetOther(id2parc[bv],bv);
                    ncauf.Union(parent,bv);
                    ancestor[ncauf.Find(bv)] = parent; //parent is the ancestor
                }
            }

            if (processed != n) {
				//fprintf (stderr, ".");
                //fprintf (stderr, "%d/%d nodes processed\n", processed, n);   
                //fprintf (stderr, "WARNING:
				//fatal ("not every node was processed");
            }

            //edges per horizontal edge
            if (verbose) fprintf (stderr, "%d edges are horizontal, total count is %d.", nhzt, count);
        
			delete [] visited;
			delete [] ancestor;
		}

        /// <summary>
        /// Add to data.hole the key path with v at the bottom.
        /// - Does not add v itself.
        /// - Returns the total cost of the keypath.
        /// </summary>
        /// <param name="v">Bottom of the key path.</param>
        /// <param name="data">Current state of computation (uses data.crucialnodes, data.hole(*), and data.dfsdata.id2parc).</param>
        /// <returns>total cost of removed path</returns>
        static EdgeCost AddPathToHole (Graph &g, int v, KeyPathData &data, bool &forbidden) {
            const bool verbose = false;

			ScreenLog (verbose, "Adding path to hole.\n");
			
			int* id2parc = data.dfsdata->id2parc;
            UniverseSet &crucial = data.crucialnodes;
            UniverseSet &hole = data.hole;

            forbidden = false;

			if (verbose) data.uf.Find(1);
			
            if (!crucial.Contains(v)) fatal ("path must begin on a crucial node");
            EdgeCost cost = 0;
            int b = v;
            while (true) {
                int e = id2parc[v];
                cost += g.GetCost(e);
                //Console.Error.WriteLine("[{0}]", e);
                if (v != b) {
                    hole.Insert(v);
                    if (data.pinned->Contains(v)) forbidden = true;
                }
                int p = g.GetOther(e, v);
                if (crucial.Contains(p)) break; //reached to key node
                v = p;
            }
			if (verbose) data.uf.Find(1);
            return cost;
        }

        /// <summary>
        /// Removes all paths P such that top[P]=v. Inserts into vdata.subgroups the 
        /// groups to which the corresponding subtrees belong.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="solution"></param>
        /// <param name="pdata"></param>
        /// <param name="vdata"></param>
        /// <returns></returns>
        static EdgeCost ProcessDownPaths (Graph &g, int v, KeyPathData &pdata, KeyVertexData &vdata, bool &forbidden) {
            int *id2parc = pdata.dfsdata->id2parc;
            int *bottom = pdata.bottom;
            UniverseSet *subgroups = vdata.subgroups;
            UnionFind &uf = pdata.uf;
			const bool verbose = false;

            forbidden = false;
			int ecount = 0;

            EdgeCost downcost = 0;
            //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
            //    int w = arc.head;
            //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				//int w = g.GetArcHead(pa);
				int w = a->head; 

                if (id2parc[w] == a->label) { //g.GetArcLabel(pa)) { //arc.label) {
					if (verbose) ecount ++;
                    int b = bottom[w];
                    subgroups->Insert(uf.Find(b)); 
					EdgeCost curdown = AddPathToHole(g, b, pdata, forbidden);
                    downcost += curdown;
					if (verbose) fprintf (stderr, "Edge %d=(%d,%d) costs %d (downcost is %d), bottom=%d group=%d.\n", a->label, v, w, a->cost, curdown, b, uf.Find(b)); 
                    if (forbidden) break;
                }
            }
			if (verbose) fprintf (stderr, "Found %d down edges costing %d.\n", ecount, downcost);
            return downcost;
        }


        /// <summary>
        /// Adds to data.hole all vertices whose voronoi bases belong to data.hole itself.
        /// (Yes, this is recursive.) Uses data.stack, data.voronoi and data.hole.
        /// </summary>
        /// <param name="data">Current local search state.</param>
        static void ExpandHole(Graph &g, KeyPathData &data) {
            UniverseSet &hole = data.hole;
            VoronoiData &voronoi = data.voronoi;
            RFWStack<int> *stack = data.stack;

            //Console.Error.Write(".");

            stack->reset();
            //foreach (int x in hole.ElementEnumerator()) stack.Push(x);
			int p, pend;
			for (hole.GetBounds(p,pend); p<pend; p++) stack->push(hole.PickPos(p));


            //stack contains hole vertices that have yet to be scanned
            while (!stack->isEmpty()) {
                int v = stack->pop();
                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                //    int w = arc.head; 
                //foreach (int pa in g.ArcPosEnumerator(v)) {
               
				//for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
                //    int w = g.GetArcHead(pa);//arc.head;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
                    if (hole.Contains(w)) continue;
                    if (hole.Contains(voronoi.GetBase(w))) {
                        hole.Insert(w);
                        stack->push(w);
                    }
                }
            }
        }


        /// <summary>
        /// Update the voronoi information associated with all vertices in hole, assuming
        /// all other vertices have the correct information. In particular, assumes no base
        /// belongs to hole.
        /// </summary>
        /// <param name="heap">Auxiliary, pre-allocated binary heap.</param>
        /// <param name="hole">Hole: the vertices that must be assigned to a different region.</param>
        /// <param name="voronoi">The voronoi diagram to be fixed.</param>
        /// <param name="uf">Partition of the original solution.</param>
        /// <param name="subgroup">Set (of the partition) containing a subtree we are interested on.</param>
        /// <param name="bestcost">(Output) cost of the best boundary path found</param>
        /// <returns>Edge representing the best boundary path found (0 if none)</returns>

		/// WARNING! SHOULD THIS BE CALLED? DOESN'T THE OTHER ONE SOLVE THIS?

        static int RepairHole(Graph &g, BinaryHeap<EdgeCost> &heap, UniverseSet &hole, VoronoiData &voronoi, UnionFind &uf, int subgroup, EdgeCost &bestcost, UniverseSet *forbidden) {
            heap.Reset();
            bool verbose = false;

			const bool USE_VORONOI_TIE_BREAKER = false;


            //traverse the adjacency list of each vertex in the hole; update info if we 
            //find a neighbor that does not belong to the hole
            int bcount = 0; //number of boundary vertices
            
            //foreach (int v in hole.ElementEnumerator()) {
			int pe, pend;
			for (hole.GetBounds(pe,pend); pe<pend; pe++) {
				int v = hole.PickPos(pe);

                voronoi.Reset(v);
                int neighbors = 0;
                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                // for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
                //    int w = g.GetArcHead(pa); //arc.head;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;

                    if (hole.Contains(w)) continue;

                    //w is not in the hole: may use edge (v,w) to update priority of v
                    EdgeCost newdist = voronoi.GetDistance(w) + a->cost; //g.GetArcCost(pa); //arc.cost;
                    if ((neighbors == 0) || (newdist < voronoi.GetDistance(v))) {
                        neighbors++;
                        //voronoi.Update(v, voronoi.GetBase(w), arc.label, newdist);
                        voronoi.Update(v, voronoi.GetBase(w), a->label, newdist); //g.GetArcLabel(pa), newdist);
                    }
                }

                //now element will have its best possible value
                if (neighbors > 0) {
                    bcount++;
                    heap.PushBack(v, voronoi.GetDistance(v));
                }
            }
            heap.Heapify();

            if (verbose) fprintf (stderr, "Hole has %d boundary vertices.\n", bcount);

            //prefer terminals
            bool USE_TIEBREAKER = USE_VORONOI_TIE_BREAKER;

            //------
            // We now perform multisource Dijkstra from here.
            //------
            int bestboundary = 0;
            bestcost = 0;
            int count = 0;
            while (!heap.IsEmpty()) {
                unsigned int v;
                EdgeCost dist;
                heap.RemoveFirst(v, dist);
                count++;

                //group to which v's base belongs
                int gbv = uf.Find(voronoi.GetBase(v));

                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
                //    int w = g.GetArcHead(pa); //arc.head;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
                    EdgeCost newdist = dist + a->cost; //g.GetArcCost(pa); //arc.cost;
                    int bw = voronoi.GetBase(w); 
         
                    //do we improve on the neighbor?
                    bool improve = false;
                    if (bw==0 || (newdist < voronoi.GetDistance(w))) {
                        improve = true;
                    } else if (USE_TIEBREAKER && (newdist == voronoi.GetDistance(w)) && hole.Contains(w)) {
                        improve = g.IsTerminal(voronoi.GetBase(v));
                        //if (improve) Console.Error.Write("!");
                    }


                    //if (bw == 0 || (newdist < voronoi.GetDistance(w))) {
                    if (improve) {
                        heap.FixKey(w, newdist); //WARNING: SHOULD BE UPDATE HERE, NO?
                        voronoi.Update(w, voronoi.GetBase(v), a->label, newdist); //g.GetArcLabel(pa), newdist);
                    } else { //if we don't improve the neighbor, this may be a boundary edge
                        int gbw = uf.Find(bw);

                        //a boundary edge links two different groups, one of which
                        //must be the relevant subtree (denoted by subgroup)
                        if ((gbv != gbw) && (gbv == subgroup || gbw == subgroup)) {
                            EdgeCost cost = voronoi.GetDistance(v) + a->cost + voronoi.GetDistance(w); //.GetArcCost(pa) + voronoi.GetDistance(w);
                            if ((bestboundary == 0) || (cost < bestcost)) {
                                //if ((forbidden != null) && (forbidden.Contains(voronoi.GetBase(v)) || forbidden.Contains(voronoi.GetBase(w)))) {
                                    //Console.Error.WriteLine("blah{0} ", g.GetArcLabel(pa));
                                //} else {
                                    bestboundary = a->label; //g.GetArcLabel(pa); //arc.label;
                                   // Console.Error.Write("BB{0}({1},{2}) ", bestboundary, v, w);
                                    bestcost = cost;
                                //}
                            }
                        }
                    }                       
                }
            }
            //Console.Error.WriteLine("Best boundary edge {0} has associated cost {1}.", bestboundary, bestcost);
            if (count != hole.Count()) fatal("invalid hole");

            // do a final check to see if the boundary edge is good enough
            // can't do before because the bases of its endpoints may change during the execution
            // BUT: CAN'T WE HAVE AN INTERMEDIATE FORBIDDEN VERTEX?
            if (bestboundary > 0 && forbidden!=NULL) {
                int x, y;
                g.GetEndpoints(bestboundary, x, y);
                if (forbidden->Contains(voronoi.GetBase(x)) || forbidden->Contains(voronoi.GetBase(y))) {
                    //Console.Error.Write("WEIRD CASE HAPPENED!");
                    //throw new Exception ("Hi");
                    bestboundary = 0;
                }
            }

            return bestboundary;
        }


		/// Repair a hole in the Voronoi diagram
        static int RepairHole(Graph &g, KeyPathData &data, UniverseSet &cedges, UniverseSet &subgroups) {
            BinaryHeap <EdgeCost> *heap = data.binheap;
            UniverseSet &hole = data.hole;
            UnionFind &uf = data.uf;
            VoronoiData &voronoi = data.voronoi;

			const bool USE_VORONOI_TIE_BREAKER = false; //this is what 

            heap->Reset();
            const bool verbose = false;
            bool USE_TIEBREAKER = USE_VORONOI_TIE_BREAKER;


            //traverse the adjacency list of each vertex in the hole; update info if we 
            //find a neighbor that does not belong to the hole
            int bcount = 0; //number of boundary vertices

            //foreach (int v in hole.ElementEnumerator()) {
			int h, hend;
			for (hole.GetBounds(h,hend); h<hend; h++) {
				int v = hole.PickPos(h);
				if (verbose) fprintf (stderr, "h%d:%d ", v, voronoi.GetBase(v));

                voronoi.Reset(v);
                int neighbors = 0;
                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                //    int w = arc.head;
                // for (int pa = g.GetStart(v); pa < g.GetEnd(v); pa++) {
                //    int w = g.GetArcHead(pa);
				SPGArc *a, *end; 
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
                    if (hole.Contains(w)) continue; //skip if neighbor in the hole as well

                    //w is not in the hole: may use edge (w,v) to update priority of v
                    EdgeCost newdist = voronoi.GetDistance(w) + a->cost; //GetArcCost(pa);//arc.cost;
                    if ((neighbors == 0) || (newdist < voronoi.GetDistance(v))) {
                        neighbors++;
                        //voronoi.Update(v, voronoi.GetBase(w), arc.label, newdist);
                        voronoi.Update(v, voronoi.GetBase(w), a->label, newdist); //g.GetArcLabel(pa), newdist);
                    }
                }

                //now element will have its best possible value
                if (neighbors > 0) {
                    bcount++;
                    heap->PushBack(v, voronoi.GetDistance(v));
                }
            }
			
            heap->Heapify(); // prepare initial heap


            if (verbose) fprintf (stderr, "Hole has %d boundary vertices.", bcount);

            //-----------------------------------------------
            // We now perform multisource Dijkstra from here
            //-----------------------------------------------
            int bestboundary = 0;
            //int bestcost = 0;
            int count = 0;
            while (!heap->IsEmpty()) {
                unsigned int v;
                EdgeCost dist;
                heap->RemoveFirst(v, dist); // out dist);
                count++;
                int gbv = uf.Find(voronoi.GetBase(v)); //group containing v's base
				if (verbose) fprintf (stderr, "Scan(v:%d,g:%d,dist:%d)\n", v, gbv, dist); 

                //foreach (WeightedGraph.Arc arc in g.ArcEnumerator(v)) {
                //    int w = arc.head;

                //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
                //    int w = g.GetArcHead(pa);
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
                    EdgeCost newdist = dist + a->cost; //g.GetArcCost(pa);//arc.cost;
                    int bw = voronoi.GetBase(w);

                    //does the neighbor get better?
                    bool improve = false;
                    if (bw == 0 || (newdist < voronoi.GetDistance(w))) {
                        improve = true;
                    } else if (USE_TIEBREAKER) { // prefer terminals   
                        improve = g.IsTerminal(voronoi.GetBase(v)) && 
                            (newdist == voronoi.GetDistance(w)) && 
                            hole.Contains(w);
                    }

                    if (improve) {
                        heap->FixKey(w, newdist); //WARNING: isn't update enough here?
                        //voronoi.Update(w, voronoi.GetBase(v), arc.label, newdist);
                        voronoi.Update(w, voronoi.GetBase(v), a->label, newdist); //g.GetArcLabel(pa), newdist);
                    } else { //if we don't improve the neighbor, this may be a boundary edge
                        int gbw = uf.Find(bw);

                        if (hole.Contains(gbw) || hole.Contains(gbv)) {fatal ("Don't know what I'm doing.");}

                        //Console.Error.Write("<{0},{1},{2},{3}>", gbv, subgroups.Contains(gbv), gbw, subgroups.Contains(gbw));

                        //a boundary edge links two different groups, one of which must not
                        //be the root tree (the root tree is not contained in subgroups)
                        if ((gbv != gbw) && (subgroups.Contains(gbv) || subgroups.Contains(gbw))) {
                            //cedges.Insert(arc.label);
                            cedges.Insert(a->label); //g.GetArcLabel(pa));
							if (verbose) fprintf (stderr, "Boundary edge %d=(%d,%d) with weight %d (%d).\n", a->label, v, w, a->cost, voronoi.GetDistance(v) + a->cost + voronoi.GetDistance(w)); 
                            
                            //Console.Error.Write("RH({0},{1}) ", g.GetFirstEndpoint(arc.label), g.GetSecondEndpoint(arc.label));
                            /* int cost = voronoi.GetDistance(v) + arc.cost + voronoi.GetDistance(w);
                            if ((bestboundary == 0) || (cost < bestcost)) {bestboundary = arc.label; bestcost = cost;}*/
                        }
                    }
                }
            }
            //Console.Error.WriteLine("Best boundary edge {0} has associated cost {1}.", bestboundary, bestcost);
            if (count != hole.Count()) fatal ("invalid hole in voronoi diagram");
            return bestboundary;
        }

        /// <summary>
        /// Clean up heap k until a valid boundary edge is found. We are looking for an edge representing 
        /// an original path with one endpoint in the group containing k, and the endpoint above the hole
        /// in the tree (the subset of solution vertices we are about to remove). While the top heap 
        /// edge does not satisfy these conditions, we remove it from the heap.
        /// 
        /// </summary>
        /// <param name="k">Heap number (the pairing heap we are cleaning).</param>
        /// <param name="subgroups">A valid path may not have its other endpoint in this group (if null, all subgroups are ok)</param>
        /// <param name="data">Current state of the local search: we use uf, voronoi, hole and heap.</param>
        /// <param name="heap_boundary">Output: the remaining top boundary edge (0 is heap empty).</param>
        /// <param name="heap_cost">Output: the cost of the path associated with this boundary edge (0 if heap empty).</param>
        static void CleanUpHeap (Graph &g, int k, UniverseSet *subgroups, KeyPathData &data, int &heap_boundary, EdgeCost &heap_cost, bool discard_forbidden) {
            bool verbose = false;

			ScreenLog (verbose, "Cleaning up heap...\n");
            int m = g.EdgeCount();

            //Console.Error.Write("!");

			ScreenLog (verbose, "About to get structures.\n");

            //get data we have to deal with
            UniverseSet &hole = data.hole;
            VoronoiData *voronoi = data.vorbkp; // WARNING: WHY IS THIS VORBKP?
            UnionFind &uf = data.uf;
			PairingHeap<EdgeCost> *heap = data.heap;
			if (verbose) {uf.Find(1);}

            if (verbose) {
				ScreenLog (verbose, "Should be printing stuff...\n");
                //foreach (int v in hole.ElementEnumerator()) {Console.Error.Write("{0} ", v);}
				int p, pend;
				// print everybody in the hole
				for (hole.GetBounds(p,pend); p<pend; p++) {fprintf (stderr, "%d ", hole.PickPos(p));}

				ScreenLog (verbose, "\nDone with hole.\n");

				if (subgroups) {
	                fprintf (stderr, "Subgroups (%d): ", subgroups->Count());
					for (subgroups->GetBounds(p,pend); p<pend; p++) {
						fprintf (stderr, "%d ", subgroups->PickPos(p));
					}
					fprintf (stderr, "\n");
					fflush (stderr);
				}
			}

			if (verbose) {
				//fprintf (stderr, "Done with verbose stuff (uf=%p, k=%d).\n", uf, k);
				fflush (stderr);
			}
            
            int kgroup = uf.Find(k); //we are looking at outgoing arcs from this group
			if (verbose) {fprintf (stderr, "Group is %d.\n", kgroup); fflush (stderr);}

            int deleted = 0;

			ScreenLog (verbose, "Entering main loop\n");

			verbose = false;
            while (true) {
				if (verbose) {fprintf (stderr, "Checking group %d (heap %p).\n", k, heap);}

				ScreenLog (verbose, "In main loop.\n");


                //empty heap: there's no original boundary edge that qualifies
                if (heap->IsEmpty(k)) {
					ScreenLog (verbose, "Empty!\n");
                    heap_boundary = 0;
                    heap_cost = 0;
                    break;
                }
				ScreenLog (verbose, "Not empty!\n");

                //pick smallest arc on the heap
                heap_boundary = heap->FindMin(k, heap_cost);

                int v, w; //arc is (v,w)
                g.GetDirectedEndpoints(heap_boundary, v, w);


                //we know v belongs kgroup; find the other
                //if (uf.Find(voronoi.GetBase(x))!=kgroup) throw new Exception("Invalid invariant!");
                int other = uf.Find(voronoi->GetBase(w));

                //criteria for valid edges
                bool valid = (other != kgroup); //must link two groups
                valid = valid && !hole.Contains(other); //not incident to the hole
                valid = valid && ((subgroups==NULL) || (!subgroups->Contains(other))); //must not be another subtree

                if (discard_forbidden && valid) {
                    if (data.forbdfs->Contains(data.voronoi.GetBase(w))) valid = false;
                    else if (data.forbdfs->Contains(data.voronoi.GetBase(v))) valid = false;
                }

                //found valid edge? stop
                if (valid) {
                    if (heap_boundary > m) heap_boundary -= m; //replace arc id by edge id                    
                    break;
                }

                //otherwise, discard edge and look for another one
                deleted++;
                heap->DeleteMin(k, heap_cost);
            }

            if (verbose) {
				fprintf (stderr, "Boundary edge %d represents path costing %d. %d edges were deleted.", heap_boundary, heap_cost, deleted);
				fflush (stderr);
			}
        }


        static void FindVerticalEdges(Graph &g, int v, KeyPathData &pdata, KeyVertexData &vdata) {
            //Console.Error.WriteLine("Looking for vertical edges out of {0}.", v);

            int *id2parc = pdata.dfsdata->id2parc;
            UniverseSet *cedges = vdata.cedges;
            int *bottom = pdata.bottom;
            UniverseSet *subgroups = vdata.subgroups;

            //foreach (WeightedGraph.Arc carc in g.ArcEnumerator(v)) {
            //    int w = carc.head;
            //    if (id2parc[w] == carc.label) {
            //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
            //    int w = g.GetArcHead(pa);
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
                //if (id2parc[w] == g.GetArcLabel(pa)) {
				if (id2parc[w] == a->label) {
                    int hboundary;
                    EdgeCost hcost;
                    CleanUpHeap(g, bottom[w], subgroups, pdata, hboundary, hcost, false);
                    if (hboundary > 0) cedges->Insert(hboundary);
                }
            }
        }


        /// <summary>
        /// Add a boundary path to the solution (that's the path between two voronoi 
        /// bases passing through a given boundary edge).
        /// </summary>
        /// <param name="e">The boundary edge that determines the path.</param>
        /// <param name="solution">The solution to which the path will be added.</param>
        /// <param name="voronoi">Current voronoi diagram.</param>
		static void AddBoundaryPath(Graph &g, int e, SteinerSolution &solution, VoronoiData &voronoi) {
            if (solution.Contains(e)) {fprintf(stderr, "o");} //there's an overlap
            //bool verbose = true;

            //add boundary edge first
            solution.Insert(e);
            int x, y;
            g.GetEndpoints(e, x, y);

            //add path to the first base
            while (voronoi.GetParentArc(x) > 0) {
                int f = voronoi.GetParentArc(x);
                solution.Insert(f);
                x = g.GetOther(f, x);
            }

            //add path to the second base
            while (voronoi.GetParentArc(y) > 0) {
                int f = voronoi.GetParentArc(y);
                solution.Insert(f);
                y = g.GetOther(f, y);
            }
        }
		


        /// <summary>
        /// Adds to 'edges' the edges on the boundary path represented by boundary
        /// edge e in 'voronoi'. Traverses the path to do so. If a nonboundary edge
        /// on the path already belongs to 'edges', this function assumes that the
        /// remainder of the path (until the corresponding base) is already there.
        /// </summary>
        /// <param name="e"></param>
        /// <param name="edges"></param>
        /// <param name="voronoi"></param>
        /// <returns></returns>
        static EdgeCost AddBoundaryPath (Graph &g, int e, UniverseSet &edges, VoronoiData &voronoi, UnionFind &uf) {
            EdgeCost icost = 0; //inremental cost
            //solution.Insert(e); //add boundary edge first
            //Console.Error.WriteLine("Trying to add {0}.", e);
            if (edges.Insert(e)) icost += g.GetCost(e);
            int x, y;
            g.GetEndpoints(e, x, y);
            bool verbose = false;
            if (verbose) fprintf (stderr, "BOUNDARY:(%d,%d) B(%d,%d) U(%d,%d)", 
                x, y, voronoi.GetBase(x), voronoi.GetBase(y), uf.Find(voronoi.GetBase(x)), uf.Find(voronoi.GetBase(y)));

            //add path to one of the bases
            while (voronoi.GetParentArc(x) > 0) {
                int f = voronoi.GetParentArc(x);
                if (edges.Insert(f)) icost += g.GetCost(f);
                else break;
                //else {Console.Error.Write("B1"); break;}
                if (verbose) {
                   int a, b;
                   g.GetEndpoints(f, a, b);
                   fprintf (stderr, "F(%d,%d) ", a, b);
                }
                x = g.GetOther(f, x);
            }

            //add path to the other base
            while (voronoi.GetParentArc(y) > 0) {
                int f = voronoi.GetParentArc(y);
                if (edges.Insert(f)) icost += g.GetCost(f);
                else break;
                //else { Console.Error.Write("B2"); break; }
                if (verbose) {
                    int a, b;
                    g.GetEndpoints(f, a, b);
                    fprintf (stderr, "S(%d,%d) ", a, b);
                }
                y = g.GetOther(f, y);
            }

            return icost;
        }



        /// <summary>
        /// Note: vdata.cedges starts with the list of candidate edges and edges with the list of edges 
        /// of the MST.
        /// </summary>
        /// <param name="root"></param>
        /// <param name="data"></param>
        /// <param name="vdata"></param>
        /// <param name="discard_forbidden"></param>
        /// <returns></returns>
        static EdgeCost ReconnectComponents(Graph &g, int root, KeyPathData &data, KeyVertexData &vdata, bool discard_forbidden) { 
            UniverseSet *subgroups = vdata.subgroups;

            bool verbose = false;
            if (verbose) fprintf (stderr, "Reconnecting components.\n");   
            int m = g.EdgeCount();
        
            //if (k == root) throw new Exception ("Should not be trying to remove root.");
            int nchildren = subgroups->Count();
            if (verbose) fprintf (stderr, "We're trying to reconnect %d children; root is %d.", nchildren, root);
            UniverseSet *cedges = vdata.cedges;

            if (nchildren == 1) {
                fprintf (stderr, "ReconnectComponents should not be called with a single child.\n");
                fatal ("Invalid call to ReconnectComponents.");
            }

            StaticBuckets *mst = vdata.mst;
            int *mstparc = vdata.mstparc;
            bool *scanned = vdata.scanned;

            //--------------------------------
            // create list of candidate edges
            //--------------------------------
            int ecount = 0;
            //foreach (int e in cedges.ElementEnumerator()) {
			int pe, pend;
			for (cedges->GetBounds(pe,pend); pe<pend; pe++) {
				int e = cedges->PickPos(pe);


                int x, y; //endpoint of the actual edges
                g.GetEndpoints(e, x, y);

                //find endpoints of this path in the actual tree
                int xtree = data.voronoi.GetBase(x);
                int ytree = data.voronoi.GetBase(y);

                //find subtrees (groups) to which these endpoints begin
                int bx = data.uf.Find(xtree);
                int by = data.uf.Find(ytree);

                if (data.hole.Contains(bx) || data.hole.Contains(by)) fatal ("Bad hole.");

                // If discard_forbidden is true, we must discard all edges incident to a forbidden 
                // region of the current solution. Using such an edge may cause a conflict with
                // previously existing edges.
                if (discard_forbidden && (data.forbdfs->Contains(xtree) || data.forbdfs->Contains(ytree))) {continue;}


                //convention: component containing root is 'root'
                if (!subgroups->Contains(bx)) bx = root;
                if (!subgroups->Contains(by)) by = root;

                //ignore edges between the same components
                if (bx == by) {
                    if (verbose) {
                        fprintf (stderr, "Skip (%d,%d) = (%d,%d) (%d,%d)\n", 
                            x, y, data.voronoi.GetBase(x), data.voronoi.GetBase(y),
                            subgroups->Contains(data.voronoi.GetBase(x)), subgroups->Contains(data.voronoi.GetBase(y)));
                    }
                    continue;
                }

                int xy, yx;
                if (bx > by) {
                    int t = bx; bx = by; by = t; //swap bx and by
                    xy = e+m; yx = e;
                } else {
                    xy = e; yx = e+m;
                }

                ecount += 2;
                
                if (verbose) {
                    fprintf (stderr, "Inserting (%d,%d) into %d and %d (bases:%d %d) %d %d\n", 
                        x, y, bx, by, data.voronoi.GetBase(x), data.voronoi.GetBase(y),
                        subgroups->Contains(data.voronoi.GetBase(x)), subgroups->Contains(data.voronoi.GetBase(y)));
                }

                mst->Insert(xy, bx);
                mst->Insert(yx, by);

            }

            if (verbose) fprintf (stderr, "There are %d candidates to consider.\n", ecount);

            cedges->Reset();       //will use this data structure keep the list of solution edges now
            EdgeCost mstcost = 0;  //cost of the mst, without shortcuts
            EdgeCost realcost = 0; //actual cost, with shortcuts (intersecting paths)
            int nscans = 0;

            //run prim's algorithm, starting from the root component
            BinaryHeap<EdgeCost> *binheap = data.binheap;
            binheap->Reset();
            //binheap.Insert(k, 0);
            binheap->Insert(root, 0);

            //int compcount = nchildren + 1;

            //Console.Error.Write("comp:{0} ", compcount);


            while (true) {
                if (binheap->IsEmpty()) break;

                unsigned int cv;
                EdgeCost curcost;
                binheap->RemoveFirst(cv, curcost);
                
                //unless this is the first arc, we are done
                if (cv != root) {
                    EdgeCost icost = AddBoundaryPath(g, mstparc[cv], *cedges, data.voronoi, data.uf);
                    if (verbose && icost != curcost) fprintf (stderr, "Incremental cost is %d<%d!\n", icost, curcost);
                    realcost += icost;
                }

                scanned[cv] = true;
                mstcost += curcost;

                if (verbose) fprintf (stderr, "\nVertex %d, weight %d, total %d, scans %d, children %d.", cv, curcost, mstcost, nscans, nchildren);
                nscans++;
                if (nscans == nchildren+1) break;

                bool subverbose = false;
                //foreach (int e in mst.ElementEnumerator(cv)) {

				// WARNING: COMPLICATED REPLACEMENT HERE
				for (int e=mst->GetFirst(cv); e>=0; e = mst->GetNext(e)) {
                    int ex, ey;
                    g.GetDirectedEndpoints(e, ex, ey);
                    int bx = data.uf.Find(data.voronoi.GetBase(ex));
                    int by = data.uf.Find(data.voronoi.GetBase(ey));

                    if (subverbose) fprintf (stderr, "\nbx0=%d by0=%d ", bx, by);
                    if (!subgroups->Contains(bx)) bx = root;
                    if (!subgroups->Contains(by)) by = root;
                    if (subverbose) fprintf (stderr, "bx1=%d by1=%d ", bx, by);


                    //we're looking at path (bx,by)
                    if (subgroups->Contains(bx) && bx != cv) {
                        fprintf (stderr, "cv=%d bx=%d by=%d x=%d y=%d basex:%d basey:%d %d\n",
                            cv, bx, by, ex, ey, data.voronoi.GetBase(ex), data.voronoi.GetBase(ey), subgroups->Contains(root));
                        fatal ("Fishy stuff!");
                    }

                    if (scanned[by]) continue;

                    int elabel = e > m ? (e - m) : e;
                    EdgeCost ecost = g.GetCost(elabel) + data.voronoi.GetDistance(ex) + data.voronoi.GetDistance(ey);
                    if (subverbose) fprintf (stderr, "(%d,%d):%d ", bx, by, ecost);

                    //update heap with the new value
                    if (binheap->Insert(by, ecost)) mstparc[by] = elabel;
                }
            }

            bool success = (nscans == nchildren + 1);

            if (verbose) {
                if (success) fprintf (stderr, "Found a solution costing %d (%d).\n", mstcost, realcost);
                else fprintf (stderr, "Could not find a viable spanning tree.\n");
            }

            //----------------------------
            // restore original solutions
            //----------------------------
            mst->Reset(root);
            scanned[root] = false;
            //foreach (int sg in subgroups.ElementEnumerator()) {
			int s, send;
			for (subgroups->GetBounds(s,send); s<send; s++) {
				int sg = subgroups->PickPos(s);
                scanned[sg] = false;
                mstparc[sg] = 0; //this is not really necessary
                mst->Reset(sg);
            }

            //Console.Error.Write("{0} ", success);

            return success ? realcost : -1;
        }


        static EdgeCost DiscardKeyPath (Graph &g, int v, SteinerSolution &solution, KeyPathData &data, bool markforb) {
            int* id2parc = data.dfsdata->id2parc;
            UniverseSet &crucial = data.crucialnodes;
            UniverseSet *forbdfs = data.forbdfs;

            //int remcount = 0;

            if (!crucial.Contains(v)) fatal ("Path must begin on a crucial node.");
            EdgeCost cost = 0;
            int b = v;
            while (true) {
                int e = id2parc[v];
                cost += g.GetCost(e);
                if (!solution.Remove(e)) fatal ("Invalid edge in DiscardKeyPath");
                
                //if (v != b) hole.Insert(v);
                int p = g.GetOther(e, v);

                if (crucial.Contains(p)) break; //reached to key node
                else if (markforb) {forbdfs->Insert(p);}
                v = p;
            }
            //Console.Error.Write("R{0} ", remcount);
            return cost;
        }


		inline static void ScreenLog (const bool verbose, const string &msg) {
			if (verbose) {
				fprintf (stderr, "%s", msg.c_str());
				fflush (stderr);
			}
		}


        /// <summary>
        /// Discards all paths descending from key vertex v
        /// </summary>
        /// <param name="v">key vertex</param>
        /// <param name="solution">current solution</param>
        /// <param name="data">key-path related data</param>
        /// <param name="markforb">forbidden set</param>
        /// <returns>total cost of all discarded paths</returns>        
        static EdgeCost DiscardDownPaths(Graph &g, int v, SteinerSolution &solution, KeyPathData &data, bool markforb) {
            int *id2parc = data.dfsdata->id2parc;
            EdgeCost cost = 0;

            //for (int pa=g.GetStart(v); pa<g.GetEnd(v); pa++) {
            //    int w = g.GetArcHead(pa);
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;

                if (id2parc[w]==a->label) {//g.GetArcLabel(pa)) {
                    cost += DiscardKeyPath (g, data.bottom[w], solution, data, markforb);
                }
            }

            return cost;
        }


        static int KeyVertexElimination (Graph &g, SteinerSolution &solution, RFWLocalRandom &random) {
            //bool perturbed = false; //(random.GetInteger(1,10)==1);
            const bool verbose = false;
			
			KeyPathData data (g);
            KeyVertexData vdata (g);
			if (verbose) {
				fprintf (stderr, "Done initializing data structures!\n");
				fflush(stderr);
			}
			//return 0;
            return KeyVertexElimination(g, solution, data, vdata, true, true, random);
        }



        /// <summary>
        /// Perform key vertex elimination.
        /// </summary>
        /// <param name="solution"></param>
        /// <param name="data"></param>
        /// <returns>Number of moves performed.</returns>
        static int KeyVertexElimination (Graph &g, SteinerSolution &solution, KeyPathData &data, KeyVertexData &vdata, bool allow_paths, bool allow_vertices, RFWLocalRandom &random) {

            if (!allow_paths && !allow_vertices) return 0; //nothing to do here

            //try all possible moves, but don't actually perform any of them;
            //this is useful to measure how many improving moves we could have
            //if there were no conflicts
            bool JUST_COUNT = false;

            // When looking for alternative paths to reconnect the solution, simply
            // ignore paths that are forbidden (i.e., that conflict with previously 
            // executed moves in the same pass). This would cause the algorithm to
            // perform suboptimal exchanges: it may replace an existing key path
            // by another that is not the shortest path between the resulting components
            // (the new path is guaranteed to be shorter than the original key path).
			bool ALLOW_SUBOPTIMAL = true; 


            bool output = false;
            const bool verbose = false;
            bool timer_verbose = true && verbose;

			if (verbose) fprintf (stderr, "Running local search based on key-vertex elimination on solution costing %d.", solution.GetCost());
            int v, n, m;
            n = g.VertexCount();
            m = g.EdgeCount();
            data.Reset();

            int minsubstep = allow_vertices ? 0 : 1;
            int maxsubstep = allow_paths ? 1 : 0;

            RFWTimer timer; //= new OptTimer();
            if (timer_verbose) fprintf(stderr, "PreTime: %.6f\n", timer.getTime());

			if (verbose) fflush (stderr);

            //make sure solution does not have degree-one nonterminals
            Prune(g, solution);

            //mark solution, key, and crucial nodes
			//fprintf (stderr, "Warning: check if dereference this makes sense.\n");
            MarkSolutionNodes(g, solution, *(data.solnodes));
            MarkKeyNodes(g, solution, data.crucialnodes, true);

            if (timer_verbose) fprintf (stderr, "Various marks: %.6f\n", timer.getTime());


            //compute voronoi diagram with all solution vertices as bases
            timer.start();
            ComputeVoronoi(g, data.voronoi, *(data.solnodes), *(data.binheap), NULL);

			if (timer_verbose) fprintf (stderr, "Voronoi time: %.6f\n", timer.getTime());
            data.vorbkp->Clone(data.voronoi);
            if (verbose) fflush (stderr);
            
            //create heap with boundary edges
            int root = PickRandomTerminal(g, random);
			//fprintf (stderr, "%d ", root);

            //initialize boundary heaps (WARNING: WOULD LIKE TO ELIMINATE INTRA-PATH EDGES FROM THE START...
            InitBoundaryHeaps(g, data.voronoi, *(data.heap)); //, data.uf);
       
            if (timer_verbose) fprintf (stderr, "heap: %.3f\n", timer.getTime());
            if (verbose) fflush (stderr);


 
            //perform DFS to sort the vertices
            int visited = DFS(g, root, solution, *(data.dfsdata), *(data.stack));
            //return 0;
			ScreenLog (verbose, "Done with initial DFS.\n");

			int *bottom = data.bottom; //bottom[v]: lowest vertex on key path containing v 
            //KeyVertexData vdata = null;
            UniverseSet *cedges = NULL; //list of edges (first boundary, on the new mst)
            UniverseSet *subgroups = NULL; //representatives for each subtree


            if (allow_vertices) {
                //vdata = new KeyVertexData(g);
                cedges = vdata.cedges; 
                subgroups = vdata.subgroups;
                for (v=0; v<=n; v++) {
                    vdata.scanned[v] = false;
                    vdata.mstparc[v] = 0;
                }
            }


            //mark bottom
            PrepareKeyPaths(g, bottom, data); //data.uf, data.crucialnodes, data.dfsdata);
			ScreenLog (verbose, "Done preparing key paths.\n");


			/*			
            if (output) {
                g.OutputGV("bottom.gv", solution, bottom);
                if (output) Console.Error.WriteLine("Solution costs {0}", solution.GetCost());
            }*/

            
            // if we plan to do vertex deletion, compute horizontal edges
            if (allow_vertices) {
                //ComputeVoronoiNCA(vdata.hzt, data.voronoi, data.dfsdata, root);
                ComputeVoronoiNCA(g, *(vdata.hzt), data.voronoi, *(data.dfsdata), root);
            }

            //return 0;
            if (timer_verbose) fprintf(stderr, "InitTime: %.6f\n", timer.getTime());
			ScreenLog (verbose, "Done computing VoronoiNCA\n");


            //--------------------------------
            // Process vertices in post-order
            //--------------------------------
            int ccount = 0;
            int totalcost = 0;
            int solsize = data.solnodes->Count();
            int impcount = 0;
            EdgeCost imptotal = 0;

            bool impverbose = true && verbose;

			if (verbose) fprintf (stderr, "----\n");

            for (int i = 1; i < solsize; i++) {
                //bool conflict = false; //is there a conflict with a previous swap?
				int k = data.dfsdata->dfs2id[i];
                if (!data.crucialnodes.Contains(k)) continue;
                ccount ++;

                if (verbose) {
					fprintf (stderr, "\nProcessing vertex %d (dfs=%d) (ccount=%d) terminal=%d.\n", k, i, ccount, g.IsTerminal(k));
					fflush (stderr);
					data.uf.Find(1);
				}

				// try to eliminate key vertex first, then parent key path
                for (int substep=minsubstep; substep<=maxsubstep; substep++) {
					if (impverbose) {fprintf (stderr, "[%d:%d] ", k, substep); fflush (stderr);}
                    EdgeCost oldcost, newcost;
                    oldcost = newcost = 0;
                    EdgeCost improvement = 0;

                    if (substep == 0) { //try to eliminate key vertex
                        //Console.Error.Write ("S0 ");
                        if (!g.IsTerminal(k) && !data.pinned->Contains(k)) {
                            //Console.Error.WriteLine("Trying to remove {0}", k);

                            data.hole.Reset();
                            subgroups->Reset();

							if (impverbose) fprintf (stderr, "Trying %d.\n", k);

                            //WARNING: SHOULDN'T I INSERT K AND THE UPPER PATH INTO THE SUBGROUPS? WHAT IS SUBGROUPS?
                            bool forbidden;
                            EdgeCost upcost = 0;
                            EdgeCost downcost = ProcessDownPaths(g, k, data, vdata, forbidden); //remove paths below k, initialize vdata.subgroups

                            if (!forbidden) {
                                //throw new Exception("WARNING: NOT PROCESSING FORBIDDEN PROPERLY");                           
                                data.hole.Insert(k); //remove k itself
                                upcost = AddPathToHole(g, k, data, forbidden); //remove the path whose bottom is k
                                oldcost = downcost + upcost;
								if (impverbose) fprintf (stderr, "upcost is %d.\n", upcost);
                            }

                            if (!forbidden) {
                                if (verbose) fprintf (stderr, "down(%d) + up(%d) = total(%d) [%d] hole:%d", downcost, upcost, oldcost, totalcost, data.hole.Count());
                                ExpandHole(g, data); //expand hole to include all vertices currently assigned to it

                                //insert into cedges all boundary edges that will be used in the computation
                                cedges->Reset();
                                RepairHole(g, data, *cedges, *subgroups); //new boundary edges
                                if (impverbose) fprintf (stderr, "There are %d candidate edges after RepairHole.\n", cedges->Count());
								FindVerticalEdges(g, k, data, vdata);   //original vertical edges
                                if (impverbose) fprintf (stderr, "There are %d candidate edges after FindVerticalEdges.\n", cedges->Count());
                               
								//foreach (int ha in vdata.hzt.ElementEnumerator(k)) {cedges.Insert(ha);} //original horizontal edges
								
								// get all edges in bucket k
								if (verbose) fprintf (stderr, "n"); //this is new...
								for (int e=vdata.hzt->GetFirst(k); e>=0; e = vdata.hzt->GetNext(e)) {
									cedges->Insert(e);
								}
								
		
                                //use cedges to reconnect the various components
                                //Console.Error.WriteLine("Solbefore: {0} cedges:{1}", solution.GetCost(), cedges.Count());
                                newcost = ReconnectComponents(g, root, data, vdata, ALLOW_SUBOPTIMAL);
								if (impverbose) fprintf (stderr, "New cost is %d.\n", newcost);

                                //maybe we couldn't find a valid solution
                                if (newcost < 0) { //could not find a valid solution
                                    if (!ALLOW_SUBOPTIMAL) fatal ("There should always be a valid solution.");
                                    forbidden = true; 
                                } else {
                                    //foreach (int f in cedges.ElementEnumerator()) {
									int pf, pfend;
									for (cedges->GetBounds(pf,pfend); pf<pfend; pf++) {
										int f = cedges->PickPos(pf);
                                        int x, y;
                                        g.GetEndpoints(f, x, y);
                                        if (data.forbdfs->Contains(x) || data.forbdfs->Contains(y)) {
                                            if (ALLOW_SUBOPTIMAL) fatal ("Should not have generated a tree with forbidden endpoints.");
                                            forbidden = true;
                                            break;
                                        }
                                    }
								}

								if (verbose) fprintf (stderr, "old=%d new=%d\n", oldcost, newcost);
						
                                
                                improvement = oldcost - newcost;
								if (!forbidden && improvement > EDGE_COST_PRECISION) {
                                    impcount++;
                                    if (!JUST_COUNT) {
                                        //Console.Error.Write("R0 ");
                                        EdgeCost discarded = 0;
                                        discarded += DiscardKeyPath(g, k, solution, data, true); //remove path with bottom k from solution, mark it in forbdfs
                                        data.stack->reset();
                                        data.stack->push(k);
                                        MarkComponent(g, *(data.stack), data.dfsdata->id2parc, *(data.forbdfs));

                                        //remove the path from the solution
                                        discarded += DiscardDownPaths(g, k, solution, data, false);

										EdgeCost added = 0;
                                        //foreach (int e in cedges.ElementEnumerator()) {
										int pe, pend;
										for (cedges->GetBounds(pe,pend); pe<pend; pe++) {
											int e = cedges->PickPos(pe);


                                            added += g.GetCost(e);
                                            solution.Insert(e);
                                            int x, y;
                                            g.GetEndpoints(e, x, y);
                                            //Console.Error.Write("+({0},{1})", x, y);
                                            data.pinned->Insert(x);
                                            data.pinned->Insert(y);

                                        }
                                    }
                                }

                                //restore voronoi diagram for next stage
                                //foreach (int h in data.hole.ElementEnumerator()) 
								int hp, hpend;
								for (data.hole.GetBounds(hp,hpend); hp<hpend; hp++) {
									int h = data.hole.PickPos(hp);
									data.voronoi.CopyFrom(h, *(data.vorbkp));
								}
                            }
                        }
                        // fix invariants for next iterations: merge the child paths to k (in uf and heap)
                        MergeChildren(g, k, data);
                    } else { //try to replace key path

                        //Console.Error.Write("S1 ");
                        const bool rpverbose = verbose;

                        if (rpverbose) {fprintf (stderr, "\nPROCESSING %d ", k); data.uf.Find(1);}

                        if (data.uf.Size(k) == 1) MergeChildren(g, k, data);
                        data.hole.Reset();

						ScreenLog (rpverbose, "MergeChildren\n");
						if (rpverbose) {data.uf.Find(1);}

                        bool double_check = false && !ALLOW_SUBOPTIMAL; 
                        EdgeCost bestcost = 0;
                        int bestedge = 0;
                        if (cedges != NULL) cedges->Reset();
						if (rpverbose) {data.uf.Find(1);}

                        //mark hole
                        bool conflict;
                        oldcost = AddPathToHole(g, k, data, conflict); //mark hole

						ScreenLog (rpverbose, "Added path to hole\n");
						if (rpverbose) {data.uf.Find(1);}

                        //if any of vertex on this path is pinned, we can't remove it
                        if (!conflict) {
                            //if path has internal vertices must expand and repair hole
                            if (data.hole.Count()>0) {
                                ExpandHole(g,data);
                                UniverseSet *flist = ALLOW_SUBOPTIMAL ? data.forbdfs : NULL;
                                bestedge = RepairHole (g, *(data.binheap), data.hole, data.voronoi, data.uf, data.uf.Find(k), bestcost, flist);
                                if (cedges!=NULL && bestedge > 0) cedges->Insert(bestedge); //are we using 
                            }

                            
                            //find best original vertical edge
                            if (rpverbose) {
                                if (cedges!=NULL && subgroups!=NULL) {
                                    fprintf (stderr, "Repair:%d %d", cedges->Count(), subgroups->Count());
                                }
                            }
							ScreenLog (rpverbose, "Done repairing\n");

                            int hboundary;
                            EdgeCost hcost;
                            CleanUpHeap(g, k, NULL, data, hboundary, hcost, ALLOW_SUBOPTIMAL);


                            if (cedges!=NULL && hboundary > 0) cedges->Insert(hboundary);

                            //pick the best among those two
                            if (rpverbose) fprintf (stderr, "k:%d hole:%d:%d heap:%d:%d ", k, bestedge, bestcost, hboundary, hcost);
							ScreenLog (rpverbose, "Heap cleaned.\n");

                            //Console.Error.Write("({0},{1})", hboundary, bestedge);

                            if ((hboundary>0) && (bestedge==0 || bestcost>hcost)) {
                                //Console.Error.Write(" UPDATED ");
                                bestedge = hboundary;
                                bestcost = hcost;
                            }

                            
                            if (double_check) {
                                newcost = ReconnectComponents(g, root, data, vdata, false);
                                if (newcost != bestcost) fprintf (stderr, "Costs differ!\n");
                            } else newcost = bestcost;

                            //WARNING: MUST ACTUALLY CHANGE THE SOLUTION

                            if (bestedge == 0) {
                                //this can happen if we eliminated forbidden edges
                                //Console.Error.Write("n");
                                conflict = true;
                            }

                            if (rpverbose) fprintf (stderr, "-> %d:%d %d", bestedge, bestcost, newcost);
                            improvement = oldcost - newcost;
							if (!conflict && (improvement > EDGE_COST_PRECISION)) {

                                //found an improvement: will change the solution
                                int x, y, e = bestedge;
                                g.GetEndpoints(e, x, y);
                                int bx = data.voronoi.GetBase(x);
                                int by = data.voronoi.GetBase(y);

                                //the edge we are trying to add is (bx,by)

                                conflict = (data.forbdfs->Contains(bx) || data.forbdfs->Contains(by));
                                if (conflict && ALLOW_SUBOPTIMAL) {
                                    fprintf (stderr, "Conflicts should have been previously detected for edge...\n"); // " + bestedge + "(" + x + "," + y + ").");
									//throw new Exception ("Conflicts should have been previously detected for edge " + bestedge + "(" + x + "," + y + ").");
                                }

                                if (!conflict) {
                                    impcount ++;
                                    if (!JUST_COUNT) {
                                        //neither endpoint can be removed
                                        data.pinned->Insert(bx); 
                                        data.pinned->Insert(by); 

                                         //remove path with bottom k from solution, mark it in forbdfs
                                        DiscardKeyPath(g, k, solution, data, true);
                                        data.stack->reset();
                                        data.stack->push(k);
                                        MarkComponent(g, *(data.stack), data.dfsdata->id2parc, *(data.forbdfs));
                                        AddBoundaryPath(g, bestedge, solution, data.voronoi);
                                    }
                                }
                            }


                            //restore the original voronoi diagram
                            //foreach (int h in data.hole.ElementEnumerator()) data.voronoi.CopyFrom(h, data.vorbkp);
							
							// restore the original Voronoi diagram
							int pe, pend;
							for (data.hole.GetBounds(pe,pend); pe<pend; pe++) {
								int h = data.hole.PickPos(pe);
								data.voronoi.CopyFrom(h, *data.vorbkp);
							}


						}
					}

                    improvement = oldcost - newcost;
                    if (improvement > EDGE_COST_PRECISION) {
                        //bool impverbose = true;
                        if (impverbose) fprintf (stderr, "%d:%d (%d->%d) ", k, improvement, oldcost, newcost);
                        //impcount ++;
                        imptotal += improvement;
                    }
                }
            }



            const bool final_verbose = false;
            if (final_verbose) fprintf (stderr, "There are %d improving moves with total improvement %d, solution costs %d.\n", impcount, imptotal, solution.GetCost());
            if (timer_verbose) fprintf (stderr, "Final time: %.6f.\n", timer.getTime());
            //moves = impcount;

			//fprintf (stderr, "%d ", impcount);
            return impcount; //number of improving moves
        }








};
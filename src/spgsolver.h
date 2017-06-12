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
#include "constructive.h"
#include "execution_log.h"
#include "LSVertexInsertion.h"
#include "LSVertexElimination.h"
#include "LSBasics.h"
#include "LSKeyPath.h"
#include "BranchBound.h"
#include "perturbation.h"
#include "preprocess.h"

using namespace std;



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
		Basics::ReportResults(stdout, "ms", mstime, mscost, bestknown); 

		fprintf (stdout, "msiterations %d\n", msit);
		fprintf(stdout, "mstype %d\n", mstype);
		//fprintf (stdout, "mssolution %.0f\n", mscost);
		//fprintf (stdout, "mstimeseconds %.12f\n", mstime);
		gbestfound = bestfound;
	}

	static void Solve (int argc, char **argv) {
		const uint64_t version = 201706010852;
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

		MULTISTART = (msit != 0);

		RFWRandom::randomize(seed);


		bool APPLY_PERTURBATION = false;
		if (APPLY_PERTURBATION) {
			fprintf (stderr, "Applying perturbation... ");
			int m = g.EdgeCount();
			vector<EdgeCost> pertcost (m+1,-1);
			RFWLocalRandom random(seed+17);
			PerturbationTools::ApplyPerturbation(g, pertcost, random, 1, 1.0001);
			g.ApplyCosts(pertcost);
			for (int i=1; i<std::min(10, m); i++) {
				fprintf (stderr, "%.5f ", (double)g.GetCost(i));
			}
			fprintf (stderr, "done.\n");
		}


		if (PREPROCESS) {
			fprintf (stderr, "Should be preprocessing.\n");
			RFWTimer preptime(true);
			Preprocessing::RunPreprocessing(g, !SAFE_PREPROCESS);
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
					BranchBound::RunBranchAndBound(thread_g, bbseed, thread_bestfound, thread_bestfound, bestknown, config, &ginfo, &executionLog);
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
						case 2: ConstructiveAlgorithms::SPH (g, solution, NULL, 1); break;
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
			BranchBound::RunBranchAndBound(g, seed, primal, bestfound, bestknown, config, NULL, &executionLog);
		}
		double walltime = timer.getTime();
		fprintf (stdout, "totalwalltimeseconds %.12f\n", walltime);
		//fprintf (stdout, "totaltimeseconds %.12f\n", walltime + second_time);
		//fprintf (stdout, "bestsolution %.0f\n", bestfound);

		Basics::ReportResults(stdout, "total", walltime + first_time, bestfound, bestknown);
		Basics::ReportResults(stdout, "totalcpu", walltime + second_time, bestfound, bestknown);

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

		int root = Basics::PickRandomTerminal(g, random);
		ConstructiveAlgorithms::SPH (g, target, &pertcost[0], root);
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
				PerturbationTools::InitPerturbation(g, pertcost, random, config);
			}
			ConstructiveAlgorithms::SPH (g, cursol, USE_PERTURBATION ? &pertcost[0] : NULL, root);
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

		if (verbose) fprintf (stderr, "%d->", solution.GetCost());
		while (failures_to_go > 0) {
			SteinerSolution *refsol = elite.GetReference(random.GetInteger(1, elite.Count()));
			CombineSolutions(combsol, solution, *refsol, random, config);
			if (!combsol.IsBetter(&solution)) {
				failures_to_go --;
			} else {
				solution.CopyFrom(&combsol);
			}
		}
		if (verbose) fprintf (stderr, "%d\n", solution.GetCost());
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
				ConstructiveAlgorithms::SPH(g, tempsol, NULL, v);
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
		const bool verbose = true;

		double precision = EDGE_COST_PRECISION;

		// run a few iterations of the local search while decaying the perturbation (towards unperturbed solution)
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

		g.ApplyCosts(original); //restore original edge costs
		solution.UpdateCost();  //recost the existing solution
		EdgeCost before = solution.GetCost(); // this is the original cost
		//fprintf (stderr, "d%.0f->", solution.GetCost());
		//SPH (g, solution, NULL, root);
		//fprintf (stderr, "[%d->", solution.GetCost());
		
		// run local search on current solution using original costs 
		MSTPrune(g,solution);
		RunLocalSearch(g, solution, random, -1, config);		
		EdgeCost after = solution.GetCost();
		if (verbose) {
			fprintf (stderr, " %.0f", after);
			//fprintf (stderr, " %3.0f", before - after);
		}
		//fprintf (stderr, "%d]", solution.GetCost());
	}

	/**
	* Generate randomized solution from scratch by perturbing edge weights, runing a constructive algorithm, then applying local search. 
	* The local search is applied to the perturbed graph initially, but the perturbation may be gradually dampened (removed). 
	*
	*/
	static void GenerateRandomizedSolution(SteinerSolution &solution, int root, vector<EdgeCost> &pertcost, RFWLocalRandom &random, SteinerConfig *config, ExecutionLog *executionLogPtr = nullptr) {
		Graph &g = *solution.g;
		int m = g.EdgeCount();
		const bool VERBOSE_STEP = false;
				
		vector<EdgeCost> original (m+1);
		//remember original costs
		g.RetrieveCosts(original);

		// find local optimum with perturbed costs
		g.ApplyCosts(pertcost);
		ConstructiveAlgorithms::SPH (g, solution, NULL, root);
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
		if (config) {
			int confpert = config->RESILIENT_PERTURBATION;
			if (confpert == 0) {RESILIENT_PERTURBATION = false;}
			else if (confpert == 1) {RESILIENT_PERTURBATION = true;}
		}
		
		fprintf (stderr, "Using resilient perturbation? %d\n", RESILIENT_PERTURBATION);
		
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
					PerturbationTools::AdaptivePerturbation(g, pertcost, elite, random);
				} else {
					//fprintf (stderr, "Here!\n");
					PerturbationTools::InitPerturbation(g, pertcost, random, config);
					//fflush (stderr);
				}
			}

			if (RESILIENT_PERTURBATION) {
				// both constructive and the local search use perturbation
				GenerateRandomizedSolution (solution, root, pertcost, random, config);
			} else {
				// non-resilient perturbation: 
				ConstructiveAlgorithms::SPH (g, solution, USE_PERTURBATION ? &pertcost[0] : NULL, root);
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
			int root = Basics::PickRandomTerminal(g, random);
			if (USE_PERTURBATION) PerturbationTools::InitPerturbation(g, pertcost, random, NULL);
			ConstructiveAlgorithms::SPH (g, solution, USE_PERTURBATION ? &pertcost[0] : NULL, root);
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
				LSVertexInsertion::VertexInsertion(g, solution, n, random);
				MSTPrune(g, solution);
				if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
			}
			//fprintf (stderr, "Starting key vertex elimination!\n");
			//fflush (stderr);
			
			if (RUN_Q) {
				LSKeyPath::KeyVertexElimination(g, solution, random);
				//fprintf (stderr, "Ending key vertex elimination!\n");
				if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
			} 
			
			if (RUN_U) {
				LSVertexElimination::VertexElimination(g, solution, random);
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
		
		ConstructiveAlgorithms::SPH (g, solution, NULL, Basics::PickRandomTerminal(g,random));
		if (verbose) fprintf (stderr, "SPH found solution %d.\n", solution.GetCost());
		int n = g.VertexCount();
		EdgeCost oldcost = solution.GetCost();
		RFWTimer timer(true);
		int rounds = 0;
		for (int i=0; i<10; i++) {
			rounds ++;
			if (RUN_V) {
				LSVertexInsertion::VertexInsertion(g, solution, n, random);
				MSTPrune(g, solution);
				if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
			}
			//fprintf (stderr, "Starting key vertex elimination!\n");
			//fflush (stderr);
			
			if (RUN_Q) {
				LSKeyPath::KeyVertexElimination(g, solution, random);
				//fprintf (stderr, "Ending key vertex elimination!\n");
				if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
			} 
			
			if (RUN_U) {
				LSVertexElimination::VertexElimination(g, solution, random);
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





	// Comput minimum spanning tree of the graph g using Prim's algorithm (ignores terminals)
	// 'solution' will contain the MST edges

	static void MSTPrim (Graph &g, SteinerSolution &solution) {
		bool verbose = false;
		int n = g.VertexCount();
		BinaryHeap<EdgeCost> heap(n); // = new BinaryHeap<ArcCost>(n);
		vector<int> parc (n+1); //not need to initialize
		unsigned int r = Basics::PickRandomTerminal(g);
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
        Basics::ComputeVoronoi(g, voronoi, baselist, heap, NULL); //pertcost);


		solution.Reset();
//        return; 
		//fprintf (stderr, "Missing boruvka!\n");
        //Boruvka(g, solution, voronoi, uf, NULL);
        Preprocessing::BoruvkaGraph(g, solution, voronoi, uf, NULL);
		//DNHPrim(g,solution,voronoi,uf);
            //Console.Error.WriteLine("t3:{0} ", timer.GetTime());
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
            Basics::MarkSolutionNodes(g, solution, svertices);
            MST(g,solution,svertices);
            Basics::Prune(g,solution);
			//fprintf (stderr, "Solution costs %d, gain is %d.\n", solution.GetCost(), original - solution.GetCost());
            //Prune(g,solution);
			//fprintf (stderr, "g%d ", original - solution.GetCost());
            return (solution.GetCost() < original) ? 1 : 0;
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
            int r = Basics::PickRandomTerminal(g);
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




};


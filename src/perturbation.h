#pragma once
#include "spgconfig.h"

	
class Perturbator {
	int PERTURBATION_MODE;
	double base;
	double range;
	double threshold;
	SteinerConfig *config;

	static void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
	}

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
			//fprintf(stderr, "r is %.10f, threshold=%.10f\n", r, threshold);
			//cerr << r << endl;
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



class PerturbationTools {
public:
	static void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
	}
	
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
			//fprintf(stderr, "<<<<");
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
		//fprintf (stderr, "USING VERTEX PERTURBATION? %d\n", USE_VERTEX_PERTURBATION);

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

};


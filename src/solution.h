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
#include <cassert>
#include "graph.h"

class SteinerSolution {
public:
	Graph *g;
private:
	vector<bool> edge;  //incident vector of edges in the solution
	vector<int> degree; //current degrees of all solution vertices
	int nleaves;        //number of nonterminal leaves
	EdgeCost cost;      //total solution cost

	inline void ProcessInsertion(int v) {
		int d = ++degree[v];
        if (d<=2 && !g->IsTerminal(v)) {
			if (d==1) nleaves ++;
            else nleaves --;
		}
	}

		// process v when removing one of its endpoints
	inline void ProcessRemoval (int v) {
		int d = --degree[v];
        if (d<2 && !g->IsTerminal(v)) {
			if (d==0) nleaves --;
            else nleaves ++;
		}
	}


public:
	inline int LeafCount() const {return nleaves;}
	inline bool Contains(int e) const {  return edge[e]; }
	inline size_t EdgeCapacity() const { return edge.size(); }

	void Reset() {
		fill_n(edge.begin(), g->EdgeCount()+1, false);
		fill_n(degree.begin(), g->VertexCount()+1, 0);
		cost = g->GetFixedCost();
		nleaves = 0;
	}

	inline bool Insert (int e) {
		if (Contains(e)) return false;
		edge[e] = true;
		cost += g->GetCost(e);
		ProcessInsertion(g->GetFirstEndpoint(e));
		ProcessInsertion(g->GetSecondEndpoint(e));
		return true;
	}

	int EdgeCount() {
		int ecount = 0;
		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			if (Contains(e)) ecount ++;
		}
		return ecount;
	}

        void Output(FILE *file) {
            fprintf (stderr, "Outputting solution..."); // to {0}... ", filename);
            //System.IO.StreamWriter file = new System.IO.StreamWriter(@filename);


			int m = g->EdgeCount();
            fprintf (file, "m %d\n", EdgeCount());
            fprintf (file, "c %.0f\n", GetCost());
            for (int e = 1; e<=m; e++) {
				if (!Contains(e)) continue;
                int v, w;
                g->GetEndpoints(e, v, w);
                fprintf (file, "e %d %d %.0f\n", v, w, g->GetCost(e));
            }
            fprintf (stderr, "done.\n");
        }

	void Output (char *prefix) {
		char filename[2048];
		sprintf (filename, "%s.%09.0lf.sol", prefix, (double)GetCost());
		FILE *file = fopen (filename, "w");
		if (!file) {
			fprintf (stderr, "Could not open <%s> for output. Will not output solution.\n", filename);
			return;
		}
		Output(file);
		fclose(file);
		fprintf (stderr, "Solution written to <%s>.\n", filename);
		fflush(stderr);
	}




	/// Remove edge e from the current solution. Returns true iff
    /// the edge was actually in the solution before. 
    inline bool Remove(int e) {
		if (edge[e]==false) return false;
		edge[e] = false;
		cost -= g->GetCost(e);
		ProcessRemoval(g->GetFirstEndpoint(e));
		ProcessRemoval(g->GetSecondEndpoint(e));
		return true;
	}


	inline int GetDegree (int v) const {return degree[v];}
	inline EdgeCost GetCost () const {return cost;}

	bool IsBetter (SteinerSolution *s) {
		return (GetCost() < s->GetCost() - EDGE_COST_PRECISION);
	}

	void UpdateCost() {
		const bool verbose = false;
		if (verbose) fprintf (stderr, "Cost updated from %d to ", cost);
		int m = g->EdgeCount();
		cost = g->GetFixedCost();
		for (int e=1; e<=m; e++) {
			if (Contains(e)) cost += g->GetCost(e);
		}
		if (verbose) fprintf (stderr, "%d.\n", cost);
	}



	SteinerSolution (Graph *_g) {
		g = _g;
		edge.resize(g->EdgeCount()+1);
		degree.resize(g->VertexCount()+1);
		Reset();
	}

	SteinerSolution (SteinerSolution *s) {
		g = s->g;
		edge.resize(g->EdgeCount()+1);
		degree.resize(g->VertexCount()+1);
		Reset();
		CopyFrom(s);
	}

	void CopyFrom(SteinerSolution *s) {
		if (g != s->g) {
			g = s->g;
		} else {
			Reset();
		}
		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			if (s->Contains(e)) Insert(e);
		}

		//fprintf (stderr, "Created solution with %d arcs.\n", this->EdgeCount());

        //foreach (int e in ((SteinerSolution)s).ElementEnumerator()) {
        //        Insert(e);
        //}
            //Console.Error.WriteLine("copy not implemented");   
	}


	/// <summary>
	/// Compute degree of difference between this solution and ds, given 
	/// by (union - intersection) / intersection. In particular, the result
	/// is zero iff the solutions are identical, and 1 if they have no edge
    /// in common.
	double GetDifference(SteinerSolution *s) {
            //SteinerSolution s = (SteinerSolution)ds;
            //int common = 0;
            //int scount = 0;

		int ucount = 0; //edges.Count(); //number of edges in the union
        int icount = 0; //number of edges in the intersection

		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			int count = (Contains(e) ? 1 : 0) + (s->Contains(e) ? 1 : 0);
			if (count == 0) continue;
			//fprintf (stderr, "%d%d ", Contains(e), s->Contains(e));
			ucount ++; //1 or 2
			if (count == 2) {
				icount ++;
			}
		}

		//fprintf (stderr, "[%d/%d] ", icount, ucount);

        return (1.0 - (double)icount/(double)ucount);
	}

};

#if 0
using ArcCost = System.Double;

namespace SteinerOpt {
    [Serializable]
    public class SteinerSolution : IOptSolution {
        private UniverseSet edges; //list of edges in the solution
        private ArcCost cost;
        private int [] degree;
        private int nleaves; //number of degree-one nonterminals

        [NonSerialized]
        WeightedGraph g;

        /// <summary>
        /// Compute degree of difference between this solution and ds, given 
        /// by (union - intersection) / intersection. In particular, the result
        /// is zero iff the solutions are identical, and 1 if they have no edge
        /// in common.
        /// </summary>
        /// <param name="ds"></param>
        /// <returns></returns>
        public double GetDifference(IOptSolution ds) {
            SteinerSolution s = (SteinerSolution)ds;
            //int common = 0;
            //int scount = 0;

            int ucount = edges.Count(); //number of edges in the union
            int icount = 0; //number of edges in the intersection

            foreach (int e in s.ElementEnumerator()) {
                if (edges.Contains(e)) {icount ++;}
                else ucount ++;
            }

            return (1.0 - (double)icount/(double)ucount);
              

            /*

            foreach (int e in s.ElementEnumerator()) {
                scount++;
                if (s.Contains(e)) common++;
            }

            int total = edges.Count() + scount;
            return (double)(total - common) / (double)total;
             */
        }

        public void OutputSolution(string filename) {
            Console.Error.Write("Outputting solution to {0}... ", filename);
            System.IO.StreamWriter file = new System.IO.StreamWriter(@filename);
            /*file.WriteLine("c <id> <inbest> <degreebest> <frac of sols as key vertex>");
            file.WriteLine("n {0}", n);
            file.WriteLine("b {0}", bestsol.GetCost());
            file.WriteLine("s {0}", nsols);
            for (v=1; v<=n; v++) {
                int degree = bestsol.GetDegree(v);
                int inbest = 0;
                if (degree > 0) inbest = 1;
                //if (bestsol.Contains(v)) {inbest = 1;}
                file.WriteLine("v {0} {1} {2} {3}", v, inbest, degree, (double)vcount[v]/(double)nsols);
            }
            file.Close();
            */

            file.WriteLine("m {0}", edges.Count());
            file.WriteLine("c {0}", GetCost());
            foreach (int e in edges.ElementEnumerator()) {
                int v, w;
                g.GetEndpoints(e, out v, out w);
                file.WriteLine("e {0} {1} {2}", e, v, w);
            }
            file.Close();
            Console.Error.WriteLine("done.\n");
        }

        public object GetValue() {
            return GetCost();
        }

        public int GetDegree(int v) {
            return degree[v];
        }

        public int LeafCount() {return nleaves;}

        public bool IsBetter(IOptSolution s) {
            return (GetCost() < ((SteinerSolution)s).GetCost());
        }

        public void CopyFrom(IOptSolution _s) {
            SteinerSolution s = (SteinerSolution)_s;

            if (g != s.GetGraph()) {
                SetGraph(s.GetGraph());
            } else {
                Reset();
            }
            foreach (int e in ((SteinerSolution)s).ElementEnumerator()) {
                Insert(e);
            }
            //Console.Error.WriteLine("copy not implemented");   
        }

        public object Clone() {
            SteinerSolution s = new SteinerSolution(g);
            s.CopyFrom(this);
            //Console.Error.WriteLine("Function not implemented");
            //return new SteinerSolution(null);
            return s;
        }

        public SteinerSolution(WeightedGraph _g) {
            SetGraph(_g);
            cost = 0;
            /*
            g = _g;
            int m = g.EdgeCount();
            edges = new UniverseSet(m);
            cost = 0;*/
        }

        public SteinerSolution() {
            g = null;
            cost = 0;
        }

        public void SetGraph(WeightedGraph _g) {
            //g = (WeightedGraph) _g;
            g = _g;
            int m = g.EdgeCount();
            edges = new UniverseSet(m);
            int n = g.VertexCount();

            //Console.Error.WriteLine("graph initialized");
            degree = new int [n+1];
            for (int i=0; i<=n; i++) degree[i] = 0;
            nleaves = 0;
        }

        public int Count() {
            return edges.Count();
        }

        public WeightedGraph GetGraph() {
            return g;
        }

        public IEnumerable<int> ElementEnumerator() {
            return edges.ElementEnumerator();
        }

        public bool Contains(int e) {
            return edges.Contains(e);
        }

        public bool Insert(int e) {
            if (edges.Insert(e)) {
                int a, b;
                g.GetEndpoints(e, out a, out b);                
                int d = (++degree[a]);
                if (d<=2 && !g.IsTerminal(a)) {
                    if (d==1) nleaves ++;
                    else nleaves --;
                }
                d = (++degree[b]);
                if (d <= 2 && !g.IsTerminal(b)) {
                    if (d == 1) nleaves++;
                    else nleaves--;
                }


                //if (((degree[a] ++)==0) && (!g.IsTerminal(a))) nleaves ++;
                //if (((degree[b] ++)==0) && (!g.IsTerminal(b))) nleaves ++;
                //degree[b] ++;
                //if (a == 1) Console.Error.WriteLine("{0}:{1} ", a, degree[a]);
                //if (b == 1) Console.Error.WriteLine("{0}:{1} ", b, degree[b]);
                cost += g.GetCost(e);
                return true;
            } else return false;
        }

        public void Reset() {
            foreach (int e in edges.ElementEnumerator()) {
                int a, b;
                g.GetEndpoints(e, out a, out b);
                //if (((--degree[a]) == 0) && (!g.IsTerminal(a))) nleaves --;
                //if (((--degree[b]) == 0) && (!g.IsTerminal(b))) nleaves ++;
                degree[a] --;
                degree[b] --;
            }
            nleaves = 0;
            edges.Reset();
            cost = 0;
        }

        /// <summary>
        /// Remove edge e from the current solution. Returns true iff
        /// the edge was actually in the solution before.
        /// </summary>
        /// <param name="e">edge to be removed</param>
        /// <returns>true iff the edge did belong to the solution</returns>
        public bool Remove(int e) {
            if (edges.Remove(e)) {
                cost -= g.GetCost(e);
                int a, b;
                g.GetEndpoints(e, out a, out b);

                int d = --degree[a];
                if (d<2 && !g.IsTerminal(a)) {
                    if (d==0) nleaves --;
                    else nleaves ++;
                }

                d = --degree[b];
                if (d<2 && !g.IsTerminal(b)) {
                    if (d==0) nleaves --;
                    else nleaves ++;
                }


                //if ((--degree[a] == 0) && (!g.IsTerminal(a))) nleaves --;
                //if ((--degree[b] == 0) && (!g.IsTerminal(b))) nleaves ++;
                return true;
            } else return false;
        }

        public void Copy(SteinerSolution s) {
            cost = 0;
            edges.Reset();
            foreach (int e in s.ElementEnumerator()) {
                edges.Insert(e);
                cost += g.GetCost(e);
            }
        }

        /// <summary>
        /// Recompute solution cost to account for (possibly new)
        /// edge weights on the graph. Useful after a perturbation 
        /// is added to the graph.
        /// </summary>
        public void UpdateCost() {
            //Console.Error.Write("Initial cost is {0}... ", cost);
            cost = 0;
            foreach (int e in this.ElementEnumerator()) {
                cost += g.GetCost(e);
            }
            ///Console.Error.WriteLine("Updated solution cost is {0}.", cost);
        }

        public ArcCost GetCost() {
            return cost;
        }
    }
}

#endif
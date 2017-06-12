#include "graph.h"

#pragma once


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

class Basics {
public:

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


	static void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
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

	// add all vertices in the current solution to solnodes
	// (vertices with incident edges)
	static void MarkSolutionNodes(Graph &g, SteinerSolution &solution, UniverseSet &solnodes) {
		int n = g.VertexCount();
		for (int v=1; v<=n; v++) {
			if (solution.GetDegree(v)>0) solnodes.Insert(v); 
		}
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


        static void CheckSolution(Graph &g, SteinerSolution &solution) {
            if (!InnerCheck(g, solution, false)) {
                fatal ("Invalid solution");
            }
        }

        static bool InnerCheck(Graph &g, SteinerSolution &solution, bool verbose) {
            int n = g.VertexCount();
            UniverseSet svertices(n);
            UnionFind uf(n);
            Basics::MarkSolutionNodes(g, solution, svertices);
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

};

#include "graph.h"
#include "solution.h"
#include "rfw_random.h"
#include "rfw_stack.h"
#include "LSBasics.h"

#pragma once

class LSVertexInsertion {
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

private:
	static void fatal (const string &msg) {
		Basics::fatal(msg);
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

	static void MarkSolutionVertices (SteinerSolution &solution, UniverseSet &svertices) {
		int n = solution.g->VertexCount();
		for (int v=1; v<=n; v++) {
			if (solution.GetDegree(v)>0) svertices.Insert(v);
		}
	}


public:		
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


			if (verbose) fprintf (stderr, "Inserting initial edges into the tree.\n");
            EdgeCost cost0 = solution.GetCost();
			if (high_verbose) fprintf(stderr, "Initial solution value is %.10f.", cost0);
            STEdgeLinear<EdgeCost> dyntree (n); // dynamic tree data structure
 
			int m = g.EdgeCount();
			for (int e=1; e<=m; e++) { // looping over all edges (could be more efficient)
				if (!solution.Contains(e)) continue;
                int a, b;
                g.GetEndpoints(e, a, b);
                EdgeCost cost = g.GetCost(e);
                if (debug_mode && (!svertices.Contains(a) || !svertices.Contains(b))) fatal ("inconsistent solution representation");
                dyntree.Evert(a); // make a the root of its component
                dyntree.Link(a, b, -cost); // add arc from a to b (b becomes a's parent) 
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

			const bool TIME_EVERTS = false; // measure time taken on "evert" operations?
			RFWTimer evertime;
			if (TIME_EVERTS) {
				evertime.start();
				evertime.pause();
			}

			// perform operations until we have n failures in a row
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
					fprintf (stderr, "Vertex %d improves the solution by %.10f.\n", v, (double)-increase);
                }
            }
            
			EdgeCost improvement = cost0 - curcost;
			
			if (high_verbose && impmoves > 0 && improvement > EDGE_COST_PRECISION) {
                fprintf (stderr, "End of local search; found %d improving moves (tried %d).", impmoves, count);
                fprintf (stderr, "Solution went from %.1f to %.1f: %.10f.\n", cost0, curcost, improvement);
            }


			if (TIME_EVERTS) fprintf (stderr, "%.2f ", 1000*evertime.getTime());

            if (improvement > EDGE_COST_PRECISION) {
				if (!changed) fatal ("inconsistent assessment");
                RebuildSolution(dyntree, solution);
            } else {
				if (changed) fprintf (stderr, "Not using changed version (improvement: %.10f / %.10f).\n", improvement, EDGE_COST_PRECISION);
			}

            if (verbose) fprintf (stderr, "Found %.2f, mst is %.2f.", (double)curcost, (double)solution.GetCost()); 

			if (solneighbors) delete [] solneighbors;
			if (perm) delete [] perm;
            return impmoves;
        }

};

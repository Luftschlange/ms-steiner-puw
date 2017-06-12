#include "graph.h"
#include "solution.h"
#include "rfw_random.h"
#include "LSBasics.h"

class LSVertexElimination {
private:
	static void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
	}

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

public:
	
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

            int root = Basics::PickRandomTerminal(g,random);
			//fprintf (stderr, "r%d ", root);

			//root = 5;
			if (verbose) fprintf (stderr, "Root is %d.\n", root);

            //profit from each possible move
            EdgeCost *profit = new EdgeCost[n + 1];
            for (v=1; v<=n; v++) profit[v] = 0;

            //CheckSolution(solution);
            int visited = Basics::DFS(g, root, solution, dfsdata, stack);
            if (verbose) fprintf (stderr, "DFS visited %d vertices from root %d.", visited, root);

            if (visited < g.TerminalCount()) {
                fprintf (stderr, "Visited %d vertices, but there are %d terminals.", visited, g.TerminalCount());
                fatal ("Something wrong.");
            }

            //find list of all original solution vertices
            Basics::MarkSolutionNodes(g, solution, svertices);


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
                            Basics::MarkComponent(g, stack, id2parc, forbidden);

							if (local_verbose) {fprintf (stderr, "%.2f  ", (double)solution.GetCost()); fflush(stderr);} 

							if (verbose) {
								fprintf (stderr, "Checking after move...\n");
								Basics::CheckSolution(g, solution);
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


};

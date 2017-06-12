#include "binheap.h"
#include "graph.h"
#include "solution.h"
#include "rfw_timer.h"
#include "uf.h"
#include "uset.h"
#include "voronoi.h"
#include "rfw_random.h"
#include "rfw_stack.h"
#include "pairheap.h"
#include "LSBasics.h"

/*
#include "drawer.h"
#include "dual.h"
#include "stedgelinear.h"
#include "buckets.h"
#include <cstring>
#include <cmath>
#include "elite.h"
#include "spgconfig.h"
#include <omp.h>
#include "execution_log.h"
#include "LSVertexInsertion.h"
#include "LSVertexElimination.h"
*/

#pragma once


class LSKeyPath {
private:
	static void fatal (const string &msg) {
		Basics::fatal(msg);
	}



public:

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
            Basics::Prune(g, solution);

            //mark solution, key, and crucial nodes
			//fprintf (stderr, "Warning: check if dereference this makes sense.\n");
            Basics::MarkSolutionNodes(g, solution, *(data.solnodes));
            MarkKeyNodes(g, solution, data.crucialnodes, true);

            if (timer_verbose) fprintf (stderr, "Various marks: %.6f\n", timer.getTime());


            //compute voronoi diagram with all solution vertices as bases
            timer.start();
            Basics::ComputeVoronoi(g, data.voronoi, *(data.solnodes), *(data.binheap), NULL);

			if (timer_verbose) fprintf (stderr, "Voronoi time: %.6f\n", timer.getTime());
            data.vorbkp->Clone(data.voronoi);
            if (verbose) fflush (stderr);
            
            //create heap with boundary edges
            int root = Basics::PickRandomTerminal(g, random);
			//fprintf (stderr, "%d ", root);

            //initialize boundary heaps (WARNING: WOULD LIKE TO ELIMINATE INTRA-PATH EDGES FROM THE START...
            InitBoundaryHeaps(g, data.voronoi, *(data.heap)); //, data.uf);
       
            if (timer_verbose) fprintf (stderr, "heap: %.3f\n", timer.getTime());
            if (verbose) fflush (stderr);


 
            //perform DFS to sort the vertices
            int visited = Basics::DFS(g, root, solution, *(data.dfsdata), *(data.stack));
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
                                        Basics::MarkComponent(g, *(data.stack), data.dfsdata->id2parc, *(data.forbdfs));

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
                                        Basics::MarkComponent(g, *(data.stack), data.dfsdata->id2parc, *(data.forbdfs));
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

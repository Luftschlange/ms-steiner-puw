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

#include "graph.h"
#include "binheap.h"
#include <omp.h>

	struct Triple {
		int u;
		int v;
		int w;
		Triple(int _u, int _v, int _w) {
			u = _u;
			v = _v;
			w = _w;
		}

		inline void Get(int &_u, int &_v, int &_w) const {
			_u = u;
			_v = v;
			_w = w;
		}

	};

struct LabelEntry {
	int hub;
	int dist;

	LabelEntry (int h, int d) {
		hub = h;
		dist = d;
	}

	LabelEntry() {
		hub=-1; dist=-1;
		//fprintf (stderr, "y");
	}
};


class Labels {
public:
	vector<LabelEntry> *labels;
	int n;

	Labels(int _n) {
		n = _n;
		labels = new vector<LabelEntry>[n+1];
		//LabelEntry *le = new LabelEntry[n+1];
		//exit(-1);
		/*
		for (int i=0; i<=n; i++) {
			LabelEntry x;
			(labels[i]).push_back(x);
			labels[i].push_back(x);
			labels[i].push_back(x);
			labels[i].push_back(x);
			fprintf (stderr, "%d ", labels[i].size());
		}*/
	}


	~Labels() {
		fprintf (stderr, "Deleting labels.\n");
		delete [] labels;
	}

	void sort () {
		for (int i=0; i<=n; i++) {
			if (labels[i].size() <= 1) continue;
			std::sort(labels[i].begin(), labels[i].end(), [](const LabelEntry &a, const LabelEntry &b) {return a.hub < b.hub;});
		}
	}

	void OutputLabel (int x) {
		vector<LabelEntry> &label = labels[x];
		fprintf (stderr, "%d: ", x);
		for (int i=0; i<label.size(); i++) {
			fprintf (stderr, "%d:%d ", label[i].hub, label[i].dist);
		}
		fprintf (stderr, "\n");
	}

	void Output() {
		for (int i=1; i<=n; i++) {OutputLabel(i);}
	}

	void AddEntry (int lab, int hub, int dist) {
		labels[lab].push_back(LabelEntry(hub,dist));
	}


};



template <class T> class Matrix {
public:
	T **dist;
	int n;

	inline void fatal (const string &msg) {
		fprintf (stderr, "%s", msg.c_str());
		exit(-1);
	}

	void Init(int _n) {
		if (dist) {
			fatal ("Cannot be initialized twice.\n");
		}
		n = _n;
		dist = new T *[n+1];
		for (int v=0; v<=n; v++) {
			dist[v] = new T [n+1];
		}
	}

	Matrix() {
		n = 0;
		dist = NULL;
	}

	Matrix (T _n) {
		dist = NULL;
		Init(_n);
	}

	inline T Get(int i, int j) const {
		return dist[i][j];
	}

	inline void Set (int i, int j, T value) {
		//if (i<1 || i>n) fatal ("invalid range");
		//if (j<1 || j>n) fatal ("invalid range");
		dist[i][j] = value;
	}

	inline void Add (int i, int j, T value) {
		dist[i][j] += value;
	}

	inline T increment (int i, int j) {
		return (++dist[i][j]);
	}

	inline T decrement (int i, int j) {
		return (--dist[i][j]);
	}

	void Output(FILE *file) {
		for (int v=1; v<=n; v++) {
			for (int w=1; w<=n; w++) {
				int d = dist[v][w];
				if (d >= 1000000000) fprintf (stderr, "<%d,%d>", v,w);
				fprintf (file, "%9d", (int)dist[v][w]);
			}
			fprintf (file, "\n");
		}
	}

	void Reset(T value) {
		for (int v=0; v<=n; v++) {
			for (int w=0; w<=n; w++) {
				dist[v][w] = value;
			}
		}
	}

	~Matrix () {
		fprintf (stderr, "Deleting matrix... ");
		fflush(stderr);
		for (int v=0; v<=n; v++) delete [] dist[v];
		delete [] dist;
		fprintf (stderr, "done.\n");
		fflush(stderr);
	}
};

class BFSData {
private:
	vector<int> inserted;
	vector<int> queue;
	int n; //maximum id
	int nextins; //position of next insertion
	int nextrem; //position of next removal

public:
	void Reset (int _n) {
		n = _n;
		inserted.resize(n+1,0);
		queue.resize(n+1,0);
		nextins = 0;
		nextrem = 0;
	}

	BFSData () {
		n = 0;
	}

	BFSData (int _n) {
		Reset(_n);
	}

	inline int InsertionCount() const {
		return nextins;
	}

	inline int GetElement(int i) const {
		return queue[i];
	}

	inline void Insert(int x) {
		queue[nextins++] = x;
		inserted[x] = true;

	}

	inline int Remove() {
		return (queue[nextrem++]);
	}


	inline bool IsInserted(int x) const {return inserted[x];}

	inline bool IsEmpty() const {
		return (nextins == nextrem);
	}

	// make the list empty
	void Reset() {
		for (int i=0; i<nextins; i++) {inserted[queue[i]] = false;}
		nextins = nextrem = 0;
	}

	~BFSData() {}
};




class CoverCounters {
public:
	int n;
	Matrix<int> sourcecount; //sourcecount(h,v): number of uncovered paths starting at v covered by h
	Matrix<int> targetcount; //targetcount(h,v): number of uncovered paths ending at v covered by h
	vector<int> fwdcount;    //number of forward labels h would be added to
	vector<int> backcount;   //number of backward labels h would be added to

	void Reset() {
		sourcecount.Reset(0);
		targetcount.Reset(0);
		for (int i=0; i<=n; i++) {
			fwdcount[i] = backcount[i] = 0;
		}
	}

	CoverCounters(int _n) {
		n = _n;
		sourcecount.Init(n);
		targetcount.Init(n);
		fwdcount.resize(n+1,0);
		backcount.resize(n+1,0);
		Reset();
	}

	inline void AddPath (int u, int v, int w) {
		if (sourcecount.Get(v,u)==0) fwdcount[v] ++; 
		sourcecount.Add(v,u,1);
		if (targetcount.Get(v,w)==0) backcount[v]++; 
		targetcount.Add(v,w,1);
	}

	inline void AddPath2 (int u, int v, int w) {
		if (sourcecount.increment(v,u)==1) fwdcount[v]++;
		if (targetcount.increment(v,w)==1) backcount[v]++;
		//sourcecount.increment(v,u); //==1) fwdcount[v]++;
		//targetcount.increment(v,w); //==1) backcount[v]++;
	}

	inline void AddPathSimple (Triple &t) {
		sourcecount.increment(t.v,t.u);
		targetcount.increment(t.v,t.w);
	}

	inline void AddPathSimple (int u, int v, int w) {
		sourcecount.increment(v,u);
		targetcount.increment(v,w);
	}

	inline void RemovePathSimple (int u, int v, int w) {
		sourcecount.decrement(v,u);
		targetcount.decrement(v,w);
	}



	inline void AddPath2 (Triple &t) {
		if (sourcecount.increment(t.v,t.u)==1) fwdcount[t.v]++;
		if (targetcount.increment(t.v,t.w)==1) backcount[t.v]++;
		//sourcecount.increment(v,u); //==1) fwdcount[v]++;
		//targetcount.increment(v,w); //==1) backcount[v]++;
	}

	void FixLabelCounters(int v) {
		int fc = 0;
		int bc = 0;
		for (int x=1; x<=n; x++) {
			if (sourcecount.Get(v,x)>0) fc ++;
			if (targetcount.Get(v,x)>0) bc ++;
		}
		fwdcount[v] = fc;
		backcount[v]= bc;
		//fprintf (stderr, "%d %d ", fc - fwdcount[v], bc - backcount[v]);
	}

	void RemovePath (int u, int v, int w) {
		sourcecount.Add(v,u,-1);
		if (sourcecount.Get(v,u)==0) fwdcount[v] --; 
		targetcount.Add(v,w,-1);
		if (targetcount.Get(v,w)==0) backcount[v]--; 
	}

	void Output(FILE *file) {
		sourcecount.Output(file);
		targetcount.Output(file);

	}


};

class LabelSolver {
public:
	int ncount;

	void DensestSubgraph(Graph &g) {
		fprintf (stderr, "Should be computing densest subgraph.\n");

		int n = g.VertexCount();

		vector<int> curdeg(n+1,0);
		vector<bool> present(n+1,true);

		BinaryHeap<double> heap(n);
		for (int v=1; v<=n; v++) {
			curdeg[v] = g.GetDegree(v);
			heap.Insert(v, curdeg[v]);
		}

		int curedges = g.EdgeCount();
		int curvertices = n;

		double maxdensity = 0;

		// QUADRATIC IMPLEMENTATION!
		int count = 0;
		bool verbose = false;
		while (!heap.IsEmpty()) {
			double curdensity = (double) curedges / (double) curvertices;
			if (curdensity > maxdensity) maxdensity = curdensity;
			if (verbose) fprintf (stderr, "%6d: Started on graph with %d vertices and %d edges, density %.2f.\n", count, curvertices, curedges, curdensity);
			count ++;
			unsigned int v;
			double deg;
			heap.RemoveFirst(v, deg);
			if (!present[v]) fprintf (stderr, "Invalid removal.\n");


			// find lowest degree vertex to remove
			/*
			int v = -1;
			for (int t=1; t<=n; t++) {
				if (!present[t]) continue;
				if (v<=0 || curdeg[t] < curdeg[v]) {v = t;}
			}*/

			if (verbose) fprintf (stderr, "Removing vertex %d with degree %d: ", v, curdeg[v]); 
			present[v] = false;
			curvertices--;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //neighbor
				if (!present[w]) continue;
				curdeg[v] --;
				curdeg[w] --;
				heap.Insert(w, curdeg[w]);
				curedges --;
            }
		}
		fprintf (stderr, "Maximum density = %.2f\n", maxdensity);
	}

	void ComputeCoverage2(Graph &g, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, bool LABEL_GREEDY, CoverCounters *covercount) {
		int n = g.VertexCount();
		RFWTimer timer(true);


		// compute curcover[v]: number of pairs that are hit by v
		for (int u=1; u<=n; u++) {
			//int duv = dt.Get(u,v);
			for (int w=1; w<=n; w++) {
				if (covered.Get(u,w)) continue; //unreachable stuff
				int duw = dt.Get(u,w);
				for (int v=1; v<=n; v++) {
					if (duw == dt.Get(u,v) + dt.Get(v,w)) {
						//path u--w covered by v
						curcover[v] ++;
						if (LABEL_GREEDY) {covercount->AddPath(u,v,w);}
					}
				}
			}
			//curcover[v] = vcount;
			if (u % 10 == 0) {
				fprintf (stderr, "%d(%.2f) ", u, n * timer.getTime() / u);
				fflush(stderr);
			}
		}
	}


	void MatrixComputeCoverage(Graph &g, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, bool LABEL_GREEDY, CoverCounters *covercount) {
		fprintf (stderr, "Matrix-based coverage computation started.\n");
		int n = g.VertexCount();
		RFWTimer timer(true);

		int threads = omp_get_max_threads();
		fprintf (stderr, "Should be working with %d threads.\n", threads);
		vector<Triple> *triples = new vector<Triple> [threads];


		// making v the external loop seems to be advantageous; w in the internal loop because 
		// the access to the distance table is sequential (we only look at rows u and v)

		// compute curcover[v]: number of pairs that are hit by v
		// Takes Theta(n^3) time because we are not using the graph---but reasonably efficient
#pragma omp parallel for schedule(dynamic)
		for (int v=1; v<=n; v++) {
			int tid = omp_get_thread_num();
			vector<Triple> &curtrip = triples[tid];
			curtrip.clear();
			int vcount = 0;
			for (int u=1; u<=n; u++) {
				int duv = dt.Get(u,v);
				for (int w=1; w<=n; w++) {
					// is path u-v covered by w?
					if (dt.Get(u,w) == duv + dt.Get(v,w)) { //awesome sequential access on u and w arrays
						if (covered.Get(u,w)) continue; //unreachable stuff
						vcount ++;
						if (LABEL_GREEDY) {
							curtrip.push_back(Triple(u,v,w));
//#pragma omp critical 
//{						
//							covercount->AddPath(u,v,w);
//}
						}
					}
				}
			}
			curcover[v] = vcount;
			int tcount = curtrip.size();
			if (tcount>0) {
#pragma omp critical
				for (int i=0; i<tcount; i++) {
					covercount->AddPath2(curtrip[i]);
				}
			}

			if (v % 10 == 0) {
				fprintf (stderr, "%d(%.2f) ", v, n * timer.getTime() / v);
				fflush(stderr);
			}
		}
	}

	//-----------------------------------------------------
	// Do search from v but only following vertices w such that:
	// - u-v-w is a shortest path
	// - u-v-w is previously uncovered
	// - w is a target we care about (BAD!) 
	// WARNING: THIS IS COMPLETELY WRONG.
	//-----------------------------------------------------------
	void TargetDAGRestricted (Graph &g, int u, int v, int bestv, Matrix<int> &dt, Matrix<int> &covered, vector<int> &targets, BFSData &bfsdata) {
		bfsdata.Reset();
		bfsdata.Insert(v);
		int n = g.VertexCount();
		int duv = dt.Get(u,v);

		while (!bfsdata.IsEmpty()) {
			int x = bfsdata.Remove();
			SPGArc *a, *end;
			int dx = dt.Get(u,x); //duv + dt.Get(v,x); //distance to x through v
			//if (dx != dt.Get(u,x)) fprintf (stderr, "!");
			for (g.GetBounds(x,a,end); a<end; a++) {
				int w = a->head;
				if (bfsdata.IsInserted(w)) continue; //already seen w in the search
				//if (!targets[w]) continue; //w is not a target
				//if (covered.Get(u,w)>0 && covered.Get(u,w)!=bestv) continue;
				//fprintf (stderr, ".");
				
				int c = covered.Get(u,w);
				if (c==bestv || c==0) {
					if (dx + a->cost == dt.Get(u,w)) bfsdata.Insert(w);
				}
				//} else badarcs ++; // arc does not lead to improved distance...
			}
		}
	}



	// run Dijkstra from v, visiting only vertices w such that u-v-w is a shortest path
	void TargetDAG (Graph &g, int u, int v, Matrix<int> &dt, Matrix<int> &covered, BFSData &bfsdata) {
		bfsdata.Reset();
		bfsdata.Insert(v);
		int n = g.VertexCount();
		int duv = dt.Get(u,v);

		while (!bfsdata.IsEmpty()) {
			int x = bfsdata.Remove();
			SPGArc *a, *end;
			int dx = dt.Get(u,x); //duv + dt.Get(v,x); //distance to x through v
			//if (dx != dt.Get(u,x)) fprintf (stderr, "!");
			for (g.GetBounds(x,a,end); a<end; a++) {
				int w = a->head;
				if (bfsdata.IsInserted(w)) continue;
				if (covered.Get(u,w)) continue;
				if (dx + a->cost == dt.Get(u,w)) bfsdata.Insert(w);
				//} else badarcs ++; // arc does not lead to improved distance...
			}
		}
	}

	void TraverseDAG (Graph &g, int s, Matrix<int> &dt, BFSData &bfsdata) {
		bfsdata.Reset();
		bfsdata.Insert(s); //start a search from s
		int n = g.VertexCount();
		int debugcount = 0;
		int badarcs = 0;

		while (!bfsdata.IsEmpty()) {
			int v = bfsdata.Remove();
			debugcount ++;
			if (debugcount > n) {
				fprintf (stderr, "Infinite loop.\n");
				exit(-1);
			}

			SPGArc *a, *end;
			int dsv = dt.Get(s,v);
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (dsv + a->cost == dt.Get(s,w)) {
					if (bfsdata.IsInserted(w)) continue;
					bfsdata.Insert(w);
				} else badarcs ++; // arc does not lead to improved distance...
			}
		}

		//fprintf (stderr, "<%.5f> ", 100 * (double)badarcs / (2.0 * g.EdgeCount())); 
	}


	void TraverseDAG(Graph &g, int s, int t, Matrix<int> &dt, BFSData &bfsdata) {
		bfsdata.Reset();
		bfsdata.Insert(s); //start a search from s
		int n = g.VertexCount();
		int dst = dt.Get(s,t);
		int debugcount = 0;

		while (!bfsdata.IsEmpty()) {
			int v = bfsdata.Remove();
			debugcount ++;
			if (debugcount > n) {
				fprintf (stderr, "Infinite loop.\n");
				exit(-1);
			}
			// v hits the u-w path
			//		curcover[v] ++;
			//int dsv = dt.Get(s,v);
			//if (dsv + dt.Get(v,t) != dst) continue;

			int dsv = dt.Get(s,v);

			int gap = dst - dsv; // this is how much is missing from the path

			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int x = a->head;
				if (bfsdata.IsInserted(x)) continue; //already seen this vertex
				//if (dt.Get(s,x) + dt.Get(x,t) == dst) bfsdata.Insert(x);

				if (a->cost + dt.Get(x,t) == gap) bfsdata.Insert(x); //WARNING: dist x--t is expensive
			}
		}
	}

	// increase curcover[v] for all vertices v that are on the shortest s-t path.
	void UpdateCoverage (vector<int> &curcover, int s, int t, Matrix<int> &dt, vector<int> &vlist, bool LABEL_GREEDY, CoverCounters *covercount) { 
		int count = vlist.size();
		int curcount = 0;

		//fprintf (stderr, "<%d> ", count);
		int dst = dt.Get(s,t);
		for (int i=0; i<count; i++) {
			int v = vlist[i];
			if (dst == dt.Get(s,v) + dt.Get(v,t)) {
				curcount ++;
				curcover[v] ++;
				if (LABEL_GREEDY) covercount->AddPath(s,v,t);
			}
		}
		//fprintf (stderr, "%.0f ", 100.0 * curcount / count);
	}


	void FancyComputeCoverage(Graph &g, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, bool LABEL_GREEDY, CoverCounters *covercount) {
		fprintf (stderr, "Computing coverage using fancy algorithm!\n");
		int n = g.VertexCount();
		BFSData fulldata(n);
		BFSData localdata(n);
		vector<int> tdone(n+1, 0);
		
		vector<int> relevant;
		vector<int> todo;
		relevant.reserve(n+1);
		todo.reserve(n+1);
		//vector<int> bestcand(n+1,0);
		
		RFWTimer timer (true);

		for (int s=1; s<=n; s++) {
			for (int v=1; v<=n; v++) tdone[v] = 0;
			TraverseDAG(g,s,dt,fulldata);
			int fcount = fulldata.InsertionCount();

			int searches = 0;
			// visit vertices in decreasing order of "distance" from s
			for (int i=fcount-1; i>=0; i--) {
				int t = fulldata.GetElement(i);
				if (tdone[t]) continue;
				tdone[t] = 1;

				searches ++;
				TraverseDAG(g,s,t,dt,localdata);
				int icount = localdata.InsertionCount();
				//fprintf (stderr, "[%d] ", icount);

				relevant.clear();
				todo.clear();

				// look at all elements visited by the search
				for (int i=0; i<icount; i++) {
					int v = localdata.GetElement(i);
					curcover[v]++;
					if (LABEL_GREEDY) covercount->AddPath(s,v,t);
					relevant.push_back(v);
					if (!tdone[v]) todo.push_back(v); //got nontrivial information about this vertex
				}

				// Let v be a vertex that hits the s-t pair.
				// Every vertex that hits s-v must also hit s-t.
				// To find the vertices that hit s-v, we only need
				// to look among those that hit s-t. This is done below.
				for (int i=todo.size()-1; i>=0; i--) {
					int v = todo[i];
					UpdateCoverage(curcover, s, v, dt, relevant, LABEL_GREEDY, covercount);
					tdone[v] = 1;
				}
			}

			//fprintf (stderr, "<%.0f:%d>\n", 100.0 * (double)searches / (double)fcount, fcount);
			if (s%10 == 0) {
				fprintf (stderr, "%d (%.2f) ", s, n * timer.getTime() / (double)s);
				fflush(stderr);
			}
		}

	}


	void GraphComputeCoverageCenter (Graph &g, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, bool LABEL_GREEDY, CoverCounters *covercount) {
		int n = g.VertexCount();
		RFWTimer timer (true);
		RFWTimer crittimer (false);

		int threads = omp_get_max_threads();
		fprintf (stderr, "Should be working with %d threads.\n", threads);
		
		// we're looking for u-v-w paths
#pragma omp parallel for schedule(dynamic)
		for (int v=1; v<=n; v++) {
			int tid = omp_get_thread_num();
			BFSData bfsdata(n);
			int vcount = 0;

			for (int u=1; u<=n; u++) {
				if (covered.Get(u,v)) continue;
				TargetDAG(g, u, v, dt, covered, bfsdata); //could do a little clean-up here...
				int wcount = bfsdata.InsertionCount();
				for (int i=0; i<wcount; i++) {
					int w = bfsdata.GetElement(i);
					vcount ++;
					if (LABEL_GREEDY) covercount->AddPathSimple(u,v,w); //this only affects v
				}
			}
			curcover[v] = vcount;
			if (LABEL_GREEDY) covercount->FixLabelCounters(v);

			if (v%10 == 0) {
				fprintf (stderr, " %d(%.2f)", v, (double)n*timer.getTime() / (double)v);
				fflush(stderr);
			}
		}
	}

	void GraphComputeCoverage(Graph &g, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, bool LABEL_GREEDY, CoverCounters *covercount) {
		int n = g.VertexCount();
		RFWTimer timer (true);
		RFWTimer crittimer (false);

		const bool FANCY = true;

		int threads = omp_get_max_threads();
		fprintf (stderr, "Should be working with %d threads.\n", threads);
		vector<Triple> *triples = new vector<Triple> [threads];

		fprintf (stderr, "Computing coverage using graph algorithm!\n");
#pragma omp parallel for schedule(dynamic)
		for (int u=1; u<=n; u++) {
			int tid = omp_get_thread_num();
			if (FANCY) {
				//triples[tid].reserve(n);
				triples[tid].clear();
			}

			BFSData bfsdata(n);
			for (int w=1; w<=n; w++) {
				if (covered.Get(u,w)) continue; //these are unreachable!
				
				//we'll visit all vertices that hit the u,w path
				TraverseDAG(g,u,w,dt,bfsdata);

				int icount = bfsdata.InsertionCount();
				int v;

				if (FANCY) {
					for (int i=0; i<icount; i++) {
						v = bfsdata.GetElement(i);
						triples[tid].push_back(Triple(u,v,w));
					} 
				} else {
//warning: this is too fine-grained
#pragma omp critical
					for (int i=0; i<icount; i++) {
						curcover[v=bfsdata.GetElement(i)]++;
						if (LABEL_GREEDY) {covercount->AddPath(u,v,w);} //this only affects stuff that have to do with v
					}
				}
			}

			if (FANCY) {
				vector<Triple> &curt = triples[tid];
				//std::sort(curt.begin(), curt.end(), [](const Triple &a, const Triple &b) {return a.v < b.v;});
#pragma omp critical
{
				crittimer.resume();
				int tcount = curt.size();
				for (int i=0; i<tcount; i++) {
					int u, v, w;
					curt[i].Get(u,v,w);

					// this will only affect memory positions indexed by v, so it could be made parallel
					//fprintf (stderr, "[%d,%d,%d] ", u, v, w);
					curcover[v]++;
					if (LABEL_GREEDY) {covercount->AddPath2(u,v,w);}
				}
				crittimer.pause();
}
			}



			if (u%10 == 0) {
				fprintf (stderr, " %d(%.2f)", u, (double)n*timer.getTime() / (double)u);
				fprintf (stderr, "[%.0f] ", 100.0 * crittimer.getTime() / timer.getTime());
				fflush(stderr);
			}
		}

		delete [] triples;
	}

	struct ThreadData {
		vector<Triple> triples;
		BFSData coveragebfs;
		vector<int> targetlist;
		int vcount;
		double totalcovered;

		void Reset(int n) {
			triples.reserve(n+1);
			targetlist.reserve(n+1);
			coveragebfs.Reset(n);
			vcount = 0;
			totalcovered = 0;
			
		}
	};

	struct DecreaseThreadData {
		//vector<Triple> triples;
		BFSData coveragebfs;
		//vector<int> targetlist;
		//int vcount;
		vector<int> target;
		vector<int> targetlist;
		double totalcovered;

		void Reset(int n) {
			//triples.reserve(n+1);
			targetlist.reserve(n+1);
			target.resize(n+1,0);
			coveragebfs.Reset(n);
			//vcount = 0;
			totalcovered = 0;
			
		}
	};

	void DecreaseCoverageMatrixCenter (Graph &g, int bestv, vector<pair<int,int>> &pairs, vector<int> &curcover, Matrix<int> &dt, CoverCounters *covercount, Matrix<int> &covered, bool LABEL_GREEDY) {
		fprintf (stderr, "DCMC ");
		int npairs = pairs.size();
		if (npairs == 0) fprintf (stderr, "!");
		
		int n = g.VertexCount();
		vector<int> ulist;
		vector<int> first (n+1,-1);
		pairs.push_back(make_pair(-1,-1)); //sentinel

		for (int p=0; p<npairs; p++) {
			int u = pairs[p].first;
			int w = pairs[p].second;
			if (p==0 || (pairs[p-1].first != u)) {
				first[u] = p;
				ulist.push_back(u);
			}
		}

		sort(ulist.begin(), ulist.end());
		fprintf (stderr, "pair%d src%d ", pairs.size(), ulist.size());

		// we're looking for u-v-w paths
//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel for schedule(dynamic)
		for (int v=1; v<=n; v++) {
			int ucount = ulist.size();
			int vcount = 0;
			for (int iu=0; iu<ucount; iu++) {
				int u = ulist[iu];
				if ((covered.Get(u,v)>0) && (covered.Get(u,v)!=bestv)) continue;

				int duv = dt.Get(u,v);

				for (int j=first[u]; pairs[j].first == u; j++) {
					int w = pairs[j].second;
					//targetlist.push_back(w);
					//targets[w] = 1;

					if (duv + dt.Get(v,w) == dt.Get(u,w)) {
						vcount ++;
						if (LABEL_GREEDY) covercount->RemovePathSimple(u,v,w); //this only affects v
					}
				}
			}
			curcover[v] -= vcount;
			if (LABEL_GREEDY) covercount->FixLabelCounters(v);
		}
	}


	//-------------------------------------------------------------------------
	// Decrease coverage of all vertices that hit a give set of shortest paths
	//-------------------------------------------------------------------------
	void DecreaseCoverageGraphCenter (Graph &g, int bestv, vector<pair<int,int>> &pairs, vector<int> &curcover, Matrix<int> &dt, CoverCounters *covercount, Matrix<int> &covered, bool LABEL_GREEDY) {

		RFWTimer timer(true);
		fflush(stderr);
		int n = g.VertexCount();
		int nthreads = omp_get_max_threads();
		fprintf (stderr, "Should be decreasing coverage with %d threads.\n", nthreads);

		DecreaseThreadData *tdata = new DecreaseThreadData[nthreads];

		for (int t=0; t<nthreads; t++) {
			tdata[t].Reset(n);
		}

		vector<int> ulist;
		vector<int> first (n+1,-1);
		int npairs = pairs.size();
		pairs.push_back(make_pair(-1,-1)); //sentinel

		for (int p=0; p<npairs; p++) {
			int u = pairs[p].first;
			int w = pairs[p].second;
			if (p==0 || (pairs[p-1].first != u)) {
				first[u] = p;
				ulist.push_back(u);
			}
		}

		sort(ulist.begin(), ulist.end());
		fprintf (stderr, "There are %d pairs from %d sources.\n", pairs.size(), ulist.size());

		// we're looking for u-v-w paths
#pragma omp parallel for schedule(dynamic)
		for (int v=1; v<=n; v++) {
			int tid = omp_get_thread_num();
			BFSData &bfsdata = tdata[tid].coveragebfs;
			vector<int> &targets = tdata[tid].target;
			vector<int> &targetlist = tdata[tid].targetlist;
			int vcount = 0;

			int ucount = ulist.size();
			for (int iu=0; iu<ucount; iu++) {
				int u = ulist[iu];
				//if (covered.Get(u,v)!=bestv) continue; 
				
				//int c = covered.Get(u,v);
				//if (c>0 && c!=bestv) continue; // path u-v has already been covered---every descendent will be covered as well

				// build incidence vecotr of targets
				targetlist.clear();
				for (int j=first[u]; pairs[j].first == u; j++) {
					int w = pairs[j].second;
					targetlist.push_back(w);
					targets[w] = 1;  
				}

				//int wcount = targetlist.size();
				TargetDAGRestricted(g, u, v, bestv, dt, covered, targets, bfsdata);
				
				int wcount = bfsdata.InsertionCount();
				//fprintf (stderr, "%d ", wcount);
				for (int i=0; i<wcount; i++) {
					int w = bfsdata.GetElement(i);
					if (targets[w]) {
						vcount ++;
						if (LABEL_GREEDY) covercount->RemovePathSimple(u,v,w); //this only affects v
					}
				}

				for (int j=first[u]; pairs[j].first == u; j++) {
					int w = pairs[j].second;
					targets[w] = 0;
				}

				for (int w=1; w<=n; w++) {
					if (targets[w]) fprintf (stderr, "B%d ", w);
				}
			}
			if (v==105) fprintf (stderr, "Updating %d by %d to %d.\n", v, vcount, curcover[v] - vcount);
			curcover[v] -= vcount;



			// WARNING: BAD!
			if (LABEL_GREEDY) covercount->FixLabelCounters(v);

			if (v%10 == 0) {
				fprintf (stderr, " %d(%.2f)", v, (double)n*timer.getTime() / (double)v);
				fflush(stderr);
			}
		}

		delete [] tdata;

	}


	//---------------------------------------------------------------------------------------------------
	// Decrease coverage of all vertices that current hit an s-t path, with t \in targetlist.
	// Does so by performing a BFS from s, but only scanning vertices that cover at least one such path.
	//---------------------------------------------------------------------------------------------------

	void DecreaseCoverageGraph (Graph &g, int s, int bestv, vector<int> &curcover, vector<int> &targetlist, Matrix<int> &dt, CoverCounters *covercount, BFSData &bfs, Matrix<int> &covered, vector<Triple> &triples) {
		int tcount = targetlist.size();
		bfs.Reset();
		bfs.Insert(s); //start starts at s


		
		//int n = g.VertexCount();
		//vector<int> dst(n+1);
		//for (int t=1; t<=n; t++) dst[t] = dt.Get(s,t);

		while (!bfs.IsEmpty()) {
			int v = bfs.Remove();
			int c = covered.Get(s,v);
			if (c>0 && c!=bestv) continue; //if the s-v path is previously covered, any superpath was already covered as well

			int dsv = dt.Get(s,v);
			bool scan = false; //does it make sense to scan this vertex?
			if (v == bestv) {scan = true;} // vertex just picked; it's already been taken care of, but we need to keep scanning
			else {
				// check for which t vertex v is on the shortest s-t path
				const bool DEBUG_COUNT = false;
				int count = 0;
				for (int i=0; i<tcount; i++) {
					int t = targetlist[i];
					if (dt.Get(s,t) == dsv + dt.Get(v,t)) { //could precompute dist(s,t)
					//if (dst[t] == dsv + dt.Get(v,t)) { //could precompute dist(s,t)
						scan = true; 
						if (DEBUG_COUNT) count ++;
						triples.push_back(Triple(s,v,t));
					}
				}
				if (DEBUG_COUNT) fprintf (stderr, "%.0f ", 100.0 * (double)count / (double)tcount);
			}

			//fprintf (stderr, "%d", scan);

			if (scan) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
					if (bfs.IsInserted(w)) continue;
					//if (dsv + a->cost != dt.Get(s,w)) continue; //not really relevant if it is a shortest path
					bfs.Insert(w);
				}
			}
		}
		//fprintf (stderr, "\n", stderr);

	}

	void DecreaseCoverageMatrix(int n, int s, vector<int> &curcover, int bestv, vector<int> &targetlist, Matrix<int> &dt, CoverCounters *covercount) {
		// for fixed s and considering all w such that s-w is newly convered by bestv,
		// we have a path u-w.
		// For each such path, we must find all other elements x that hit u-w as well
		int tcount = targetlist.size();
		for (int x=1; x<=n; x++) { // try all possible middle vertices
			if (curcover[x] == 0) continue; //vertex no longer relevant
			if (x==bestv) continue;
			int dux = dt.Get(s,x); // prefix

			for (int i=0; i<tcount; i++) {
				int w = targetlist[i];
				if (dt.Get(s,w) == dux + dt.Get(x,w)) {
					curcover[x] --;
					if (covercount) {covercount->RemovePath(s,x,w);}
				}
			}
		}
	}

	int SetActualThreads(int maxthreads, bool verbose) {
		int available_threads = omp_get_max_threads();
		int threads;
		if (maxthreads == -1) threads = available_threads;
		else threads = min(maxthreads, available_threads);

		/*if (maxthreads != -1) {cores = min(cores, numCores);}*/
		omp_set_num_threads(threads);

		if (verbose) {
			fprintf (stderr, "threadsmax %d\n", maxthreads);
			fprintf (stderr, "threadsavailable %d\n", available_threads);
			fprintf (stderr, "threadsactual %d\n", threads);
			fflush (stderr);
		}

		return threads;
	}

	void BFS (Graph &g, int r, BFSData &bfsdata) {
		bfsdata.Reset();
		bfsdata.Insert(r);

		while (!bfsdata.IsEmpty()) {
			int v = bfsdata.Remove();

			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (!bfsdata.IsInserted(w)) bfsdata.Insert(w);
			}
		}
	}

	void OutputMETIS(FILE *file, Graph &g, vector<int> &old2new) {
		int newn = 0;
		int oldn = g.VertexCount();
		fprintf (stderr, "Outputting metis file.\n");
		fflush(stderr);

		for (int v=1; v<=oldn; v++) {
			newn = std::max(old2new[v], newn);
		}
		fprintf (stderr, "Outputting graph with %d vertices.\n", newn);
		fflush(stderr);

		int acount = 0;
		for (int step=0; step<=1; step++) {
			if (step==1) {
				fprintf (file, "%d %d 0\n", newn, acount/2);
			}
			for (int v=1; v<=oldn; v++) {
				if (old2new[v] <= 0) continue;
				fprintf (stderr, "<%d> ", old2new[v]);
				int count = 0;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = old2new[a->head];
					if (w <= 0) continue;
					if (step==0) acount ++;
					else {
						if (count++ > 0) fprintf (file, " ");
						fprintf (file, "%d", w);
					}
				}
				if (step==1) fprintf (file, "\n");
			}
			if (step==0) {
				fprintf (stderr, "It appears the graph has %d arcs (%d edges) and %d vertices.\n", acount, acount/2, newn);
			}
		}
	}

	void ExtractGiantComponent(Graph &g) {
		int n = g.VertexCount();
		BFSData bfsdata(n);
		vector<int> component(n+1,-1);
		vector<int> compsize(n+1,0);

		int gc = 0;

		fprintf (stderr, "Extracting giant component.\n");
		fflush(stderr);

		for (int v=1; v<=n; v++) {
			if (component[v]>=0) continue;
			BFS(g,v,bfsdata);
			int size = bfsdata.InsertionCount();
			for (int i=0; i<size; i++) {
				int w = bfsdata.GetElement(i);
				component[w] = v;
			}
			compsize[v] = size;
			if (gc==0 || size>compsize[gc]) gc = v;
			if (size > 10) fprintf (stderr, "Component %d has size %d.\n", v, size);
		}

		fprintf (stderr, "Biggest component is %d, with size %d.\n", gc, compsize[gc]);
		fflush(stderr);

		vector<int> old2new(n+1,-1);
		int nextid = 1;
		for (int v=1; v<=n; v++) {
			if (component[v]==gc) {
				old2new[v] = nextid++;
			}
		}
		fprintf (stderr, "About to output.\n");
		fflush(stderr);
		FILE *file = fopen ("giantcc.graph", "w");
		if (!file) {
			fprintf (stderr, "could not open file\n");
			exit(-1);
		}
		OutputMETIS(file, g, old2new);
		fclose(file);
	}


	void UpdateCoverageInfo(Graph &g, int bestv, vector<int> &curcover, Matrix<int> &covered, Matrix<int> &dt, vector<int> &forward, vector<int> &backward, 
		bool LABEL_GREEDY, CoverCounters *covercount, Labels &flabels, Labels &blabels, int &addcount, 
		double &fulltotalcovered, int &fullvcount, int origcount, bool DECREASE_COUNTERS) {


		int n = g.VertexCount();
		int threads = omp_get_max_threads();
		//fprintf (stderr, "Should be updating using %d threads.\n", threads);
		/*
		vector<Triple> *triples = new vector<Triple> [threads];
		BFSData coveragebfs(n);
		vector<int> targetlist;
		targetlist.reserve(n+1);*/

		vector<pair<int,int>> pairs; //these are the pairs of vertices covered by adding v
		ThreadData *tdata = new ThreadData [threads];
		for (int t=0; t<threads; t++) {
			tdata[t].Reset(n);
		}

		bool CENTER_BASED_UPDATE = false;

		RFWTimer localtimer(true);
		int timerstep = 500000;
		int timertarget = timerstep;
		
#pragma omp parallel for schedule(dynamic)
		for (int u=1; u<=n; u++) {
			if (covered.Get(u,bestv)) continue; //speedup! already covered this part 
			int duv = dt.Get(u,bestv);

			int tid = omp_get_thread_num();
			vector<Triple> &triples = tdata[tid].triples;
			vector<int> &targetlist = tdata[tid].targetlist;
			double &totalcovered = tdata[tid].totalcovered;
			//int &vcount = tdata[tid].vcount;
			int vcount = 0;
			triples.clear();
			targetlist.clear();

			// do we really need to loop over everything to mark?
			for (int w=1; w<=n; w++) {
				if (covered.Get(u,w)) continue; // fprintf (stderr, "!");
				if (dt.Get(u,w) == duv + dt.Get(bestv,w)) {
					//if (LABEL_GREEDY) {covercount->RemovePath(u,bestv,w);}
					triples.push_back(Triple(u,bestv,w)); //found a path we want to cover
					targetlist.push_back(w);
					// decrease conters of everybody that covers this path

					// could save something here by having two (pre-)sorted lists at x
					// the first has distances to x, the second distances from x
					// we can do binary search for the first match, then just walk within the desired range;
					// could delete entries as paths are covered
					totalcovered --;
					vcount ++;
					covered.Set(u,w,bestv); //only one thread will try to updated u ---- no races
				}
			}

			// now remember all paths starting at u
			bool SIMPLE = false;
			if (SIMPLE) {
				fprintf (stderr, "NOT SUPPORTED ANYMORE.\n");
				exit(-1);
				int tcount = targetlist.size();
				for (int i=0; i<tcount; i++) {
					int w = targetlist[i];
					int duw = dt.Get(u,w); //must find all x that hit the u-w path
					for (int x=1; x<=n; x++) {
						if (duw == dt.Get(u,x) + dt.Get(x,w)) {
							if (curcover[x] == 0) continue; //x doesn't cover anything anyway
							if (x==bestv) continue;
							curcover[x] --;
							if (LABEL_GREEDY) {covercount->RemovePath(u,x,w);}
						}
					}
				}
			} else {
				//DecreaseCoverageMatrix(n, u, curcover, bestv, targetlist, dt, &covercount);
				if (DECREASE_COUNTERS && !CENTER_BASED_UPDATE) {
					DecreaseCoverageGraph(g, u, bestv, curcover, targetlist, dt, covercount, tdata[tid].coveragebfs, covered, triples);
					//DecreaseCoverageMatrix(g.VertexCount(), u, curcover, bestv, targetlist, dt, covercount); //, tdata[tid].coveragebfs, covered, triples);
				}

				int tcount = triples.size();
				if (CENTER_BASED_UPDATE) {
#pragma omp critical 
					{
						if (!forward[u]) {
							forward[u] = true; 
							flabels.AddEntry(u,bestv,dt.Get(u,bestv));
							addcount ++;
						}
						fullvcount += vcount;
						int wcount = targetlist.size();
						for (int i=0; i<wcount; i++) {
							int w = targetlist[i];
							pairs.push_back(make_pair(u,w));
							if (!backward[w]) {
								backward[w] = true; 
								blabels.AddEntry(w,bestv,dt.Get(bestv,w));
								addcount ++;
							}
						}
						//fprintf (stderr, "Pairs has %d vertices.\n", pairs.size());
					}
					continue;
				}
#pragma omp critical
{
				fullvcount += vcount;
				for (int i=0; i<tcount; i++) {
					int u = triples[i].u;
					int v = triples[i].v;
					int w = triples[i].w;
					curcover[v] --;
					if (LABEL_GREEDY) {covercount->RemovePath(u,v,w);}
					if (v == bestv) {
						if (!forward[u]) {
							forward[u] = true; 
							flabels.AddEntry(u,bestv,dt.Get(u,bestv));
							addcount ++;
						}
						if (!backward[w]) {
							backward[w] = true; 
							blabels.AddEntry(w,bestv,dt.Get(bestv,w));
							addcount ++;
						}
					}
				}
				if (fullvcount > timertarget) {
					while (fullvcount > timertarget) {timertarget += timerstep;}
					double t = localtimer.getTime();
					fprintf (stderr, "%.0f/%.0f ", t, (double)origcount * t / (double)fullvcount); 
					fflush(stderr);
				}
}				

			}
		}

		for (int t=0; t<threads; t++) {
			fulltotalcovered += tdata[t].totalcovered;
			//fullvcount += tdata[t].vcount;
		}

		if (DECREASE_COUNTERS && CENTER_BASED_UPDATE) {
			fprintf (stderr, "CTR ");
			DecreaseCoverageGraphCenter(g, bestv, pairs, curcover, dt, covercount, covered, LABEL_GREEDY);
			//DecreaseCoverageMatrixCenter(g, bestv, pairs, curcover, dt, covercount, covered, LABEL_GREEDY);
			//DecreaseCoverageGraphCenter (g, u, bestv, curcover, targetlist, dt, covercount, tdata[tid].coveragebfs, covered, triples);
		}
		delete [] tdata;
	}

	void OutputDetailedLabel(FILE *file, Graph &g, Labels &flabels, Labels &blabels, int v, vector<int> &rank) {		
		vector<LabelEntry> &label = flabels.labels[v];

		int hcount = label.size();
		//fprintf (stderr, "%d:%d:%d ", v, g.GetDegree(v), hcount);

		for (int i=0; i<hcount; i++) {
			int h = label[i].hub;
			int dvh = label[i].dist;

			int rh = rank[h]; //that's the rank of h
			
			int p = v;
			// find the parent
			for (int j=0; j<hcount; j++) {
				int x = label[j].hub;
				if (x == h) continue; 
				if (rank[x] < rank[p]) {
					int dvx = label[j].dist;
					int dxh = Query(flabels.labels[x], blabels.labels[h]);
					if (dvx + dxh == dvh) {p = x;}
				}
			}

			if (v==h) p = 2147483647;

			//v hub distance parent shortcut
			fprintf (file, "%d %d %d %d %d\n", v, h, dvh, p, 0);

		}

	}

	// weighted: take sum of the degrees of the neighbors
	// adaptive: once a vertex is picked, adjust priorities of neighbors
	// neighbors: not supported yet 
	void ComputeSimpleOrders (char *prefix, Graph &g, bool weighted, bool adaptive, int neighbors) {
		int n = g.VertexCount();

		vector<double> origpriority (n+1,0);
		BinaryHeap<double> heap (n);
		vector<int> rank2vertex(n+1,-1);

		for (int v=1; v<=n; v++) {
			double priority;

			if (weighted) {
				double sumneighbors = 0;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
					sumneighbors += g.GetDegree(w);
				}
				priority = sumneighbors;
			} else {
				priority = g.GetDegree(v);
			}
			heap.Insert(v,-priority);
			origpriority[v] = priority;
		}

		int rank = 0;
		while (!heap.IsEmpty()) {
			double value;
			unsigned int v;
			heap.RemoveFirst(v, value);
			value = -value;
			//fprintf (stderr, "%.0f:%.0f ", abs(value), origpriority[v]);
			rank ++;
			rank2vertex[rank] = v; 
			if (adaptive) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;
					if (heap.Contains(w)) {
						double oldvalue = heap.GetElementValue(w);
						heap.RemoveElement(w);
						double delta = weighted ? g.GetDegree(v) : 1.0; 
						double newvalue = oldvalue + delta;
						//fprintf (stderr, "<%.0f> ", newvalue);
						heap.Insert(w, newvalue);

					}
				}
			}
		}

		char filename[256];
		sprintf (filename, "%s-%s-%s.order", prefix, weighted ? "weighted" : "plain", adaptive ? "adaptive" : "fixed");
		OutputOrder(rank2vertex, n, filename);
	}

	void OutputOrder (vector<int> &rank2vertex, int n, char *filename) {
		FILE *file = fopen (filename, "w");
		fprintf (stderr, "\nOutputting order to file %s... ", filename);
		if (!file) {
			fprintf (stderr, "Could not open order file for writing.\n");
		} else {
			for (int r=1; r<=n; r++) {
				fprintf (file, "%d\n", rank2vertex[r]);
			}
			fclose(file);
		}
		fprintf (stderr, "done.\n");
	}




	void Solve (int argc, char **argv) {
		fprintf (stderr, "Should be solving labels.\n");


		bool LABEL_GREEDY = true; //false;
		bool PERTURB_EDGES = false;
		const int INFINITY = 1000000000;

		Graph g;
		if (argc <= 2 || argc % 2 != 1) {
			fprintf (stderr, "invalid number of parameters");
			exit(-1);
		}

		bool DIMACS_FORMAT = false;
		char *filename = NULL;
		char *prefix = NULL;
		int threads = 1;
		int seed = 1;
		//bool PERTURB_EGDGES = false;

		for (int i=1; i<argc; i+=2) {
			char *key = argv[i];
			char *value = argv[i+1];

			if (strcmp(key, "-dimacs")==0) {
				DIMACS_FORMAT = true;
				filename = value;
			} else if (strcmp(key, "-metis")==0) {
				DIMACS_FORMAT = false;
				filename = value;
			} else if (strcmp(key, "-order")==0) {
				prefix = value;
			} else if (strcmp(key, "-seed")==0)  {
				seed = atoi(value);
			} else if (strcmp(key, "-threads")==0)  {
				threads = atoi(value);
			} else if (strcmp(key, "-pert")==0) {
				PERTURB_EDGES = (atoi(value)!=0);
			} else if (strcmp(key, "-pathgreedy")==0) {
				LABEL_GREEDY = (atoi(value)==0);
			} else {
				fprintf (stderr, "Parameter [%s] not recognized.\n", key);
			}
		}
		RFWRandom::randomize(seed);

		if (DIMACS_FORMAT) {
			g.ReadDimacs(filename, PERTURB_EDGES);
		} else {
			g.ReadMETIS(filename, PERTURB_EDGES);
		}

		if (prefix != NULL) {
			ComputeSimpleOrders(prefix, g, false, false, 0);
			ComputeSimpleOrders(prefix, g, false, true, 0);
			ComputeSimpleOrders(prefix, g, true, false, 0);
			ComputeSimpleOrders(prefix, g, true, true, 0);
			exit(-1);
		}
		/*
		if (PERTURB_EDGES) {
			int n = g.VertexCount(); 
			vector<int> old2new(n+1,0);
			for (int v=1; v<=n; v++) old2new[v] = v;
			FILE *outgraph = fopen ("pert.graph", "w");
			OutputMETIS (outgraph, g, old2new); //WARNING: NOT WEIGHTED
			fclose(outgraph);
		}*/
		fprintf (stderr, "Read graph with %d vertices and %d edges.\n", g.VertexCount(), g.EdgeCount());
		//exit(-1);
		//ExtractGiantComponent(g);
		fflush(stderr);
		//exit(-1);
		//exit(-1);
		//g.ReadMETIS(argv[1], PERTURB_EDGES);
		//g.ReadSTP(argv[1]);

		//int maxthreads = -1;
		int maxthreads = threads; //atoi(argv[2]);
		SetActualThreads(maxthreads,true);

		fprintf (stderr, "Should be using %d threads.\n", maxthreads);
		
		fprintf (stderr, "perturb %d\n", PERTURB_EDGES);
		fprintf (stderr, "pathgreedy %d\n", !LABEL_GREEDY);
		fprintf (stderr, "threads %d\n", maxthreads);
		fprintf (stderr, "seed %d\n", seed);


		//DensestSubgraph(g);
		//exit(-1);

		double timedist, timecover, timeadd;
		int n = g.VertexCount();
		CoverCounters covercount(n); //statistics about what paths each candidate hub covers

		// now we have the distance table among everybody
		fprintf (stderr, "Creating matrix...\n");
		fflush (stderr);
		Matrix<int> dt(n); //distances

		fprintf (stderr, "Matrix created...\n");
		fflush (stderr);
		RFWTimer timer (true);


#pragma omp parallel for schedule(dynamic)
		for (int s=1; s<=n; s++) {
			Dijkstra(g,s,dt.dist[s]);
			double t = timer.getTime();
			if (s % 100 == 0) {
				fprintf (stderr, "%.2f:%.2f ", t, (double)n * (double)t / (double)s);
				fflush(stderr);
			}
		}
		timedist = timer.getTime();
		fprintf (stderr, "Distance table created in %.2f milliseconds.\n", 1000.0*timedist);
		timer.start();

		//dt.Output(stderr);
		fprintf (stderr, "Creating labels... (%.2f)\n", 1000 * timer.getTime());
		fflush (stderr);
		Labels flabels(n);
		Labels blabels(n);

		double totalcovered = 0;
		Matrix<int> covered(n);
		covered.Reset(false);
		for (int u=1; u<=n; u++) {
			for (int w=1; w<=n; w++) {
				if (dt.Get(u,w) >= INFINITY) {covered.Set(u,w,n+1);}
				else totalcovered ++;
			}
		}


		vector<int> rank(n+1,-1); //rank[v]: rank of vertex v
		vector<int> rank2vertex(n+1,1);

		//BinaryHeap<int> heap(n); //cou
		vector<int> curcover(n+1,0);

		//MatrixComputeCoverage(g, curcover, covered,  dt, LABEL_GREEDY, &covercount);
		GraphComputeCoverageCenter(g, curcover, covered,  dt, LABEL_GREEDY, &covercount);
		//FancyComputeCoverage(g, curcover, covered,  dt, LABEL_GREEDY, &covercount);
		timecover = timer.getTime();
		fprintf (stderr, "Initial coverage took %.3f seconds.\n", timecover);

		double sumcover = 0;
		for (int v=1; v<=n; v++) {
			sumcover += (double)(curcover[v]);
		}

		fprintf (stderr, "The average path is covered by %.4f vertices.\n", sumcover / ((double)n * (double)n));
		fprintf (stderr, "\n");
		fflush(stderr);
		//fprintf (stderr, "Exiting.\n");
		//exit(-1);

		timer.start();
		vector<int> forward(n+1,false);
		vector<int> backward(n+1,false);

		double totalsize = 0;


		int hubcount = 0;
		while (1) {

			//---------------------------------
			// pick vertex with the best score
			//---------------------------------
			int bestv = -1;
			double bestscore = 0;
			for (int v=1; v<=n; v++) {
				if (curcover[v] == 0) continue;
				double curscore = curcover[v];
				if (LABEL_GREEDY) {
					curscore /= (double)(covercount.backcount[v] + covercount.fwdcount[v]);
				} 

				//curscore = g.GetDegree(v);

				if (bestv<0 || curscore > bestscore) {
					bestv = v;
					bestscore = curscore;
				}
			}
			if (bestv<0) {
				fprintf (stderr, "Done with all vertices.\n");
				fflush(stderr);
				break;
			}
			hubcount ++;
			rank[bestv] = hubcount;
			rank2vertex[hubcount] = bestv;
			fprintf (stderr, "Hub %d (%d) covers %d paths.", hubcount, bestv, curcover[bestv]);
			int origcount = curcover[bestv];
			int vcount = 0;

			if (LABEL_GREEDY) {fprintf (stderr, "[%d:%d] ", covercount.fwdcount[bestv], covercount.backcount[bestv]);}
			for (int v=0; v<=n; v++) {forward[v] = backward[v] = false;}

			int addcount = 0;

			vector<int> targetlist;
			targetlist.reserve(n+1);

			BFSData coveragebfs(n);

			RFWTimer localtimer(true);
			int timerstep = 500000;
			int timertarget = timerstep;

			//-------------------------------------------
			// identify all uncovered pairs hit by bestv
			//-------------------------------------------
			bool FROM_SCRATCH = true;
			double target =(double)n*(double)n;
			fprintf (stderr, "[%.0f pairs] ", target);
			FROM_SCRATCH = ((double)origcount / target > 0.1);
			//FROM_SCRATCH = false;

			fprintf (stderr, "s%d ", FROM_SCRATCH);

			UpdateCoverageInfo(g, bestv, curcover, covered, dt, forward, backward, LABEL_GREEDY, &covercount, flabels, blabels, addcount, totalcovered, vcount, origcount, !FROM_SCRATCH);

			if (FROM_SCRATCH) {
				for (int x=1; x<=n; x++) curcover[x] = 0;
				if (LABEL_GREEDY) covercount.Reset();
				GraphComputeCoverageCenter(g, curcover, covered,  dt, LABEL_GREEDY, &covercount);
			}

#if 0
			//#pragma omp parallel for schedule(dynamic)
			for (int u=1; u<=n; u++) {
				if (covered.Get(u,bestv)) continue; //speedup! already covered this part 
				int duv = dt.Get(u,bestv);

				targetlist.clear();
				// do we really need to loop over everything to mark?
				for (int w=1; w<=n; w++) {

					//fprintf (stderr, "<%d,%d> ", u, w);
					if (covered.Get(u,w)) continue; // fprintf (stderr, "!");
					//fflush(stderr);
					if (dt.Get(u,w) == duv + dt.Get(bestv,w)) {
						if (!forward[u]) {
							forward[u] = true; 
							flabels.AddEntry(u,bestv,dt.Get(u,bestv));
							addcount ++;
						}
						if (!backward[w]) {
							backward[w] = true; 
							blabels.AddEntry(w,bestv,dt.Get(bestv,w));
							addcount ++;

						}

						if (LABEL_GREEDY) {covercount.RemovePath(u,bestv,w);}

						targetlist.push_back(w);
						// decrease conters of everybody that covers this path

						// could save something here by having two (pre-)sorted lists at x
						// the first has distances to x, the second distances from x
						// we can do binary search for the first match, then just walk within the desired range;
						// could delete entries as paths are covered
						
						
						/*
						bool UPDATE_COUNT = true;
						if (UPDATE_COUNT) {
							int duw = dt.Get(u,w); //must find all x that hit the u-w path
							for (int x=1; x<=n; x++) {
								if (duw == dt.Get(u,x) + dt.Get(x,w)) {
									if (curcover[x] == 0) continue; //x doesn't cover anything anyway
									if (x==bestv) continue;
									curcover[x] --;
									if (LABEL_GREEDY) {covercount.RemovePath(u,x,w);}
								}
							}
						}*/
						totalcovered --;
						vcount ++;
						covered.Set(u,w,bestv);
					}
				}

				// now remember all paths starting at u

				bool SIMPLE = false;
				if (SIMPLE) {
					int tcount = targetlist.size();
					for (int i=0; i<tcount; i++) {
						int w = targetlist[i];
						int duw = dt.Get(u,w); //must find all x that hit the u-w path
						for (int x=1; x<=n; x++) {
							if (duw == dt.Get(u,x) + dt.Get(x,w)) {
								if (curcover[x] == 0) continue; //x doesn't cover anything anyway
								if (x==bestv) continue;
								curcover[x] --;
								if (LABEL_GREEDY) {covercount.RemovePath(u,x,w);}
							}
						}
					}
				} else {
					//DecreaseCoverageMatrix(n, u, curcover, bestv, targetlist, dt, &covercount);
					DecreaseCoverageGraph(g, u, bestv, curcover, targetlist, dt, &covercount, coveragebfs, covered);


					
					// for fixed u and considering all w such that u-w is newly convered,
					// we have a path u-w.
					// For each such path, we must find all other elements x that hit u-w as well
					/*
					int tcount = targetlist.size();
					for (int x=1; x<=n; x++) { // try all possible middle vertices
						if (curcover[x] == 0) continue; //vertex no longer relevant
						if (x==bestv) continue;
						int dux = dt.Get(u,x); // prefix

						for (int i=0; i<tcount; i++) {
							int w = targetlist[i];
							if (dt.Get(u,w) == dux + dt.Get(x,w)) {
								curcover[x] --;
								if (LABEL_GREEDY) {covercount.RemovePath(u,x,w);}
							}
						}
					}*/
				}

				if (vcount > timertarget) {
					while (vcount > timertarget) {
						timertarget += timerstep;
					}
					double t = localtimer.getTime();
					fprintf (stderr, "%.0f/%.0f ", t, (double)origcount * t / (double)vcount); 
					fflush(stderr);
				}

				//fflush(stderr);
			}
#endif
			totalsize += addcount;

			fprintf (stderr, "Counted %d; added %d label entries; %.2f pct uncovered; %.2f seconds\n", vcount, addcount, 100 * totalcovered / ((double)n*(double)n), timer.getTime());
			fflush(stderr);
			curcover[bestv] = 0;
			//fflush(stderr);
		}
		timeadd = timer.getTime();

		fprintf (stderr, "time_distance %.2f\n", timedist);
		fprintf (stderr, "time_cover %.2f\n", timecover);
		fprintf (stderr, "time_add %.2f\n", timeadd);
		fprintf (stderr, "time_total %.2f\n", timedist + timecover + timeadd);
		fprintf (stderr, "Average label size is %.2f\n", (double)totalsize / (double)(2*n)); 
		fflush(stderr);
		fflush(stdout);

		//if (LABEL_GREEDY) {covercount.Output(stderr);}
		bool OUTPUT_LABELS = true;
		if (OUTPUT_LABELS) {
			fprintf (stderr, "Outputting labels... ");
			blabels.sort();
			flabels.sort();
			FILE *file = fopen ("output.labels", "w");
			if (!file) {
				fprintf (stderr, "Could not open labels file for writing.\n");
			} else {
				for (int v=1; v<=n; v++) {
					OutputDetailedLabel(file, g, flabels, blabels, v, rank);
				}
				fclose(file);
			}
			fprintf (stderr, "done.\n");
		}

		bool OUTPUT_ORDER = true;
		if (OUTPUT_ORDER) {
			OutputOrder(rank2vertex, n, "output.order");
			/*
			FILE *file = fopen ("output.order", "w");
			if (!file) {
				fprintf (stderr, "Could not open labels file for writing.\n");
			} else {
				for (int r=1; r<=n; r++) {
					fprintf (file, "%d\n", rank2vertex[r]);
				}
				fclose(file);
			}
			fprintf (stderr, "done.\n");*/
		}



		//flabels.Output();
		bool CHECK_LABELS = true;
		if (CHECK_LABELS) {
			fprintf (stderr, "Checking correctness of all labels...\n");
			fflush(stderr);
			blabels.sort();
			flabels.sort();
			for (int v=1; v<=n; v++) {
				for (int w=1; w<=n; w++) {
					int labdist = Query(flabels.labels[v], blabels.labels[w]);
					int dist = dt.Get(v,w);
					//fprintf (stderr, "(%d,%d): %d,%d\n", v, w, labdist, dist);
					if (dist != labdist) {
						flabels.OutputLabel(v);
						blabels.OutputLabel(w);
						exit(-1);
					}
				}
			}
			fprintf (stderr, "Labels are valid.\n");
			fflush(stderr);
		}


		// must compute the 

	}

	int Query(vector<LabelEntry> &flabel, vector<LabelEntry> &blabel) {
		int bestval = 1000000000;
		int i = 0;
		int j = 0;
		int fsize = flabel.size();
		int bsize = blabel.size();

		while (1) {
			if (i == fsize) break;
			if (j == bsize) break;

			//fprintf (stderr, "%d:%d ", flabel[i].hub, blabel[j].hub);

			int diff = flabel[i].hub - blabel[j].hub;

			if (diff == 0) {
				//fprintf (stderr, "m%d ", flabel[i].hub);
				int d = flabel[i].dist + blabel[j].dist;
				if (d < bestval) {bestval = d;}
				i++; j++;
			} else {
				if (diff > 0) {j++;}
				else {i++;}
			}
		}
		return bestval;
	}

	void Dijkstra (Graph &g, int s, int *dist) {
		bool verbose = false;
		if (verbose) fprintf (stderr, "Dij "); //RUNNING DIJKSTRA!\n");
		int n = g.VertexCount();
		BinaryHeap<EdgeCost> heap(n); // = new BinaryHeap<ArcCost>(n);
		int INFINITY = 1000000000;
		for (int v=0; v<=n; v++) dist[v] = INFINITY;
		int nscanned = 0;

		heap.Insert(s, 0);
		dist[s] = 0;

		while (!heap.IsEmpty()) {
			unsigned int v;
			EdgeCost vdist;
			heap.RemoveFirst(v,vdist); //v, out acost);

			//scan outgoing arcs
            nscanned ++;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				int w = a->head; //neighbor
				int dw = vdist + a->cost;
				if (dw < dist[w]) {
					dist[w] = dw;
					heap.Insert(w, dw);
				}
            }
        }


		if (verbose) {
			//if (nscanned != n) {
				fprintf (stderr, "Scanned %d/%d vertices.\n", nscanned, n);
				fflush(stderr);
			//}
		}
	}




};
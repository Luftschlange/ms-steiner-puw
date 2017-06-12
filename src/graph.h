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
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include "binheap.h"
#include "rfw_random.h"
using std::vector;
typedef double EdgeCost;
//THE LINE BELOW WAS SET BY THOMAS
#define EDGE_COST_PRECISION 0.0000001 
//#define EDGE_COST_PRECISION 0.0001
#define INFINITE_COST 1e32
#define DEBUG_VERTEX -1


class Coordinates {
private:
	int maxid;
	int dim;
	vector <double> coord;

	void fatal (const string &msg) {
		fprintf (stderr, "ERROR::Coordinates: %s.\n", msg.c_str());
		exit(-1);
	}

	void AllocateCoords() {
		if (CoordinatesAvailable()) fatal ("cannot allocate coordinates twice");
		if (dim>10) fatal ("too many dimensions");
		//fprintf (stderr, "Allocating coordinates of size %d on instance %p.\n", maxid, this);
		coord.resize(dim*(maxid+1));
		//coord = new double[];
	}

	bool CheckRange(int v) {
		if (v<0 || v>maxid) {
			fprintf (stderr, "element out of range");
			return false;
		} else {
			return true;
		}
	}

public:
	inline void SetCoordinates(int v, double x, double y) {
		if (!CoordinatesAvailable()) AllocateCoords();
		if (CheckRange(v)) {
			coord[2*v] = x;
			coord[2*v+1] = y;
		}
	}
	inline void GetCoordinates(int v, double &x, double &y) {
		if (!CoordinatesAvailable()) {x = y = 0;}
		else {
			CheckRange(v);
			//fprintf (stderr, "<%p:%d> ", this, v);
			//fflush(stderr);
			x = coord[2*v];
			y = coord[2*v+1];
		}
	}

	Coordinates() {
		//coord = NULL;
		//fprintf (stderr, "Created coordinates with size %d.\n", coord.size());
		//fflush(stderr);
		dim = 2;
		maxid = 0;
	}

	void SetMaxId(int n) {
		maxid = n;
	}

	~Coordinates() {
		//delete [] coord;
	}

	inline bool CoordinatesAvailable() {
		return (coord.size() != 0);
		//return (coord!=NULL);
	}
};

class EdgeDescriptor {
public:
	int v;
	int w;
	EdgeCost cost;

	EdgeDescriptor (int _v, int _w, EdgeCost _cost) {
		if (_v < _w) {
			v = _v;
			w = _w;
		} else {
			v = _w;
			w = _v;
		}
		cost = _cost;
	}
};

class GraphDescriptor {
private:
	int m;

public:
	//int m;
	int n;
	//int t;
	int tcount;
	vector<EdgeDescriptor> edges;
	vector<bool> terminal;
	double fixed;

	GraphDescriptor () {
		// add dummy edge at position 0; edge IDs start at 1
		edges.push_back(EdgeDescriptor(-1,-1,-1));
		n = m = -1;
		tcount = 0;
		fixed = 0;
	}

	inline void SetFixed(double f) {fixed = f;}
	inline void IncFixed(double f) {fixed += f;}


	inline int GetEdgeCount() const {
		return (int)(edges.size() - 1);//m;
	}

	inline void CheckRange (int v) {
		if (v<=0 || v>n) {
			fprintf (stderr, "Vertex %d is out of range [1,%d].\n", v, n);
			exit(-1);
		}
	}

	void SetVertices(int _n) {n = _n; terminal.resize(n+1,false);}
	void SetEdges(int _m) {
		//m = _m;}
		//fprintf (stderr, "Ignoring edge counts.\n");
	}

		//void SetTerminals(int _t) {t = _t;}
	int AddEdge(int v, int w, EdgeCost cost) {
		int cursize = (int)edges.size();
		//if (cursize >= (m+1)) {fprintf (stderr, "too many edges"); exit(-1);}
		CheckRange(v);
		CheckRange(w);
		edges.push_back(EdgeDescriptor(v,w,cost));
		return cursize;
	}
	void MakeTerminal(int v) {CheckRange(v); if (!terminal[v]) {terminal[v]=true; tcount++;}}
	void UnmakeTerminal(int v) {
		CheckRange(v); if (terminal[v]) {terminal[v] = false; tcount--;}
	}

	int TerminalCount() {return tcount;}
};


struct SPGArc {
	int label; //arc label
	int head;
	EdgeCost cost;

	inline void Set (int h, EdgeCost c, int l) {head = h; cost = c; label = l;}
};


// dynamic graph, with new edges added as they are created
class DynamicGraph {
	int m, n, tcount;
	vector <SPGArc> arclist; //this will not change
	vector <int> first;
	vector <int> last;
	GraphDescriptor gd;

	void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		exit(-1);
	}

public:
	Coordinates coord;
	void OutputAdjacencyList (FILE *file, int v) {
		fprintf (stderr, "%4d:", v);
		SPGArc *a, *end;
		for (GetBounds(v,a,end); a<end; a++) {
			fprintf (stderr, " %d:%.0f:%d ", a->head, a->cost, a->label);
		}
		fprintf (stderr, "\n");
	}
	
	void OutputGraph (FILE *file) {
		for (int v=1; v<=n; v++) OutputAdjacencyList(file, v);
	} 

	inline double GetFixedCost() {return gd.fixed;}
	inline void SetFixedCost(EdgeCost f) {gd.SetFixed(f);}
	inline void IncFixedCost(EdgeCost f) {gd.IncFixed(f);}
	inline EdgeCost GetCost (int e) const {return gd.edges[e].cost;}
	inline int VertexCount() const {return n;}
	inline int EdgeCount() const {return gd.GetEdgeCount();}
	inline int TerminalCount() const {return gd.tcount;}
	inline int GetDegree(int v) const {return last[v]-first[v];}
	inline int GetReservedDegree(int v) const {return first[v+1] - first[v];}
	inline int GetFirstEdgeLabel(int v) const { 
		return arclist[first[v]].label;
	}

	inline bool IsTerminal(int v) const {return gd.terminal[v];}

	inline void GetBounds (int v, SPGArc * &start, SPGArc * &end) {
		start = &arclist[first[v]];
		end = &arclist[last[v]];
	}
	inline void SetVertices(int _n) {gd.SetVertices(_n);}
	inline void SetEdges(int _m) {gd.SetEdges(_m);} //is this needed
	inline int BatchAddEdge(int v, int w, EdgeCost cost) {return gd.AddEdge(v,w,cost);}
	inline void MakeTerminal(int v) {gd.MakeTerminal(v);}
	inline void UnmakeTerminal(int v) {gd.UnmakeTerminal(v);}
	inline void GetEndpoints (int e, int &v, int &w) {
		v = gd.edges[e].v;
		w = gd.edges[e].w;
	}

	// insert an edge on a committed graph, acutally changing
	// the appropriate adjacency lists
	inline int InsertEdge (int _v, int _w, EdgeCost cost) {
		int newid = gd.AddEdge(_v,_w,cost);
		const bool verbose = false;

		int endpoint[2] = {_v,_w};
		if (verbose) fprintf (stderr, "Inserting (%d,%d):%.0f\n", _v, _w, cost);

		// insert new arc at the back of each adjacency list
		for (int i=0; i<2; i++) {
			int v = endpoint[i];
			int w = endpoint[1-i];
			if (verbose) fprintf (stderr, "Before: ");
			if (verbose) OutputAdjacencyList(stderr, v);
			if (last[v] >= first[v+1]) fatal ("not enough allocated space for new arc");
			arclist[last[v]].Set(w,cost,newid);
			last[v]++;

			if (verbose) fprintf (stderr, "After: ");
			if (verbose) OutputAdjacencyList(stderr, v);
		}

		if (verbose) fprintf (stderr, "\n");
		return newid;
	}

	inline void HideEdge (int e) {
		const bool verbose = false;
		int neighbor[2];
		GetEndpoints(e,neighbor[0],neighbor[1]);
		if (neighbor[0] == neighbor[1]) {
			fatal ("found self-loop");
		}
		for (int i=0; i<2; i++) {
			int v = neighbor[i];
			if (verbose) fprintf (stderr, "Before: ");
			if (verbose) OutputAdjacencyList(stderr, v);
			int p;
			for (p = first[v]; p<last[v]; p++) {
				if (arclist[p].label == e) break;
			}
			if (p == last[v]) fatal ("tried to remove nonexistent edge");
			// now overwrite stuff

			last[v] --;
			for ( ; p<last[v]; p++) {arclist[p] = arclist[p+1];}
			if (verbose) fprintf (stderr, "After: ");
			if (verbose) OutputAdjacencyList(stderr, v);
		}
		if (verbose) fprintf (stderr, "Done hiding edge %d=(%d,%d).\n", e, neighbor[0], neighbor[1]);
	}

	//--------------------------------------------------
	// Create adjacency lists from plain lists of edges
	//--------------------------------------------------
	void Commit () {
		const bool verbose = false;
		n = gd.n;
		int ecount = (int)gd.edges.size() - 1; //we start at 1...
		
		/*
		if (ecount!=gd.m) {
			gd.m = ecount;
			//fprintf (stderr, "Graph actually has %d edges.\n", gd.m);
		}*/

		m = gd.GetEdgeCount();

		//fprintf (stderr, "Creating subgraph with %d edges.\n", m);
		tcount = gd.tcount;
		arclist.resize(2*m+1); //adding 1 to avoid stupid vector out of range in debug mode
		first.resize(n+2,0);
		last.resize(n+1,0);

		int v, e;
		// now first will correspond to the degree
		for (e=1; e<=m; e++) {
			first[gd.edges[e].v]++;
			first[gd.edges[e].w]++;
		}

		// convert degrees into first positions
		int acc = 0; //accumulated degree up to this point
		for (v=1; v<=n; v++) {
			if (verbose) fprintf (stderr, " %4d", first[v]);
			int temp = first[v];
			first[v] = acc;
			acc += temp;
			if (verbose) fprintf (stderr, " %4d\n", first[v]);
		}

		// insert the actual arcs
		for (e=1; e<=m; e++) {
			int v = gd.edges[e].v;
			int w = gd.edges[e].w;
			EdgeCost cost = gd.edges[e].cost;
			arclist[first[v]++].Set(w,cost,e);
			arclist[first[w]++].Set(v,cost,e);
		}

		// restore first pointers (including sentinel)
		for (v = n+1; v > 0; v--) {first[v] = first[v-1];}

		// saved last for now
		for (v=1; v<=n; v++) {last[v] = first[v+1];}

		if (verbose) OutputGraph (stderr);

		//for (v=0; v<=n+1; v++) {fprintf (stderr, "%d ", first[v]);}
		// now determine first outgoing vertex
		if (verbose) {
			fprintf (stderr, "Instance committed.\n");
			fflush(stderr);
		}
	}

};

class Graph {
private:
	int m, n, tcount;
	vector <SPGArc> arclist;
	vector <int> first;
	//vector <bool> terminal;
	GraphDescriptor gd;

	void fatal (const string &msg) {
		fprintf (stderr, "ERROR: %s.\n", msg.c_str());
		exit(-1);
	}

public:
	Coordinates coord;

	Graph () {
		n = m = tcount = 0;
		//arclist = NULL;
		//first = NULL;
	}


	inline double GetFixedCost() {return gd.fixed;}
	inline void SetFixedCost(double f) {gd.SetFixed(f);}
	inline void IncFixedCost(double f) {gd.IncFixed(f);}

	void ReadSTP (char *filename) {
		fprintf (stderr, "Should be reading STP file %s.\n", filename);
		FILE *file = fopen (filename, "r");
		if (!file) fatal ("could not open input file");
		ReadSTP(file);
		fclose (file);
	}

	void ReadDimacs (char *filename, bool disturb) {
		fprintf (stderr, "Should be reading dimacs file %s.\n", filename);
		FILE *file = fopen (filename, "r");
		if (!file) fatal ("could not open input file");
		ReadDimacs(file, disturb);
		fclose (file);
	}

	void ReadFile (char *filename, bool perturb) {
		fprintf (stderr, "Checking format of file %s.\n", filename);
		FILE *file = fopen (filename, "r");

		fclose(file);


	}

	void ReadMETIS (char *filename, bool perturb) {
		fprintf (stderr, "Should be reading METIS file %s.\n", filename);
		FILE *file = fopen (filename, "r");
		fprintf (stderr, "Opened file %d.\n", file);
		fflush(stderr);
		if (!file) fatal ("could not open input file");
		fflush(stderr);
		ReadMETIS(file, perturb);
		fclose (file);
	}

	inline void SetVertices(int _n) {gd.SetVertices(_n);}
	inline void SetEdges(int _m) {gd.SetEdges(_m);}
	//inline void SetTerminals(int _t) {gd.SetTerminals(_t);}
	inline int AddEdge(int v, int w, EdgeCost cost) {return gd.AddEdge(v,w,cost);}
	inline void MakeTerminal(int v) {gd.MakeTerminal(v);}
	inline void UnmakeTerminal(int v) {gd.UnmakeTerminal(v);}

	//--------------------------------------------------
	// Create adjacency lists from plain lists of edges
	//--------------------------------------------------
	void Commit () {
		const bool verbose = false;
		n = gd.n;
		int ecount = (int)gd.edges.size() - 1; //we start at 1...
		
		/*
		if (ecount!=gd.m) {
			gd.m = ecount;
			//fprintf (stderr, "Graph actually has %d edges.\n", gd.m);
		}*/

		m = gd.GetEdgeCount();

		//fprintf (stderr, "Creating subgraph with %d edges.\n", m);
		tcount = gd.tcount;
		arclist.resize(2*m+1); //adding 1 to avoid stupid vector out of range in debug mode
		first.resize(n+2,0);

		int v, e;
		// now first will correspond to the degree
		for (e=1; e<=m; e++) {
			first[gd.edges[e].v]++;
			first[gd.edges[e].w]++;
		}

		// convert degrees into first positions
		int acc = 0; //accumulated degree up to this point
		for (v=1; v<=n; v++) {
			if (verbose) fprintf (stderr, " %4d", first[v]);
			int temp = first[v];
			first[v] = acc;
			acc += temp;
			if (verbose) fprintf (stderr, " %4d\n", first[v]);
		}

		// insert the actual arcs
		for (e=1; e<=m; e++) {
			int v = gd.edges[e].v;
			int w = gd.edges[e].w;
			EdgeCost cost = gd.edges[e].cost;
			arclist[first[v]++].Set(w,cost,e);
			arclist[first[w]++].Set(v,cost,e);
		}

		// restore first pointers (including sentinel)
		for (v = n+1; v > 0; v--) {first[v] = first[v-1];}

		if (verbose) OutputGraph (stderr);

		//for (v=0; v<=n+1; v++) {fprintf (stderr, "%d ", first[v]);}
		// now determine first outgoing vertex
		if (verbose) {
			fprintf (stderr, "Instance committed.\n");
			fflush(stderr);
		}
	}

	inline int GetArcTail (int a) const {
		return (a<=m) ? GetFirstEndpoint(a) : GetSecondEndpoint(a-m);
	}

	inline void GetArcEndpoints(int a, int &v, int &w) {
		if (a <= m) {
			GetEndpoints(a,v,w); //[1..m]: same as undirected edge (smaller first)
		} else {
			GetEndpoints(a-m,w,v); //high: oppositive as undirected edge (higher first)
		}
	}

    inline int GetOutgoingLabel (int v, SPGArc *a) const {
		int w = a->head;
		return (v<w) ? a->label : a->label + m;
	}

	inline int GetIncomingLabel (int v, SPGArc *a) const {
		int w = a->head;
		return (v<w) ? a->label+m : a->label;
	}

    inline int GetIncomingLabel (int v, int w, int alabel) const {
		return (v<w) ? alabel + m : alabel;
	}



	inline int VertexCount() const {return n;}
	inline int EdgeCount() const {return m;}
	inline int TerminalCount() const {return gd.tcount;}
	inline int GetDegree(int v) const {return first[v+1]-first[v];}

	inline bool IsTerminal(int v) const {return gd.terminal[v];}

	inline void GetBounds (int v, SPGArc * &start, SPGArc * &end) {
		start = &arclist[first[v]];
		end = &arclist[first[v+1]];
	}

	// copy current edge costs to "costs"
	void RetrieveCosts (vector<EdgeCost> &costs) {
		if (costs.size() != m+1) fatal ("invalid array in cost retrieval");

		costs.clear();
		costs.resize(m+1);
		int esize = (int)gd.edges.size();
		for (int e=0; e<=m; e++) {
			costs[e] = gd.edges[e].cost;
		}
	}

	// update all costs in the graph according 
	void ApplyCosts (vector<EdgeCost> &newcost) {
		if (newcost.size() != m+1) fatal ("invalid array in apply costs");
		for (int v=1; v<=n; v++) {
			SPGArc *a, *end;
			for (GetBounds(v,a,end); a<end; a++) {
				a->cost = newcost[a->label];
			}
		}
		for (int e=0; e<=m; e++) {
			gd.edges[e].cost = newcost[e];
		}
	}

	void OutputGraph (FILE *file) {
		for (int v=1; v<=n; v++) {
			fprintf (stderr, "%4d:", v);
			SPGArc *a, *end;
			for (GetBounds(v,a,end); a<end; a++) {
				fprintf (stderr, " %d:%d:%d ", a->head, a->cost, a->label);
			}
			fprintf (stderr, "\n");
		}		
	}

	void ReadDimacs (FILE *file, bool perturb) {
		if (gd.n>=0) fatal ("cannot read file twice");
		const bool verbose = false;
		const int BUFSIZE = 65534; //1048574;
		char buffer [BUFSIZE+2];
		int n, m, t, tcount, ecount;
		t = n = m = -1;
		tcount = ecount = 0;
		int tdeclared = -1;
		int selfedges = 0;

		int curvertex = -1;
		int a, b, c;

		fprintf (stderr, "About to start reading file in DIMACS format.\n");
		fflush (stderr);

		while (fgets(buffer, BUFSIZE, file)!=0) {
			//fprintf (stderr, ">> %s", buffer);
			if (n < 0) {
				if (sscanf (buffer, "p sp %d %d", &a, &b)>0) {
					n = a;
					m = b;
					tdeclared = 0;
					fprintf (stderr, "Reading graph with %d vertices and %d edges.\n", n, m);
					fflush(stderr);
					curvertex = 0;
					SetVertices(n);
					SetEdges(m);
				}
			} else {
				curvertex ++;
				if (sscanf (buffer, "a %d %d %d", &a, &b, &c)>0) {
					if (a==b) selfedges ++; //fprintf (stderr, "Found self (%d,%d) edge.\n", a, b);
					if (a<b) {
						ecount ++;
						int weight = c;
						if (perturb) weight = 1000000 - RFWRandom::getInteger(1,100);
						AddEdge(a, b, weight);
					}
				}
			}
		}
		fprintf (stderr, "Done with main loop.\n");
		fflush(stderr);


		if (tdeclared != gd.tcount) {
			fprintf (stderr, "%d terminals declared, %d in descriptor\n", tdeclared, gd.tcount);
			fatal ("invalid number of terminals");
		}
		if (selfedges > 0) fprintf (stderr, "Found %d self edges.\n", selfedges);
		if (ecount != m) {
			fprintf (stderr, "ecount:%d m:%d\n", ecount, m);
			if (ecount + selfedges / 2 == m) {
				fprintf (stderr, "Redefining the number of edges in the graph.\n");
			} else {
				fatal ("invalid number of edges");
			}
		}
		fprintf (stderr, "Done reading compact description of instance (n=%d m=%d).\n", gd.n, gd.GetEdgeCount());
		fflush(stderr);
		Commit();
		//Init(gd);
	}

	void ReadMETIS (FILE *file, bool perturb) {
		if (gd.n>=0) fatal ("cannot read stp twice");

		const bool verbose = false;
		const int BUFSIZE = 65534; //1048574;
		char buffer [BUFSIZE+2];
		//int v, w, 
		//int cost;
		int n, m, t, tcount, ecount;
		//double x, y;
		t = n = m = -1;
		tcount = ecount = 0;
		int tdeclared = -1;

		int curvertex = -1;
		int a, b, c;

		fprintf (stderr, "About to start reading file.\n");
		fflush (stderr);

		while (fgets(buffer, BUFSIZE, file)!=0) {
			//fprintf (stderr, "!");
			//fflush(stderr);
			//fprintf (stderr, "%d > %s >\n", n, buffer);
			//continue;

			//fprintf (stderr, ">> %s", buffer);
			if (n < 0) {
				if (sscanf (buffer, "%d %d %d", &a, &b, &c)>0) {
					n = a;
					m = b;
					tdeclared = 0;
					if (c != 0) fatal ("Cannot read weighted metis graphs.\n");
					fprintf (stderr, "Reading graph with %d vertices and %d edges.\n", n, m);
					fflush(stderr);
					//exit(-1);
					curvertex = 0;
					SetVertices(n);
					SetEdges(m);
				}
			} else {
				curvertex ++;
				char *token = strtok (buffer, "\t\n ");
				while (token) {
					int w = atoi(token);
					if (curvertex < w) {
						ecount ++;
						int weight = 1;
						if (perturb) weight = 10000000 + RFWRandom::getInteger(1,100000);
						AddEdge(curvertex, w, weight);
						//fprintf (stderr, "(%d,%d) ", curvertex,w);
					}
					token = strtok(NULL, "\t\n ");
				}
				//fflush(stderr);
			}
			/*
			if (sscanf (buffer, "E %d %d %d", &v, &w, &cost)>0) {
				if (verbose) fprintf (stderr, "(%d,%d) has cost %d\n", v, w, cost); 
				ecount ++; 
				AddEdge(v,w,cost);
			} else if (sscanf (buffer, "T %d", &v)>0) {
				if (verbose) fprintf (stderr, "%d is a terminal.\n", v); 
				tcount ++; 
				MakeTerminal(v);
			} else if (sscanf (buffer, "DD %d %lg %lg", &v, &x, &y)>0) {
				if (verbose) fprintf (stderr, "Vertex %d has coordinates (%.0f,%.0f).\n", v, x, y); 
				coord.SetCoordinates(v, x, y);
			} else if (sscanf (buffer, "Terminals %d", &t)>0) {
				if (verbose) fprintf (stderr, "There are %d terminals.\n", t); 
				tdeclared = t;//SetTerminals(t);
			} else if (sscanf (buffer, "Nodes %d", &n)>0) {
				if (verbose) fprintf (stderr, "There are %d nodes.\n", n);
				coord.SetMaxId(n);
				SetVertices(n);
			} else if (sscanf (buffer, "Edges %d", &m)>0) {
				if (verbose) fprintf (stderr, "There are %d edges.\n", m); 
				SetEdges(m);
			}*/
		}

		if (tdeclared != gd.tcount) {
			fprintf (stderr, "%d terminals declared, %d in descriptor\n", tdeclared, gd.tcount);
			fatal ("invalid number of terminals");
		}
		if (ecount != m) {
			fprintf (stderr, "ecount:%d m:%d\n", ecount, m);
			fatal ("invalid number of edges");
		}
		fprintf (stderr, "Done reading compact description of instance (n=%d m=%d t=%d).\n", gd.n, gd.GetEdgeCount(), gd.tcount);
		fflush(stderr);
		Commit();
		fprintf (stderr, "Done committing.\n");
		fflush(stderr);
		//Init(gd);
	}




	void ReadSTP (FILE *file) {
		if (gd.n>=0) fatal ("cannot read stp twice");

		int zerocount = 0;
		const bool verbose = false;
		const int BUFSIZE = 2046;
		char buffer [BUFSIZE+2];
		char section [BUFSIZE+2];
		int v, w;
		double cost;
		int n, m, t, tcount, ecount;
		double x, y;
		t = n = m = -1;
		tcount = ecount = 0;
		int tdeclared = -1;
		double fixed; //total fixed cost
		bool READING_TERMINALS = true;
		bool FOUND_PRESOLVE = false;
		int SCALING_FACTOR = 1;
		if (SCALING_FACTOR != 1) fprintf (stderr, "USING SCALING FACTOR %d.\n", SCALING_FACTOR);
		bool KILL_ZERO = true;

		double mincost = 1.0 / 1000000;

		while (fgets(buffer, BUFSIZE, file)!=0) {
			//fprintf (stderr, ">> %s", buffer);
			if (FOUND_PRESOLVE) {
				if (sscanf(buffer, "Fixed %lg\n", &fixed) == 1) {
					fprintf (stderr, "Read fixed cost: %.3f\n", fixed);
					gd.IncFixed(fixed);
				}
				continue; //ignore edges and terminals in presolve section
			} 



			if ((sscanf (buffer, "E %d %d %lf", &v, &w, &cost)>0) || (sscanf (buffer, "A %d %d %lf", &v, &w, &cost)>0)) {
				if (verbose) fprintf (stderr, "(%d,%d) has cost %f\n", v, w, (double)cost); 
				//cost *= 1000000;
				ecount ++; 
				if (SCALING_FACTOR != 1) cost = (cost + SCALING_FACTOR - 1) /SCALING_FACTOR;
				if (cost == 0) {
					if (zerocount == 0) {
						fprintf(stderr, "Found zero edge!\n");
					}
					zerocount ++;
					cost = mincost;
				}
				const bool PERTURB = false;
				if (PERTURB) {cost -= 0.0001;}
				AddEdge(v,w,(EdgeCost)cost);
			} else if (READING_TERMINALS && sscanf (buffer, "T %d", &v)>0) {
				if (verbose) fprintf (stderr, "%d is a terminal.\n", v); 
				tcount ++; 
				MakeTerminal(v);
			} else if (sscanf (buffer, "DD %d %lg %lg", &v, &x, &y)>0) {
				if (verbose) fprintf (stderr, "Vertex %d has coordinates (%.0f,%.0f).\n", v, x, y); 
				coord.SetCoordinates(v, x, y);
			} else if (sscanf (buffer, "Terminals %d", &t)>0) {
				if (verbose) fprintf (stderr, "There are %d terminals.\n", t); 
				tdeclared = t;//SetTerminals(t);
			} else if (sscanf (buffer, "Nodes %d", &n)>0) {
				if (verbose) fprintf (stderr, "There are %d nodes.\n", n);
				coord.SetMaxId(n);
				SetVertices(n);
			} else if (sscanf (buffer, "Edges %d", &m)>0) {
				if (verbose) fprintf (stderr, "There are %d edges.\n", m); 
				SetEdges(m);
			} //else if (sscanf (buffer, "SECTION Terminals")==0) {
				//fprintf (stderr, "Matched.\n");
			else if (sscanf (buffer, "SECTION %s", section)==1) {
				if (strcmp(section, "Presolve")==0) {
					FOUND_PRESOLVE = true;
				}
			}
		}


		if (tdeclared != gd.tcount) fatal ("invalid number of terminals");
		if (ecount != m) {

			fprintf (stderr, "WARNING: expected m=%d edges, found only ecount=%d.\n", m, ecount);
			//fatal ("invalid number of edges");
		}
		if (zerocount > 0) fprintf (stderr, "Found %d zero edges.\n", zerocount);
		fprintf (stderr, "Done reading compact description of instance (n=%d m=%d t=%d).\n", gd.n, gd.GetEdgeCount(), gd.tcount); fflush(stderr);
		Commit();
		//fprintf (stderr, "Done committing.\n"); fflush(stderr);
	}

	const void CheckEdgeRange(int e) {
		if (e < 1 || e>m) fatal ("edge out of range");
	}

	inline EdgeCost GetCost (int e) const {
		//CheckEdgeRange(e);
		return gd.edges[e].cost;
		//fprintf (stderr, "Should be getting cost.\n");
		//return 0;
	}

	inline int GetFirstEndpoint(int e) const {return gd.edges[e].v;}
	inline int GetSecondEndpoint(int e) const {return gd.edges[e].w;} 


	inline EdgeCost GetMinCost() {
		if (this->EdgeCount() < 1) fatal ("Undefined minimum edge cost.\n");
		EdgeCost mincost = GetCost(1);
		for (int e=2; e<=m; e++) {
			if (GetCost(e) < mincost) mincost = GetCost(e);
		}
		return mincost;
	}

	inline void GetEndpoints (int e, int &v, int &w) {
		v = gd.edges[e].v;
		w = gd.edges[e].w;
	}

	inline void GetDirectedEndpoints(int a, int &v, int &w) {
		if (a > m) {
			a -= m;
            v = gd.edges[a].w;
            w = gd.edges[a].v;
        } else {
			v = gd.edges[a].v;
            w = gd.edges[a].w;
			//fprintf (stderr, "<%d,%d> ", v, w); fflush(stderr);
        }
    }


	inline int GetOther (int e, int v) {
		int a, b;
		GetEndpoints(e,a,b);
		return (v==a) ? b : a;
	}
};

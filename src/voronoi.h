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

/*----------------------------------------------------
 | Information necessary to maintain Voronoi diagrams
 *---------------------------------------------------*/


#include "graph.h"

struct VoronoiUnit {
public:
	int vorbase;   //base of voronoi region
    int parc;      //label of the parent edge (WARNING: CHANGE NAME TO PEDGE)
    EdgeCost dist; //distance from the base

	inline void Reset() { 
		vorbase = parc = 0;
		dist = 0; 
	}

    inline void Reset(EdgeCost d) {
		vorbase = parc = 0;
        dist = d;
	}

    inline void Update(int b, int p, EdgeCost d) {
		vorbase = b;
        parc = p;
        dist = d;
	}

    inline void CopyFrom(VoronoiUnit x) {
		vorbase = x.vorbase;
        parc = x.parc;
        dist = x.dist;
    }
};


class VoronoiData {
private:
	VoronoiUnit *vor;
    int n;

public:
    void Reset() {
		for (int v = 0; v <= n; v++) {vor[v].Reset();}
	}

    void CopyFrom(int v, VoronoiData &other) {
		Update(v, other.GetBase(v), other.GetParentArc(v), other.GetDistance(v));
	}

    void Clone(VoronoiData &other) {
		for (int v=0; v<=n; v++) CopyFrom(v, other);
	}

    inline void Reset(int v) {vor[v].Reset();}

    VoronoiData(int _n) {
		n = _n;
        vor = new VoronoiUnit[n + 1];
        Reset();
    }

	~VoronoiData() {delete [] vor;}

    inline void Update(int v, int b, int p, EdgeCost d) {vor[v].Update(b, p, d);}

	// create voronoi region for v
    inline void MakeBase(int v) {vor[v].Update(v,0,0);}

	// query functions
	inline int GetBase(int v) const {return vor[v].vorbase;}
    inline EdgeCost GetDistance(int v) const {return vor[v].dist;}
	inline int GetParentArc(int v) const {return vor[v].parc;}
};


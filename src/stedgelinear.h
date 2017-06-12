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

/// <summary>
/// Linear-time implementation of link-cut trees with costs on edges.
/// </summary>
template <class CostType> class STEdgeLinear {
private:
	struct Node {
        CostType cost;
        int parent;

        void SetParentCost(int p, CostType c) {
            cost = c;
            parent = p;
        }
    };

    int n; //maximum item id
    int m; //edge count
    Node* node;

	void fatal (const string &msg) {
		fprintf (stderr, "STEdgeLinear::Error: %s\n", msg.c_str());
		exit(-1);
	}

public:

    STEdgeLinear(int _n) {
        n = _n;
        node = new Node[n + 1];
        for (int v = 1; v <= n; v++) {
            node[v].parent = 0;
        }
    }

	~STEdgeLinear() {
		delete [] node;
	}

    /// <summary>
    /// Return the root of the tree containing v.
    /// </summary>
    /// <param name="v">Input vertex.</param>
    /// <returns></returns>
    inline int GetRoot(int v) {
		const bool verbose = false;
		if (verbose) {
			int p, count = 0;
			while ((p=node[v].parent) != 0) {v = p; ++count;}
			if (count>0) fprintf (stderr, "%d ", count);
			return v;

		} else {
			int p;
			while ((p=node[v].parent) != 0) {v = p;}
			return v;
		}
    }

    /// <summary>
    /// Return the cost of the arc between v and its parent. Throws 
    /// an exception if v is the root.
    /// </summary>
    /// <param name="v">Input vertex.</param>
    /// <returns></returns>
    inline CostType GetCost(int v) {
        if (v == 0) fatal ("Invalid node id.");
        if (node[v].parent == 0) fatal("root has no parent edge");
        return node[v].cost;
    }

    
    /// <summary>
    /// Return the parent of vertex v (zero if it's a root).
    /// </summary>
    /// <param name="v">Input vertex.</param>
    /// <returns></returns>
    inline int GetParent(int v) const {
        return node[v].parent;
    }


    /// <summary>
    /// Return the id of the bottom vertex of the minimum-cost edge on the path from v to 
    /// the root. Returns 0 if v is the root.
    /// </summary>
    /// <param name="v">Input vertex.</param>
    /// <returns></returns>
    inline int GetMinCost(int v) {
        int w = 0;
        CostType mincost = node[v].cost; //could be undefined, don't care

        int p = node[v].parent;
        while (p != 0) {
            if (node[v].cost <= mincost) {
                w = v;
                mincost = node[v].cost;
            }
            v = p;
            p = node[v].parent;
        }
        return w;
    }

    /// <summary>
    /// Increase by 'cost' the cost of all arcs on the path from v to the root.
    /// </summary>
    /// <param name="v">Starting vertex.</param>
    /// <param name="cost">Amount to increment.</param>
    inline void AddCost(int v, CostType cost) {
        int p = GetParent(v);
        while (p != 0) {
            node[v].cost += cost;
            v = p;
            p = node[v].parent;
        }
    }

    /// <summary>
    /// Insert edge (v,w) with cost c. Does nothing if v is not a root, 
    /// or if v and w are the same component.
    /// </summary>
    /// <param name="v">First endpoint.</param>
    /// <param name="w">Second endpoint.</param>
    /// <param name="cost">Edge cost.</param>
    inline void Link(int v, int w, CostType c) {
        if (node[v].parent != 0) return;
        if (GetRoot(w) == v) return;
        m++; //one more edge!
        node[v].parent = w;
        node[v].cost = c;
    }

    /// <summary>
    /// Cut the edge between v and its parent and returns its cost. 
    /// Throws an exception if v is a root.
    /// </summary>
    /// <param name="v">Input vertex.</param>
    /// <returns>Cost of parent edge.</returns>
    inline CostType cut(int v) {
        if (node[v].parent == 0) fatal ("root cannot be cut");
        m--;
        node[v].parent = 0;
        return node[v].cost; //no need to actually change the cost
    }

    /// <summary>
    /// Reverse the direction of every edge on the path from v to w
    /// </summary>
    /// <param name="v">Starting vertex.</param>
    inline void Evert(int v) {
        int s = v; //source

		const bool verbose = false;

        CostType curcost = node[v].cost;
        int curparent = node[v].parent;

		int ecount = 0;

        while (curparent!=0) {
			Node &curnode = node[curparent];

			/*
            CostType nextcost = node[curparent].cost;
            int nextparent = node[curparent].parent;
            node[curparent].SetParentCost(v,curcost);
			*/
			CostType nextcost = curnode.cost;
			int nextparent = curnode.parent;
			curnode.SetParentCost(v,curcost);

			if (verbose) ecount ++;

            v = curparent;
            curcost = nextcost;
            curparent = nextparent;
        }
		if (verbose) fprintf (stderr, "%d ", ecount);
        node[s].parent = 0;
    }

    /// <summary>
    /// Get number of vertices in the forest.
    /// </summary>
    /// <returns>Number of vertices in the forest.</returns>
    int GetVertexCount() { return n; }

    
    /// <summary>
    /// Find number of edges in the forest.
    /// </summary>
    /// <returns>Number of edges in the tree.</returns>
    int GetEdgeCount() { return m; }
};
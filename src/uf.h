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

class UnionFind {
private:
//#define USE_UF_STRUCT 
#ifdef USE_UF_STRUCT	
	struct UFUnit {
        int p;
        int size;
		inline void SetPS(int _p, int _s) {p = _p; size = _s;}
    };
	inline int GetSize(int i) const {return uf[i].size;}
	inline int GetParent(int i) const {return uf[i].p;}
	inline void SetSize(int i, int s) {uf[i].size = s;}
	inline void SetParent(int i, int _p) {uf[i].p = _p;}
	inline void IncSize(int i, int s) {uf[i].size += s;}
	inline void ResetParentSize(int i) {uf[i].SetPS(i,1);}
	inline void AllocateStructures() {uf = new UFUnit[n+1];}
	inline void FreeStructures() {delete [] uf;}
    UFUnit *uf; 
#else
	int *ufmemory;
	int *parent;
	int *size;
	inline int GetSize(int i) const {return size[i];}
	inline int GetParent(int i) const {return parent[i];}
	inline void SetSize(int i, int s) {size[i] = s;}
	inline void SetParent(int i, int p) {parent[i] = p;}
	inline void IncSize(int i, int s) {size[i] += s;}
	inline void ResetParentSize(int i) {parent[i]=i; size[i]=1;;}
	inline void AllocateStructures() {
		//ufmemory = new int [n + n + 2];
		//parent = ufmemory;
		//size = &ufmemory[n+1];
		parent = new int [n+1];
		size = new int [n+1];
	}
	inline void FreeStructures() {
		//delete [] ufmemory;
		delete [] size;
		delete [] parent;
	}

#endif



	int n; //data structure supports numbers from 0 to maxid


public:
	// total capacity of the data structure
    inline int Elements() const {return n;}

	// number of elements in a set
	inline int ComponentSize(int i) {return GetSize(Find(i));}
	//{return uf[Find(i)].size;}

	// number of independent groups 
    //inline int GroupCount(int i) const {return groups;}
    	
	void Output() {
        for (int i=1; i<=n; i++) {
            int r = Find(i);
            fprintf(stderr, "%d:%d (s%d)", i, r, GetSize(r));
        }
		fprintf (stderr, "\n");
    }

    /// <summary>
    /// Constructor.
    /// </summary>
    /// <param name="_n">Number of elements (=maximum id)</param>
    UnionFind(int maxid) {
		n = maxid;
		AllocateStructures();
		Reset();
	} 
	
	~UnionFind() {FreeStructures();}

    /// Partition elements into singletons (in linear time).
    void Reset() {
		for (int i=0; i<=n; i++) ResetParentSize(i);
		
		//{uf[i].SetPS(i,1);}
    }

    /// <summary>
    /// Function that resets a single element; BE CAREFUL WHEN USING IT!
    /// </summary>
    /// <param name="i">the group to be reset</param>
    inline void Reset(int i) {
		// maybe ResetParentSize is enough

		if (GetParent(i)!=i) {
			SetParent(i,i);
			//groups --;
		}
		SetSize(i,1);
		
		/*
		if (uf[i].p != i) {
            uf[i].p = i;
            groups--;
        }
        uf[i].rank = 0;
        uf[i].size = 1;*/
    }
      
    
    // Find representative of the set containing x.
	// x: an element of the set
	// output: set representative
    inline int Find (int x) {
		//fprintf (stderr, "Finding %d (this=%p,parent=%p,size=%p).\n", x, this, parent, size);
		//fflush (stderr);
		//optimizing the common case: parent is already the root
		int p = GetParent(x); //parent
		int r = GetParent(p); //grandparent
		if (r==p) return r; //common case!

		//keep looking for the root
		while ((p=GetParent(r))!=r) {r = p;}

		//compress path
		while ((p=GetParent(x))!=r) {
			SetParent(x,r);
			x = p;
		}

        return r;
    }

	// get the size of set containing x
    inline int Size(int x) {
		return GetSize(Find(x));
    }

    // Join sets containing x and y (if needed), returns true iff successful 
    inline bool Union (int x, int y) {
		// pick roots of both sets
		x = Find(x); y = Find(y);
    	if (x == y) return false; //same set: nothing to do

		// heaviest element becomes new root of the other
		int sx = GetSize(x);
        int sy = GetSize(y);

        //x root
        if (sx > sy) { // x becomes root
			IncSize(x,sy);
			SetParent(y,x);
        } else { // y becomes root
			IncSize(y,sx);
			SetParent(x,y);
		}
		return true;


		/*
		// heaviest element becomes new root of the other
		int sx = uf[x].size;
        int sy = uf[y].size;

        //x root
        if (sx > sy) { // x becomes root
            uf[x].size += sy;
            uf[y].p = x;
        } else { // y becomes root
			uf[y].size += sx;
			uf[x].p = y;
		}
		return true;
		*/
    }
};











/*

class UnionFindRank {
private:
    struct UFUnit {
        int p;
        int rank;
        int size;

		inline void SetPRS(int _p, int _r, int _s) {
			p = _p;
			rank = _r;
			size = _s;
		}
    };

    UFUnit *uf; 
	int n;      //number of vertices
	int groups; //number of groups

public:
	// total capacity of the data structure
    inline int Elements() const {return n;}

	// number of elements in a set
	inline int ComponentSize(int i) {return uf[Find(i)].size;}

	// number of independent groups 
    inline int GroupCount(int i) const {return groups;}
    	
	void Output() {
        for (int i=1; i<=n; i++) {
            int r = Find(i);
            fprintf(stderr, "%d:%d (r%d,s%d)", i, r, uf[r].rank, uf[r].size);
        }
		fprintf (stderr, "\n");
    }

    /// <summary>
    /// Constructor.
    /// </summary>
    /// <param name="_n">Number of elements (=maximum id)</param>
    UnionFindRank(int _n) {
       n = _n;
       uf = new UFUnit[n+1];
       Reset();
   } 

   ~UnionFindRank() {delete [] uf;}

    /// Partition elements into singletons (in linear time).
    void Reset() {
        for (int i=1; i<=n; i++) {
            uf[i].p = i;
            uf[i].rank = 0;
            uf[i].size = 1;
        }
        groups = n;
    }

    /// <summary>
    /// Function that resets a single element; BE CAREFUL WHEN USING IT!
    /// </summary>
    /// <param name="i">the group to be reset</param>
    //public void Reset(int i) {
    //    if (uf[i].p != i) {
    //        uf[i].p = i;
    //        groups--;
    //    }
    //    uf[i].rank = 0;
    //    uf[i].size = 1;
    //}
      
    
    // Find representative of the set containing x.
	// x: an element of the set
	// output: set representative
    int Find (int x) {
        int p, r = x; 
        while ((p=uf[r].p) != r) {r = p;}
        
        //compress path
        while ((p=uf[x].p) != r) {
            uf[x].p = r;
            x = p;
        }
        return r;
    }

	// get the size of set containing x
    inline int Size(int x) {
        return uf[Find(x)].size;
    }

    // Join sets containing x and y (if needed), returns true iff successful 
    inline bool Union (int x, int y) {
		// pick roots of both sets
		x = Find(x); 
        y = Find(y);
    	if (x == y) return false; //same set: nothing to do

		// pick highest rank
		groups --;
        int rx = uf[x].rank;
        int ry = uf[y].rank;

        //x root
        if (rx > ry) { // x root
            uf[x].size += uf[y].size;
            uf[y].p = x; //x is the new root
        } else { // y root
	        if (rx == ry) {uf[y].rank++;} //increase rank if needed
			uf[y].size += uf[x].size;
			uf[x].p = y;
		}
		return true;
    }

	~UnionFindRank() {delete [] uf;}
};
*/

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

template <class T> class PairingNode {
	public: 
        int child;
        int sibling;
        int parent; //0:isolated >0:actual parent <0:has root '-parent'
        T nodevalue;

		void SetSiblingParentValue(int s, int p, T v) {
            sibling = s;
            parent = p;
            nodevalue = v;
        }
        
        //node[x].Sibling = nil; //no siblings
        //node[x].Value = value; //assign value to x
        //node[x].Parent = nil;

		/*
        public int Sibling {
            get { return sibling; }
            set { sibling = value; }
        }*/

        bool IsBetter (PairingNode x) {
            //return (nodevalue.CompareTo(x.Value) < 0);
			return (nodevalue < x.nodevalue);
        }

        //return (node[x].Value.CompareTo(node[y].Value) < 0);

		/*
        inline int Child {
            get { return child; }
            set { child = value; }
        }

        inline T Value {
            get { return nodevalue; }
            set { nodevalue = value; }
        }

        inline int Parent {
            get { return parent; }
            set { parent = value; }
        }*/

        void Reset(int nil, int p) {
            child = sibling = nil;
            parent = p;
        }
    };


template <class T> class PairingHeap {
public:

private:
	void fatal (const string &msg) {
		fprintf (stderr, "ERROR::PairingHeap: %s.\n", msg.c_str());
		fflush(stderr);
		exit(-1);
	}

	bool verbose;
    int emax; //maximum id used in the heap (from 1 to n)
    int hmax; //maximum heap name (1 to hmax)
    vector <PairingNode<T> > node; //* node;
    //int *name2root; //name2root[h]: root of heap named h. May be nil.
	vector<int> name2root;
    int nil; //null pointer (should perhaps replace this by zero)

public:
    //get parent node
    inline int GetParent(int i) const {return node[i].parent;}
    inline int GetSibling(int i) const {return node[i].sibling;}
    inline T GetValue(int i) const {return node[i].nodevalue; }
    inline int GetChild(int i) const {return node[i].child;}

    void Reset() {
        for (int h = 1; h <= hmax; h++) name2root[h] = nil;
        for (int e = 1; e <= emax; e++) node[e].Reset(nil,nil);
		name2root[0] = -1;
		node[0].Reset(nil,nil);
    }

    /// <summary>
    /// Create collection of _hmax heaps with _emax elements in total.
    /// </summary>
    /// <param name="_emax">Number of elements.</param>
    /// <param name="_hmax"></param>
public:
    PairingHeap(int _emax, int _hmax) {
        verbose = false;
        emax = _emax; //number of elements to represent (maximum id) 
        hmax = _hmax;
        nil = -1;
        //return;
        //node = new PairingNode<T>[emax + 1]; //each element has a node
        //name2root = new int[hmax + 1]; //pointer from each heap name to its root node
		node.resize(emax+1);
		name2root.resize(hmax+1);
        //stack = new int [hmax];
        Reset();
        if (verbose) fprintf (stderr, "Created meldable structure with %d elements and %d heaps.\n", emax, hmax);
    }

	~PairingHeap() {
		//fprintf (stderr, "Deleting pairing heap... ");
		//fflush(stderr);
		//delete [] name2root;
		//fprintf (stderr, "*"); fflush(stderr);
		//delete [] node;
		//fprintf (stderr, "done deleting pairing heap.\n");
		//fflush(stderr);
	}

    /// <summary>
    /// Create an empty heap.
    /// </summary>
    /// <param name="h">Heap name.</param>
    /*
    public void MakeHeap(int h) {
        if (name2root[h] != nil) {
            Console.Error.WriteLine("WARNING: CREATING PREVIOUSLY EXISTING HEAP.");
        }
        name2root[h] = nil;
    }*/


    /// <summary>
    /// Return the top element of a heap, without changing it.
    /// </summary>
    /// <param name="h">Heap name.</param>
    inline int FindMin(int h, T &value) {
        int r = name2root[h];
        if (r == nil) fatal ("Cannot take minimum from empty heap");
        value = GetValue(r);
        return r;
    }

    inline bool Isolated(int x) const {
        return (GetParent(x) == nil);

    }

private:
    inline bool IsBetter(int x, int y) {
        return node[x].IsBetter(node[y]); //SHOULD MAKE THIS SIMPLER!
        //return (node[x].Value.CompareTo(node[y].Value) < 0);
        //return (GetValue(x).CompareTo(GetValue(y)) < 0);
    }


    /// <summary>
    /// Link heaps rooted at r1 and r2. Return the root of the resulting heap.
    /// </summary>
    /// <param name="h1">One heap, and the name of the resulting heap.</param>
    /// <param name="h2">Second heap</param>
    int Link(int r1, int r2) {
        //simple cases: at least one is empty
        if (r1 == nil) return r2;
        if (r2 == nil) return r1;

        int p, c; //parent, child
        if (IsBetter(r2, r1)) {
        //if (node[r2].IsBetter(node[r1])) {
            p = r2; c = r1;
        } else {
            p = r1; c = r2;
        }

		if (p<1 || p>emax) fatal ("parent out of range");
		if (c<1 || c>emax) {
			fprintf (stderr, "p=%d c=%d\n", p, c);	
			fatal ("child out of range");
		}

        int oldchild = node[p].child;

        node[c].sibling = oldchild;
        if (oldchild != nil) node[oldchild].parent = c;
        node[p].child = c;
        node[c].parent = p;
        return p;
    }

    /// <summary>
    /// Assign name h to the heap rooted at r.
    /// </summary>
    void SetName(int h, int r) {
        if (r != nil) node[r].parent = -h; //make root point to heap
        name2root[h] = r; //keep pointer to root
    }

    void FreeName(int h) {
        int r = name2root[h];
        if (r != nil) {
            name2root[h] = nil;
            node[r].parent = nil;
        } //else throw new Exception ("FreeName cannot be called on empty heap.");
    }

	public:
    void Insert(int h, int x, T value) {
        //if (!Isolated(x)) {}; //throw new Exception("Only isolated elements can be inserted");
		if (h<1 || h>hmax) fatal ("heap out of range");
		if (x<1 || x>emax) fatal ("element out of range");

        if (node[x].parent != nil) {
			fprintf (stderr, "Node %d has parent %d.\n", x, node[x].parent);	
			fatal ("Only isolated elements can be inserted");
		}

        node[x].sibling = nil; //no siblings
        node[x].nodevalue = value; //assign value to x
        node[x].parent = nil; //THIS IS PROBABLY NOT NECESSARY
        //node[x].SetSiblingParentValue(nil, nil, value);
        int r = name2root[h];

        r = Link(r, x);
        SetName(h, r);
        //if (r != nil) node[r].Parent = -h; //make root point to heap
        //name2root[h] = r; //keep pointer to root
    }


public:
    /// <summary>
    /// Merge two heaps, h1 and h2, creating a new heap h1
    /// </summary>
    void Merge(int h1, int h2) {
        if (h1 == h2) return;
        //Console.Error.WriteLine("Merging {0}:{1}.", h1, h2);
        int r1 = name2root[h1];
        int r2 = name2root[h2];
        //Console.Error.WriteLine("Done with names.");
        FreeName(h1);
        FreeName(h2);

        int r = Link(r1, r2);
        //if (r != nil) {
        SetName(h1, r); //do this even if r is null...
        //else throw new Exception ("The result of a merge cannot be null.");
    }

private:
    /// <summary>
    /// Create a single heap from the sibling list starting at v.
    /// </summary>
    /// <param name="v"></param>
    /// <returns></returns>
    int PairUp (int current) {
        int firstdone = nil;
        int next = 0;

        //Console.Error.WriteLine("Pairing up, starting at {0}.", current);

        //state: (<<<firstdone<<<) (>>>current next nextnext>>>) 
        //(the firstdone list is reversed)

        //FIRST PASS: pair elements from left to right while 
        //reversing the newly created list.
        while (true) {
            next = node[current].sibling;
            int nextnext = (next!=nil) ? node[next].sibling : nil;
            int temp = Link(current, next);
            node[temp].sibling = firstdone;
            firstdone = temp;

            if (nextnext == nil) break;
            current = nextnext; //advance
        }

        //Second pass: just link everybody
        current = firstdone;
        next = node[current].sibling;
        while (next!=nil) {
            int nextnext = node[next].sibling;
            current = Link(current, next); //join the first two trees
            node[current].sibling = nextnext; //is this really necessary?
            next = nextnext; //advance pointers
        }

        return current;
    }

public:
	void DecreaseKey(int h, int x, T value) {
        if (Isolated(x)) fatal("Cannot decrease key of isolated node.");
        fprintf (stderr, "Warning: not checking for improvement.");
        //if (node[x].Value < value) throw new Exception("DecreaseKey cannot increase key.");
        node[x].Value = value;

        int p = node[x].parent;
        if (p == nil) return;

        //remove subtree rooted at x
        int s = node[x].sibling;

        //there is a parent in the binary tree
        if (node[p].sibling == x) { //actually a sibling
            node[p].sibling = s;
        } else { //really a child
            node[p].child = s;
        }

        if (s != nil) {
            node[s].parent = p;
        }

        node[x].sibling = nil;
        node[x].parent = nil;

        //combine the subtree with the current root
        int r = name2root[h];
        r = Link(r, x);
        SetName(h, r);
    }

public:
    void Output(int h) {
        fprintf (stderr, "Outputting heap %d.\n", h);
        int r = name2root[h];
        fprintf (stderr, "Root is %d.\n", r);

        for (int e = 1; e <= emax; e++) {
            fprintf (stderr, "%d: p%d s%d c%d", e, node[e].parent, node[e].sibling, node[e].child);
        }
    }

	/*
    public delegate T Transform(T x);


    /// <summary>
    /// Traverse the entire heap h, applying a transformation to each element.
    /// </summary>
    /// <param name="h">The heap.</param>
    /// <param name="transf">Function from T to T.</param>
    public void Traverse(int h, Transform transf) {
        int v = name2root[h];      
        int from = 0; //0: parent; 1:child; 2:sibling
        int travcount = 0;

        bool verbose = false;
        Console.Error.WriteLine("Traversing heap {0}, root {1}... ({2} {3})", h, v, node[v].Parent, nil);

        while (v >= 0) {
            int p = node[v].Parent;
            int c = node[v].Child;
            int s = node[v].Sibling;

            if (from == 0) {
                if (verbose) {
                    Console.Error.WriteLine("v={0} p={1} c={2} s={3} ", v, p, c, s);
                    Console.Error.Write("BEFORE{0}:{1} -> ", v, node[v].Value);
                }
                node[v].Value = transf(node[v].Value);
                if (verbose) Console.Error.Write(" -> AFTER{0}:{1}\n", v, node[v].Value);
                
                
                travcount++;
            }

            //go to child if possible
            if ((from == 0) && (c != nil)) {
                v = c; //from remains zero
            } else { //after child, go to sibling
                if ((from != 2) && (s != nil)) {
                    v = s;
                    from = 0;
                } else { //after child and sibling, go back up
                    //going to parent (from doesn't matter if v is the root)
                    from = 0;
                    if (p > 0) from = (node[p].Child == v) ? 1 : 2;
                    v = p;
                }
            }
        }
        Console.Error.WriteLine(" ({0} elements traversed).", travcount);

    }*/

public:
    int DeleteMin(int h, T &value) {
        int r = name2root[h];
        if (r == nil) fatal("Cannot delete from an empty heap.");

        value = node[r].nodevalue;
        int c = node[r].child;
        node[r].parent = node[r].child = nil; //make node r isolated

        if (c != nil) c = PairUp(c);
        SetName(h,c);

        return r;
    }


    inline bool IsEmpty(int h) const {
		/*
		fprintf (stderr, "Checking if %d is empty.\n", h);
		fflush (stderr);
		fprintf (stderr, "Checking if %d is empty (%d).\n", h, name2root[h]);
		fflush (stderr);*/

        return (name2root[h] == nil);
    }
};

/*
public class TraversablePairingHeap<T> : PairingHeap<T> where T : IComparable<T> {

    public TraversablePairingHeap (int _emax, int _hmax) : base (_emax, _hmax) {
    
    }


}*/

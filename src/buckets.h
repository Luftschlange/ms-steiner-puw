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
/// StaticBuckets: maps elements (0...maxid) to buckets (0...nbuckets)
/// Once in a bucket, an element cannot be moved or deleted.
/// And it's impossible to know which bucket an element belongs to.
/// </summary>
class StaticBuckets {
private:
    int maxid; //universe is 0...maxid
    int maxbucket; //buckets 0...maxbucket

    int *next;  //next[v]: successor of v in its list (-2:none, -1:end of list) 
    int *first; //first[b]: first element in b's bucket (-1: none)

	void fatal (const string &msg) {
		fprintf (stderr, "StaticBuckets::Error: %s.\n", msg.c_str());
		exit(-1);
	}

public:
    void Reset() {
        int i;
        for (i = 0; i <= maxid; i++) next[i] = -2; //element is no list
        for (i = 0; i <= maxbucket; i++) first[i] = -1; //all buckets empty
    }

    void CheckEmpty() {
        fprintf (stderr, "Checking if StaticBuckets is empty... ");
        int i;
        for (i=0; i<=maxid; i++) if (next[i] != -2) fatal ("StaticBuckets.next has wrong value");
        for (i=0; i<=maxbucket; i++) if (first[i] != -1) fatal ("StaticBuckets.first has wrong value");
        fprintf (stderr, "done.\n");
    }

    inline bool IsEmpty(int b) { return (first[b] == -1); }

    StaticBuckets (int _maxid, int _maxbucket) {
        maxid = _maxid;
        maxbucket = _maxbucket;

        next = new int[maxid + 1];
        first = new int[maxbucket + 1];

        Reset();
    }

	~StaticBuckets() {
		delete [] next;
		delete [] first;
	}

    void Reset(int b) {
        int i = first[b];
        if (i > -1) {
            do {
                int j = next[i];
                next[i] = -2;
                i = j;
            } while (i > 0);
        }
        first[b] = -1;
    }

    /// <summary>
    /// Insert an element into a bucket. Element must not be in any bucket already.
    /// </summary>
    /// <param name="e">the element</param>
    /// <param name="b">the bucket</param>
    void Insert(int e, int b) {
        if (next[e] >= -1) fatal ("Elements cannot be reassigned.");
        next[e] = first[b];
        first[b] = e;
    }

    /// <summary>
    /// Enumerate all elements in bucket b.
    /// </summary>
    /// <param name="b">the bucket</param>
    /// <returns></returns>
    
	/*
	void GetBounds (int b, int &p, int &pend) {

	}*/

	inline int GetFirst(int b) {return first[b];}
	inline int GetNext(int e) {return next[e];}
	
	/*
	public IEnumerable<int> ElementEnumerator(int b) {
        for (int e = first[b]; e >= 0; e = next[e]) { yield return e; }
    }*/
};

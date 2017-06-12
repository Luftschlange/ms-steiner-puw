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

/// Trivial queue used for BFS-like applications. May store elements from
/// 0 to maxid, and supports no more than maxid insertions (it does not 
/// wrap-up).

#include <cstdio>
#include <cstdlib>

class SimpleQueue {
private:
	int maxid;
    int nextins; //position of next insertion
    int nextrem; //position of next removal
    int *queue;

	void fatal (const string &msg) {
		fprintf (stderr, "ERROR::SimpleQueue: %s.\n", msg.c_str());
		exit(-1);
	}

public:
    /// Create queue with of a given max size.
    SimpleQueue(int _maxid) {
        maxid = _maxid;
        queue = new int [maxid+1];
        Reset();
    }

	~SimpleQueue() {delete [] queue;}

    inline void Reset() {
        nextins = 0;
        nextrem = 0;
    }

    inline int Count() const {
        return (nextins - nextrem);
    }

    inline void Insert (int x) {
        if (nextins > maxid) fatal ("insertion limit exceeded");
        queue[nextins++] = x;
    }

    inline bool IsEmpty() const {
        return (nextrem >= nextins);
    }

    inline int Remove() {
        if (IsEmpty()) fatal ("cannot remove from empty queue");
        return queue[nextrem++];
    }

    inline void Output() {
		fprintf (stderr, "Queue: [");
        for (int i=nextrem; i<nextins; i++) {
            fprintf(stderr, " %d", queue[i]);
        }
        fprintf(stderr, " ]");
    }
};
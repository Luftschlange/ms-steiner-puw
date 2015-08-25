/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

//#ifndef UNIVERSESET_H
//#define UNIVERSESET_H
//#include "rfw_random.h"

#pragma once

#include <cstdio>
#include <cstdlib>


class UniverseSet {
	private:
		int *perm;
		int *pos; //pos[v]: position of element v
		int maxid;  //maximum allowed identifier
		int nextpos; //position of next insertion
		int nil;

    
		/// <summary>
		/// Swap element in position pi into position nextpos.
		/// </summary>
		/// <param name="pi">current element position</param>
		void MakeNext(int pi) {
			int i = perm[pi]; //element originally in position i
			int j = perm[nextpos]; //element originally in nextpos

			//i goes to nextpos
			perm[nextpos] = i;
			pos[i] = nextpos;

			//j goes to pi
			perm[pi] = j;
			pos[j] = pi;
		}

		void fatal (const string &msg) {
			fprintf (stderr, "%s\n", msg.c_str());
			exit(-1);
		}

		void die (const string &msg) {
			fprintf(stderr, "(UniverseSet): %s.\n", msg.c_str());
			fatal(msg);
		}


	public:
		/// <summary>
		/// Create a set representing elements 0...n.
		/// </summary>
		/// <param name="n"></param>
		UniverseSet (int _maxid) {
			maxid = _maxid;
			nil = maxid + 1;
			perm = new int [_maxid+1];
			pos = new int [_maxid+1];
			HardReset();
		}


		/*
		void Permute(OptRandom random) {
			for (int i=nextpos-1; i>0; i--) {
				int j = random.GetInteger(0, i);

				int t = perm[i];
				perm[i] = perm[j];
				perm[j] = t;

				pos[perm[i]] = i;
				pos[perm[j]] = j;
			}
		}*/

		inline int PickAny() {
			if (IsEmpty()) die ("Cannot pick element from an empty set.");
			return (perm[0]);
		}

		/// <summary>
		/// Empty the set.
		/// </summary>
		void HardReset()
		{
			//Console.WriteLine("Resetting...\n");
			for (int i=0; i<=maxid; i++) {
				pos[i] = i;
				perm[i] = i;
			}
			nextpos = 0;
		}

		void Reset() {
			while (!IsEmpty()) {Remove(perm[0]);}
		}


		inline bool Contains(int i) const {return (pos[i]<nextpos);}
		inline bool IsEmpty() {return (nextpos==0);}

		inline bool Insert(int i) {
			int pi = pos[i];
			if (pi < nextpos) return false; //already there: nothing to do
			MakeNext(pi);
			nextpos++;
			return true;
		}

		inline bool Remove (int i) {
			int pi = pos[i];
			if (pi >= nextpos) return false;
			nextpos--;
			MakeNext(pi);
			return true;
		}


	/*
	int PickRandom(RFWLocalRandom &r) {
		if (IsEmpty()) die("Cannot pick element from an empty set.");
		return (perm[r.getInteger(0,nextpos-1)]);
	}
		
    int PickRandom() {
        if (IsEmpty()) die("Cannot pick element from an empty set.");
		return (perm[RFWRandom::getInteger(0,nextpos-1)]);
    }*/

    inline int Count() {return (nextpos);}

    inline int MaxId() {return maxid;}

    inline void Copy(UniverseSet *s) {
        if (s->MaxId() != MaxId()) die("Cannot copy from incompatible set");

        Reset(); //make this set empty
		Unite(s); //take the union with the other set

		/*
		int i, end;
		for (s->GetBounds(i, end); i<end; i++) {
			int v = s->PickPos(i);

        //foreach (int v in s.ElementEnumerator()) {
            Insert(v);
        }*/
    }

	inline void Unite (UniverseSet *s) {
		if (s->MaxId() > MaxId()) die ("Cannot unite potentially larger set.\n");
		int i, end;
		for (s->GetBounds(i,end); i<end; i++) {
			Insert(s->PickPos(i));
		}

	}

	inline void GetBounds(int &start, int &end) {
		start = 0;
		end = nextpos;
	}

	inline int PickPos(int p) {
		if (p<0 || p>nextpos) die ("picking from invalid position");
		return perm[p];
	}



    void Output() {
        int i, end;
        for (i=0; i<=maxid; i++) {
            if (Contains(i)) fprintf (stderr, "#");
            else fprintf (stderr, "_");
        }
		for (GetBounds(i,end); i<end; i++) {
			fprintf (stderr, " %d", PickPos(i));
		}
		fprintf (stderr, "\n");

        
        CheckConsistency();
    }

    void CheckConsistency() {
        fprintf (stderr, "Check consistency.");

        for (int i=0; i<=maxid; i++) {
            if (perm[pos[i]] != i) fprintf(stderr, "Inconsistency at position %d.\n", i);
            if (pos[perm[i]] != i) fprintf(stderr, "Inconsistency for element %d.\n", i);
        }
    }

	~UniverseSet() {
		delete [] perm;
		delete [] pos;
	}
};


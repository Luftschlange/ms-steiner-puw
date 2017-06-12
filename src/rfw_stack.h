/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*-------------------------------------------------------
 | RFWStack: stack of elements of type T. Upper bound
 |           n on maximum size must be known in advance.
 |           Initialization done in O(n) time, all other
 |           operations in constant time.
 |
 | author: Microsoft Corporation
 *------------------------------------------------------*/

#ifndef RFW_STACK_H
#define RFW_STACK_H
#include <cassert>

template <class T, bool debug=false> class RFWStack {
	private:
		T *stack;
		int top, size;

	public:
		inline bool isFull() {return (top==size);}

		inline void push (T i) {
			if (debug) assert (top<size);
			stack[++top] = i;
		}

		inline T peekTop() {
			return stack[top];
		}

		inline T pop() {
			if (debug) assert (!isEmpty());
			return (stack[top--]);
		};

		inline bool isEmpty() {return (top==0);};

		inline int getNElements() {return top;}

		// peek element at position i; first allowed position is 1
		inline T peek (int p) {
			if (debug) assert (p>0 && p<=top);
			return (stack[p]);
		}

		RFWStack (int s) {
			if (debug) assert (s>0);
			top = 0;
   			size = s;
			stack = new T [size+1];
		};

		inline void reset() {top = 0;}

		~RFWStack () {delete [] stack;};
};

#endif

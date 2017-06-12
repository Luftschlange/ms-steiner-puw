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
using namespace std;

/// <summary>
/// Binary heap holding integral keys between 1 and maxid, with
/// values of type T. Lower values have higher priority.
/// </summary>
template <class T> class BinaryHeap {

	private:
		/// <summary>
		/// Basic heap element structure supporting generic types
		/// </summary
		typedef struct {
			unsigned label; //key of the element being stored
			T value; //its value
		} HeapUnit;

		// Heap elements
		int lastpos;   //last position in which insertion is allowed
		HeapUnit *h;  //heap elements (elements will be stored in positions h[1] to h[lastpos])
		int nextpos;   //position of next element to be inserted
		int *heappos; //position of each element in the heap (0 is not yet inserted)

		/// <summary>
		/// Swaps elements in positions a and b.
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		inline void Swap(int a, int b) {
			HeapUnit tmp = h[a];
			h[a] = h[b];
			h[b] = tmp;
			heappos[h[a].label] = a;
			heappos[h[b].label] = b;
		}

		/// <summary>
		/// Standard 'up' operation performed on the element in the n-th position
		/// </summary>
		/// <param name="n">Position.</param>
		inline void Up(int n) {
			//while ((n > 1) && (GetValue(n / 2).CompareTo(GetValue(n)) > 0)) {
			while ((n>1) && IsBetter(GetValue(n), GetValue(n/2))) {
				Swap(n, n/2);
				n = n/2;
			}
		}

		inline T GetValue(int pos) const { return h[pos].value; }
		inline int GetLabel(int pos) const { return (h[pos].label); }

		/// <summary>
		/// Move the element currently in position n down to an appropriate position.
		/// </summary>
		/// <param name="n">Position.</param>
		void Down(int n) {
			HeapUnit tmp;
			T minchild;
			int two_n;

			if ((two_n = 2*n) < nextpos) {
				tmp = h[n]; //Save the original value
				do {
					bool first = true;
					minchild = GetValue(two_n); //guess first child is the best
					//if (((two_n + 1) < nextpos) && (minchild.CompareTo(GetValue(two_n + 1)) > 0)) {
					if ((two_n+1 < nextpos) && IsBetter(GetValue(two_n+1), minchild)) { //unless the second is better
						minchild = GetValue(two_n + 1);
						first = false;
					}
					//if (tmp.value.CompareTo(minchild) > 0) { //tmp is higher than the child
					if (IsBetter(minchild, tmp.value)) { //is the child better than the current value?
						//if (minchild.CompareTo(GetValue(two_n)) == 0) {
						if (first) { //minchild was the first
							h[n] = h[two_n];
							heappos[h[n].label] = n;
							n = two_n;
						} else { //minchild is the second
							h[n] = h[two_n + 1];
							heappos[h[n].label] = n;
							n = two_n + 1;
						}
						two_n = 2 * n;
					}
					else {
						break;
					}
				} while (two_n < nextpos);
				h[n] = tmp;
				heappos[h[n].label] = n;
			}
		}


		/// <summary>
		/// Find out if 'a' is strictly better than 'b'
		/// </summary>
		/// <param name="a">first value</param>
		/// <param name="b">second value</param>
		/// <returns>true iff a is strictly better (smaller) than b</returns>
		inline bool IsBetter(T a, T b) const {
			//return (a.CompareTo(b) < 0);
			return (a<b);
		}
	public:
		/// <summary>
		/// Creates a new heap that can hold elements from 1 to n.
		/// </summary>
		/// <param name="n">Maximum id (minimum is 1).</param>
		BinaryHeap(int n) {
			nextpos = 1;
			lastpos = n;

			h = new HeapUnit[n+1];
			heappos = new int[n+1];

			for (int v = 0; v <= n; v++) heappos[v] = 0;
		}

		/*
		* Destructor
		*/
		virtual ~BinaryHeap() {
			delete[] h;
			delete[] heappos;
		};


		/// <summary>
		/// Remove the element with the highest priority from the heap. Its label and value are returned.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		void RemoveFirst(unsigned &label, T &value) {
			label = h[1].label;  
			value = GetValue(1);

			h[1] = h[nextpos-1];	
			heappos[h[1].label] = 1;
			
			heappos[label] = 0;

			nextpos --;					
			Down(1);					
		}

		T PeekTopValue() {return GetValue(1);}

		/*
		void PeekFirst(unsigned &label, T &value) {
			label = h[1].label;  
			value = GetValue(1);
		}*/

		T MinValue(){
			return GetValue(1);
		}

		/// <summary>
		/// Heapify operation
		/// </summary>
		void Heapify() {
			int i;
			for (i = (nextpos - 1) / 2; i >= 1; i--) {
				Down(i);
			}
		}


		/// <summary>
		/// Removes the element with the specified label from the heap
		/// </summary>
		/// <param name="label"></param>
		void RemoveElement(unsigned label) {
			//get original position, check if element is actually in the heap
			int pos = heappos[label];
			if (pos!=0) {

				//heappos[label] = lastpos+1;  // Element will no longer belong to the heap
				heappos[label] = 0;
				nextpos--;  // One fewer element

				// If the heap is not empty and this is not the last element, update the heap
				if ((nextpos != 1) && (nextpos != pos)) { 
					h[pos] = h[nextpos];  // Last element replaces the removed one
					heappos[h[pos].label] = pos;  // Let the element know its new position

					//unless parent is strictly better, go up
					//if ((pos > 1) && (h[pos].value.CompareTo(h[pos / 2].value) <= 0)) {
					if ((pos>1) && !IsBetter(h[pos/2].value, h[pos].value)) {
						Up(pos);
					} else {
						Down(pos);
					}
				}
			}
		}


		/// <summary>
		/// Insert and/or decrease the value of an element in the heap. Returns
		/// true iff successful (false if the new value is worse).
		/// </summary>
		/// <param name="label">Label of the new element.</param>
		/// <param name="value">Value of the new element.</param>
		/// <returns>True iff the element is inserted or updated.</returns>
		bool Insert(unsigned label, T value) {
			int prevpos = heappos[label];
			if (prevpos != 0) {  // Already in the heap
				//if (h[prevpos].value.CompareTo(value) < 0) { //do not allow demotions
				if (IsBetter(h[prevpos].value, value)) { //do not allow demotions
					return false;
				} else { //positive updates
					h[prevpos].value = value;
					Up(prevpos);
					return true;
				}
			}

			//not in the heap: insert it
			h[nextpos].label = label;
			h[nextpos].value = value;
			heappos[label] = nextpos;
			Up(nextpos);
			nextpos++;
			return true;
		}

		/// <summary>
		/// Change the key of element 'label' to value. If the element is
		/// not already in the heap, insert it; otherwise, just update it.
		/// The key will be changed regardless of whether the new key has
		/// higher or lower priority.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		void FixKey (unsigned label, T value) {
			int current_pos = heappos[label];
			if (current_pos == 0) { //element is not there: insert it
				h[nextpos].label = label;
				h[nextpos].value = value;
				heappos[label] = nextpos;
				Up(nextpos);
				nextpos++;
			} else {
				T old_value = h[current_pos].value;
				//int sign = value.CompareTo(current_value);
				//if (sign < 0) { //strictly better

				h[current_pos].value = value;

				//if the old value was better, move element down; otherwise, bubble up as needed
				if (IsBetter(old_value, value)) Down(current_pos);
				else Up(current_pos);
			}
		}


		/// <summary>
		/// Add (label,value) as the last element of the current heap;
		/// MAY VIOLATE HEAP ORDER. The heap will be inconsistent
		/// until 'heapify' is called. This operation takes constant time.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		/// <returns></returns>
		bool PushBack(unsigned label, T value) {
			if (Contains(label)) return false;
			h[nextpos].label = label;
			h[nextpos].value = value;
			heappos[label] = nextpos;
			nextpos++;
			return true;
		}

		inline T GetElementValue(int v) { return (GetValue(heappos[v])); }
		inline bool Contains(int v) { return (heappos[v] != 0); }
		inline bool IsEmpty() { return (nextpos == 1); }
		inline int GetSize() { return (nextpos-1); }
		inline int GetMaxSize() { return (lastpos); }

		/// <summary>
		/// Remove all elements from the heap.
		/// (Time proportional to the number of elements.)
		/// </summary>
		void Reset() {
			bool full = false;
			if (full) {
				//Console.Error.WriteLine("Resetting heap (expensive).");
				cerr << "Resetting heap (expensive).\n";
				nextpos = 1;
				for (int i = 0; i < lastpos; i++) {heappos[i] = 0;}
			} else {
				for (int p = 1; p < nextpos; p++) heappos[h[p].label] = 0;
				nextpos = 1;
			}
		}




		/// <summary>
		/// Sanity checks the heap invariants and internal mappings.
		/// DEBUG FUNCTION ONLY.
		/// </summary>
		/// <returns>true iff the data structure passes all tests</returns>
		bool Check() {
			fprintf (stderr, "Checking heap (%d elements)... ", nextpos);
			for (int i = nextpos - 1; i > 1; i--) {
				//if (h[i].value.CompareTo(h[i / 2].value) < 0) {
				if (IsBetter(h[i].value, h[i/2].value)) {
					//Console.WriteLine("{0} ({1}) is child of {2} ({3}). There are {{4} elements.",
					cerr << h[i].value << " (" << i << ") is child of " << h[i/2].value << " (" << i/2 << "). There are " << nextpos-1 << "elements.\n";
						//h[i].value, i, h[i / 2].value, i / 2, nextpos - 1);
					return false;
				}
			}

			for (int i = 1; i < nextpos; i++) {
				if (heappos[h[i].label] != i) {
					//Console.WriteLine("Element {0} thinks it is in position {1}, but it isn't.", h[i].label, i);
					cerr << "Element " << h[i].label << " things it's in position " << i << ", but it isn't.\n";
					return false;
				}
			}

			for (int i = 1; i <= lastpos; i++) {
				int pos = heappos[i];
				if ((pos > 0) && (pos <= lastpos)) {
					if (pos >= nextpos) {
						//Console.WriteLine("Element {0}'s reported position ({1}) is out of bounds (there are only {2} elements).\n", i, pos, nextpos - 1);
						cerr << "Element " << i << "'s reported position (" << pos << ") is out of bounds (there are only " << nextpos-1 << " elements).\n";
						return false;
					}
					if (h[pos].label != i) {
						//Console.WriteLine("Element {0} is not in the position he thinks he is ({1}); {2} is there.\n", i, pos, h[pos].label);
						cerr << "Element " << i << " is not in the position it things it is (" << pos << "); " << h[pos].label << " is there.\n";
						//i, pos, h[pos].label
						return false;
					}
				}
			}

			fprintf (stderr, "passed.\n");

			return true;
		}
};

//#endif
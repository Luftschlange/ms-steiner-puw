/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "rfw_random.h"

const unsigned long int RFWRandom::maxvalue = (unsigned long int)0xffffffff;// (unsigned long int)(-1);

const unsigned long int RFWRandom::MATRIX_A = 0x9908b0dfUL;   /* constant vector a */
const unsigned long int RFWRandom::UPPER_MASK = 0x80000000UL; /* most significant w-r bits */
const unsigned long int RFWRandom::LOWER_MASK = 0x7fffffffUL; /* least significant r bits */

unsigned long RFWRandom::mt[N]; /* the array for the state vector  */
int RFWRandom::mti=N+1; /* mti==N+1 means mt[N] is not initialized */
			

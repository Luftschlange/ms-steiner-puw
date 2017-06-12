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
/* C++ Wrapper written by Microsoft Corporation on original C code by 
   Nishimura and Matsumoto. Original header follows. */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#ifndef RFW_RANDOM_H
#define RFW_RANDOM_H

#include <cstdlib>
#include <cstdio>

class RFWRandom {
	public:
		static const unsigned long int maxvalue;

		/* Period parameters */  
		enum {N=624, M=397};
		static const unsigned long int MATRIX_A;   /* constant vector a */
		static const unsigned long int UPPER_MASK; /* most significant w-r bits */
		static const unsigned long int LOWER_MASK; /* least significant r bits */
	private:

		static unsigned long mt[N]; /* the array for the state vector  */
		static int mti; /* mti==N+1 means mt[N] is not initialized */

			
		/* initializes mt[N] with a seed */
		static void init_genrand(unsigned long s)
		{
			mt[0]= s & 0xffffffffUL;
			for (mti=1; mti<N; mti++) {
				mt[mti] = 
				(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
				/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
				/* In the previous versions, MSBs of the seed affect   */
				/* only MSBs of the array mt[].                        */
				/* 2002/01/09 modified by Makoto Matsumoto             */
				mt[mti] &= 0xffffffffUL;
				/* for >32 bit machines */
			}
		}


		/* generates a random number on [0,0xffffffff]-interval */
		static unsigned long genrand_int32(void)
		{
			unsigned long y;
			static unsigned long mag01[2]={0x0UL, MATRIX_A};
			/* mag01[x] = x * MATRIX_A  for x=0,1 */

			if (mti >= N) { /* generate N words at one time */
				int kk;

				if (mti == N+1)   /* if init_genrand() has not been called, */
					init_genrand(5489UL); /* a default initial seed is used */

				for (kk=0;kk<N-M;kk++) {
					y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
					mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				for (;kk<N-1;kk++) {
					y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
					mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
				mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

				mti = 0;
			}
		  
			y = mt[mti++];

			/* Tempering */
			y ^= (y >> 11);
			y ^= (y << 7) & 0x9d2c5680UL;
			y ^= (y << 15) & 0xefc60000UL;
			y ^= (y >> 18);

			return y;
		}

	public:
		//constructors
		RFWRandom () {randomize(1);}

		//randomize procedures
		static void randomize (unsigned long s) {
			if (s==0) s = 1;
			init_genrand(s);
		}
		
		static unsigned long getRand() {return genrand_int32();}

		//pick an integer uniformly at random between inf and sup (both inclusive)
		static int getInteger (int inf, int sup) {
			if (sup<=inf) return inf;
			unsigned long range, minallowed, u;

			range = (unsigned long)(sup-inf+1);    //number of values allowed
			minallowed = (maxvalue % range) + 1;   //restrict search space to avoid small numbers
			if (minallowed==range) minallowed = 0;
			do {u = getRand();}     //repeat until a good number is found
			while (u < minallowed);

			return (inf + (int)(u % range)); //return a number in the range
		}



		static float getFloat () {return (float)getDouble();} //get a float number in [0;1]
		static double getDouble() {return getDoubleClosed();} //double in the range [0;1]
		static double getDoubleClosed() {return ((double)getRand()/(double)maxvalue);}  //double in the range [0;1]
		static double getDoubleOpen() {return ((double)getRand()/((double)(maxvalue)+1.0));} //double in the range [0;1)
		static bool getBool () {return (getRand() & 1);}

};







class RFWLocalRandom {
	private:
		//static const unsigned long int maxvalue;

		/* Period parameters */  
		//enum {N=624, M=397};
		//static const unsigned long int MATRIX_A;   /* constant vector a */
		//static const unsigned long int UPPER_MASK; /* most significant w-r bits */
		//static const unsigned long int LOWER_MASK; /* least significant r bits */

		unsigned long mt[RFWRandom::N]; /* the array for the state vector  */
		int mti; /* mti==N+1 means mt[N] is not initialized */

			
		/* initializes mt[N] with a seed */
		void init_genrand(unsigned long s)
		{
			mt[0]= s & 0xffffffffUL;
			for (mti=1; mti<RFWRandom::N; mti++) {
				mt[mti] = 
				(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
				/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
				/* In the previous versions, MSBs of the seed affect   */
				/* only MSBs of the array mt[].                        */
				/* 2002/01/09 modified by Makoto Matsumoto             */
				mt[mti] &= 0xffffffffUL;
				/* for >32 bit machines */
			}
		}


		/* generates a random number on [0,0xffffffff]-interval */
		unsigned long genrand_int32(void)
		{
			unsigned long y;
			unsigned long mag01[2]={0x0UL, RFWRandom::MATRIX_A};
			/* mag01[x] = x * MATRIX_A  for x=0,1 */

			if (mti >= RFWRandom::N) { /* generate N words at one time */
				int kk;

				if (mti == RFWRandom::N+1)   /* if init_genrand() has not been called, */
					init_genrand(5489UL); /* a default initial seed is used */

				for (kk=0;kk<RFWRandom::N-RFWRandom::M;kk++) {
					y = (mt[kk]&RFWRandom::UPPER_MASK)|(mt[kk+1]&RFWRandom::LOWER_MASK);
					mt[kk] = mt[kk+RFWRandom::M] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				for (;kk<RFWRandom::N-1;kk++) {
					y = (mt[kk]&RFWRandom::UPPER_MASK)|(mt[kk+1]&RFWRandom::LOWER_MASK);
					mt[kk] = mt[kk+(RFWRandom::M-RFWRandom::N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				y = (mt[RFWRandom::N-1]&RFWRandom::UPPER_MASK)|(mt[0]&RFWRandom::LOWER_MASK);
				mt[RFWRandom::N-1] = mt[RFWRandom::M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

				mti = 0;
			}
		  
			y = mt[mti++];

			/* Tempering */
			y ^= (y >> 11);
			y ^= (y << 7) & 0x9d2c5680UL;
			y ^= (y << 15) & 0xefc60000UL;
			y ^= (y >> 18);

			return y;
		}

	public:
		//constructors
		RFWLocalRandom () {Randomize(1);}
		RFWLocalRandom (unsigned long s) {Randomize(s);}


		void CriticalRandomize() {
//#pragma omp critical
			{
				//fprintf (stderr, "cr");
				Randomize((unsigned int)RFWRandom::getInteger(0,2000000000));
			}
		}

		void Randomize() {
			Randomize((unsigned int)RFWRandom::getInteger(0,2000000000));
		}

		//randomize procedures
		void Randomize (unsigned long s) {
			if (s==0) s = 1;
			init_genrand(s);
		}
		
		unsigned long GetRand() {return genrand_int32();}

		//pick an integer uniformly at random between inf and sup (both inclusive)
		int GetInteger (int inf, int sup) {
			if (sup<=inf) return inf;
			unsigned long range, minallowed, u;

			range = (unsigned long)(sup-inf+1);    //number of values allowed
			minallowed = (RFWRandom::maxvalue % range) + 1;   //restrict search space to avoid small numbers
			if (minallowed==range) minallowed = 0;
			do {u = GetRand();}     //repeat until a good number is found
			while (u < minallowed);

			return (inf + (int)(u % range)); //return a number in the range
		}



		float GetFloat () {return (float)GetDouble();} //get a float number in [0;1]
		double GetDouble() { return GetDoubleClosed(); } //double in the range [0;1]
		//double GetDouble() { double r = GetDoubleClosed(); fprintf(stderr, "<<< %.10f : %.10f : %d >>>", r, (double)RFWRandom::maxvalue, sizeof(RFWRandom::maxvalue));  return r; } //double in the range [0;1]
		double GetDoubleClosed() { return ((double)GetRand() / (double)RFWRandom::maxvalue); }  //double in the range [0;1]
		double GetDoubleOpen() {return ((double)GetRand()/((double)(RFWRandom::maxvalue)+1.0));} //double in the range [0;1)
		bool GetBool () {return (GetRand() & 1);}

};




#endif

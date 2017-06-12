/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/**********************************************************
 *
 * RFWTimer
 * - measures running time in second
 * - supports four different facilities (clock, rusage, htime, hrtimer) 
 * - define RFW_USAGE, RFW_CLOCK, RFW_HTIME, RFW_HRTIMER
 * - if none is the defined, a very crude heuristic is used
 *   to determine which one is to be used; don't rely on
 *   it too much
 **********************************************************/

#pragma once

//-----------------------------------------------------
// heuristic to determine which timing facility to use
//-----------------------------------------------------
#ifndef RFW_RUSAGE
#ifndef RFW_CLOCK
#ifndef RFW_HTIME
#ifndef RFW_HRTIMER
#ifndef RFW_STL

#ifdef _MSC_VER
	#define RFW_HRTIMER
#else 
	#define RFW_CLOCK
#endif

#endif
#endif
#endif
#endif 
#endif

#ifdef RFW_HRTIMER
	#include <Windows.h>
	#undef min
	#undef max
#endif


#ifdef RFW_RUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef RFW_HTIME 
#include <sys/time.h>
#endif

#ifdef RFW_CLOCK
#include <stdio.h>
#include <time.h>
#endif

#ifdef RFW_STL
#include <chrono>
#endif

//------------------
// the class itself
//------------------
class RFWTimer {
	private:
		bool running; //is it running now? (false -> paused)
		double base_time; //time of previous runs since last reset
		double max_time;  //reference time  

	#ifdef RFW_CLOCK 
		double getUserTime();
		double start_time;
	#endif

	#ifdef RFW_HRTIMER
		typedef struct {
			LARGE_INTEGER start;
			LARGE_INTEGER stop;
		} stopWatch;
		stopWatch timer;
		double frequency;
		//LARGE_INTEGER frequency;
	#endif

	#ifdef RFW_HTIME
		long start_time, end_time;
	#endif

	#ifdef RFW_RUSAGE
		struct rusage ru;
		struct timeval start_time, end_time, sample_time;
	#endif
		
	#ifdef RFW_CHRONO
		// Wall clock
		std::chrono::high_resolution_clock::time_point start_time;
	#endif
		
		void setBaseTime (double bt); 

		//facility-dependent functions
		double getElapsedTime(); //time since last resume/start
		void startTiming();      //store time for future comparison
		void init();

	public:
		//basic functions
		RFWTimer (bool start=false);
		double getTime(); //return current time
		double pause();   //pause and return current time
		double resume();  //continue if paused, start if reset; return time before resuming
		double start();   //reset and resume; return time before reset
		double reset();   //reset timer and pause (at zero); return time before reset

		//auxiliary functions
		void setMaxTime (double mt);
		double getMaxTime();
		double getTimeToExpire();
		bool isTimeExpired ();
};

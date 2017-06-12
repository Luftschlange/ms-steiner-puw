/*
Algorithm for Steiner Problem in Graphs

Copyright (c) Microsoft Corporation

All rights reserved. 

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*************************************
 *
 * RFWTimer (used for time measured)
 *
 *************************************/

#include "rfw_timer.h"

#define INFINITE_TIME 10e100

/********************
 * Generic functions
 ********************/

RFWTimer::RFWTimer (bool s) {
	init();
	base_time = 0.0;
	max_time = 0.0;
	if (s) start();
	else running = false;
}

double RFWTimer::getTime() {
	if (running) return (getElapsedTime() + base_time);
	else return base_time;
}

double RFWTimer::start() {
	double current = getTime();
	base_time = 0.0;
	startTiming();
	running = true;
	return current;
}

void RFWTimer::setMaxTime (double mt) {max_time = mt;}
double RFWTimer::getMaxTime () {return max_time;}

bool RFWTimer::isTimeExpired () {
	if (getMaxTime() == 0) return false;
	bool time_expired = (getTime() >= getMaxTime());
	return time_expired;
}

double RFWTimer::getTimeToExpire () { //may be negative!
	if (getMaxTime() == 0) return INFINITE_TIME;
	else return (getTime() - getMaxTime());
}

void RFWTimer::setBaseTime (double bt) {base_time = bt;}

double RFWTimer::reset() {
	double current = getTime();
	running = false;
	base_time = 0.0;
	return current;
}

double RFWTimer::pause() {
	base_time = getTime();
	running = false;
	return base_time;
}

double RFWTimer::resume() {
	if (running) return getTime();
	else {
		running = true;
		startTiming();
		return base_time;
	}
}

//------------
// RFW_RUSAGE
//------------

#ifdef RFW_RUSAGE

void RFWTimer::init(){};

void RFWTimer::startTiming() {
  getrusage(RUSAGE_SELF, &ru);
  start_time = ru.ru_utime;
}

double RFWTimer::getElapsedTime() {
  double t;
  getrusage(RUSAGE_SELF, &ru);
  end_time = ru.ru_utime;
  if (end_time.tv_usec < start_time.tv_usec){
    end_time.tv_usec += 1000000;
    end_time.tv_sec -= 1;
  }
  t = 100.0*(double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_usec - start_time.tv_usec) / (double)10000.0;
  return ((double)t/(double)100);
}

#endif


//-----------
// RFW_CLOCK
//-----------

#ifdef RFW_CLOCK

void RFWTimer::init(){};

void RFWTimer::startTiming() {start_time = getUserTime();}

double RFWTimer::getElapsedTime() {
	return (getUserTime() - start_time);
}

double RFWTimer::getUserTime() {
	double msecs = (double) clock() / CLOCKS_PER_SEC;
	if (msecs > 0) return msecs; 
	else return 0.0; //sometimes msecs is -0.000 (go figure...)
}

#endif


/*************
 *RFW_HTIME
 *************/

#ifdef RFW_HTIME

void RFWTimer::init() {};

void RFWTimer::startTiming() {start_time = (long)gethrvtime();}

double RFWTimer::getElapsedTime() {
	return (double)(gethrvtime()-start_time)/(double)10e8;
}

#endif


#ifdef RFW_HRTIMER

void RFWTimer::init() {
	LARGE_INTEGER f;
	QueryPerformanceFrequency(&f) ;
	frequency = (double)f.QuadPart;
}

void RFWTimer::startTiming() {
	QueryPerformanceCounter(&timer.start);
}

double RFWTimer::getElapsedTime() {
	QueryPerformanceCounter(&timer.stop);
	/*
	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;*/
	//LARGE_INTEGER frequency;
	//QueryPerformanceFrequency(&frequency) ;
	//return ((double)time.QuadPart / (double)frequency.QuadPart);
	return ((double)(timer.stop.QuadPart - timer.start.QuadPart) / frequency);
}

#endif


/*************
 *RFW_CHRONO
 *************/

#ifdef RFW_CHRONO

void RFWTimer::init() {};

void RFWTimer::startTiming() { start_time = std::chrono::high_resolution_clock::now(); }

double RFWTimer::getElapsedTime() {
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start_time;
	return elapsed.count();
}

#endif
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


class SteinerConfig {
public:
	double ELITE_DENOMINATOR;
	int FIRST_BRANCH;
	bool OUTPUT_INCUMBENT;
	int PERTURBATION_MODE;
	bool STRONG_BRANCHING;
	int SB_LOW;
	int SB_HIGH;
	int SB_GAP;
	int VERBOSE_DEPTH;
	double EARLY_STOP_BOUND;
	bool AGGRESSIVE_COMBINATION;
	int LS_PERT_ROUNDS; // maximum number of rounds of local search on perturbed instance
	double LS_PERT_EXPONENT; // additive perturbation p will be set to (p)^LS_PERT_DECAY_FACTOR
	double PERT_VERTEX; // vertex-centric perturbation
	double PERT_RANGE; //expected multiplier after perturbation
	double PERT_EXTRA;
	double PERT_FLOOR;
	double OUTPUT_THRESHOLD; //output solution if better than this value
	int ROOT_COMP_MODE; // selection mode for root components
	int MAX_COMB_FAIL;
	char LSTYPE[1024];
	int DEPTH_LIMIT;
	string LOG_FILENAME;
	double TIME_LIMIT; // the maximum time limit in seconds; 0 = none.
	int RESILIENT_PERTURBATION; //0:off 1:on 2:auto (true if using perturbation)
	
	SteinerConfig () {
		STRONG_BRANCHING = false;
		FIRST_BRANCH = 1;
		SB_LOW = 10;
		SB_HIGH = 18;
		SB_GAP= 250;
		VERBOSE_DEPTH = 20;
		OUTPUT_INCUMBENT = true;
		ELITE_DENOMINATOR = 2.0;
		PERTURBATION_MODE = 7;
		EARLY_STOP_BOUND = 0; //stop if you a find a solution with this value
		AGGRESSIVE_COMBINATION = true;
		LS_PERT_ROUNDS = 3;
		LS_PERT_EXPONENT = .5; //no decay!
		PERT_VERTEX = .5;
		PERT_FLOOR = 1.0;
		PERT_EXTRA = .75;
		PERT_RANGE = .25;
		OUTPUT_THRESHOLD = -1;
		MAX_COMB_FAIL = 3;
		ROOT_COMP_MODE = 2;
		sprintf(LSTYPE, "vq");
		DEPTH_LIMIT  = 0x7FFFFFFF;
		LOG_FILENAME = "";
		TIME_LIMIT = 0.0;
		RESILIENT_PERTURBATION = 2; //auto
	}

	void Output (FILE *file) {
		fprintf (file, "elitedenominator %f\n", ELITE_DENOMINATOR);
		fprintf (file, "firstbranch %d\n", FIRST_BRANCH);
		fprintf (file, "outputincumbent %d\n", (int)OUTPUT_INCUMBENT);
		fprintf (file, "pertmode %d\n", PERTURBATION_MODE);
		fprintf (file, "sb %d\n", STRONG_BRANCHING);
		fprintf (file, "sbgap %d\n", SB_GAP);
		fprintf (file, "sblow %d\n", SB_LOW);
		fprintf (file, "sbhigh %d\n", SB_HIGH);
		fprintf (file, "verbdepth %d\n", VERBOSE_DEPTH);
		fprintf (file, "earlystopbound %.0f\n", EARLY_STOP_BOUND); 
		fprintf (file, "aggressivecomb %d\n", AGGRESSIVE_COMBINATION);
		fprintf (file, "lspertrounds %d\n", LS_PERT_ROUNDS);
		fprintf (file, "lspertdecay %.2f\n", LS_PERT_EXPONENT);
		fprintf (file, "pertvertex %.4f\n", PERT_VERTEX);
		fprintf (file, "pertfloor %.4f\n", PERT_FLOOR);
		fprintf (file, "pertrange %.4f\n", PERT_RANGE);
		fprintf (file, "pertextra %.4f\n", PERT_EXTRA);
		fprintf (file, "outputthreshold %.4f\n", OUTPUT_THRESHOLD);
		fprintf (file, "lstype %s\n", LSTYPE);
		fprintf (file, "maxcombfail %d\n", MAX_COMB_FAIL);
		fprintf (file, "rootcompmode %d\n", ROOT_COMP_MODE);
		fprintf (file, "depthlimit %d\n", DEPTH_LIMIT);
		fprintf (file, "logfilename %s\n", LOG_FILENAME.c_str());
		fprintf (file, "timelimit %f\n", TIME_LIMIT);
		fprintf (file, "resilientperturbation %d\n", RESILIENT_PERTURBATION);
	}

	void ReadParameter (const char *key, const char *value) {
		if (strcmp(key,"-outputincumbent")==0) {OUTPUT_INCUMBENT = (atoi (value) != 0);}
		if (strcmp(key,"-sb")==0) {STRONG_BRANCHING = (atoi (value) != 0);}
		if (strcmp(key,"-sblow")==0) {SB_LOW = (atoi (value));}
		if (strcmp(key,"-sbhigh")==0) {SB_HIGH = (atoi (value));}
		if (strcmp(key,"-sbgap")==0) {SB_GAP = (atoi (value));}
		if (strcmp(key,"-verbdepth")==0) {VERBOSE_DEPTH = (atoi(value));}
		if (strcmp(key,"-firstbranch")==0) {FIRST_BRANCH = atoi(value);}
		if (strcmp(key,"-maxcombfail")==0) {MAX_COMB_FAIL = atoi(value);}
		if (strcmp(key,"-elitedenominator")==0) {ELITE_DENOMINATOR = atof(value);}
		if (strcmp(key,"-pertmode")==0) {PERTURBATION_MODE = atoi(value);}
		if (strcmp(key,"-earlystopbound")==0) {EARLY_STOP_BOUND = atof(value);}
		if (strcmp(key,"-aggressivecomb")==0) {AGGRESSIVE_COMBINATION = (atoi (value) != 0);}
		if (strcmp(key,"-lspertdecay")==0) {LS_PERT_EXPONENT = atof(value);}
		if (strcmp(key,"-lspertrounds")==0) {LS_PERT_ROUNDS = atoi(value);}
		if (strcmp(key,"-pertvertex")==0) {PERT_VERTEX = atof (value);}
		if (strcmp(key,"-pertfloor")==0) {PERT_FLOOR = atof(value);}		
		if (strcmp(key,"-pertrange")==0) {PERT_RANGE = atof(value);}
		if (strcmp(key,"-pertextra")==0) {PERT_EXTRA = atof(value);}
		if (strcmp(key,"-outputthreshold")==0) {OUTPUT_THRESHOLD = atof(value);}
		if (strcmp(key,"-lstype")==0) {strcpy(LSTYPE, value);}
		if (strcmp(key,"-rootcompmode")==0) {ROOT_COMP_MODE = atoi(value);}
		if (strcmp(key,"-depthlimit")==0) {DEPTH_LIMIT = atoi(value);}
		if (strcmp(key, "-logfilename") == 0) { LOG_FILENAME = value; }
		if (strcmp(key, "-timelimit") == 0) { TIME_LIMIT = atof(value); }
		if (strcmp(key, "-resilientperturbation")==0) {RESILIENT_PERTURBATION = atoi(value);}
	}

};


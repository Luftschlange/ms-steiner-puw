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
#include "graph.h"
#include <cmath>

class GraphDrawer {
private:
	static void fatal (const string &msg) {
		fprintf (stderr, "GraphDrawer::ERROR: %s.\n", msg.c_str());
		exit(-1);
	}

	typedef enum {GDE_STANDARD, GDE_SOLUTION, GDE_DELETED} EdgeType;
	typedef enum {GDV_STANDARD, GDV_DELETED} VertexType;


	static void OutputHeader(FILE *file) {


		fprintf (file, "add style vertex #vertex\n");
		fprintf (file, "set style vertex #size 1\n");
		fprintf (file, "set style vertex #shape 0\n");
		fprintf (file, "set style vertex #color_r 0.5\n");
		fprintf (file, "set style vertex #color_g 0.5\n");
		fprintf (file, "set style vertex #color_b 0.5\n");
		fprintf (file, "set style vertex #label_height 1\n");
		fprintf (file, "set style vertex #label_font Helvetica\n");
		fprintf (file, "add style delvertex vertex\n");
		fprintf (file, "set style delvertex #color_r 0.9\n");
		fprintf (file, "set style delvertex #color_g 0.9\n");
		fprintf (file, "set style delvertex #color_b 0.9\n");
		fprintf (file, "proxy vertex #label #name string\n");
		fprintf (file, "add style terminal vertex\n");
		fprintf (file, "set style terminal #size 2\n");
		fprintf (file, "set style terminal #shape 1\n");
		fprintf (file, "add style edge #edge\n");
		fprintf (file, "set style edge #width 0.25\n");
		fprintf (file, "set style edge #color_r 0.75\n"); 
		fprintf (file, "set style edge #color_g 0.75\n"); 
		fprintf (file, "set style edge #color_b 0.75\n"); 
	}

	template <class GRAPH> static double FindScalingFactor(GRAPH &g,vector<EdgeType> &etype) {
		fprintf (stderr, "Finding scaling factor.\n");

		vector<double> lengths;

		int m = g.EdgeCount();
		for (int e=1; e<=m; e++) {
			if (etype[e]==GDE_DELETED) continue;
			int v,w;
			g.GetEndpoints(e,v,w);
			
			double vx,vy,wx,wy;
			g.coord.GetCoordinates(v,vx,vy);
			g.coord.GetCoordinates(w,wx,wy);

			//fprintf (stderr, "%.2f %.2f\n", vx, vy);

			double deltax = (wx-vx);
			double deltay = (wy-vy);

			double sqdist = deltax*deltax + deltay*deltay;
			lengths.push_back(sqdist);
		}

		/*
		for (int i=0; i<m; i++) {
			fprintf (stderr, "%.2f ", lengths[i]);
		}*/


		sort (lengths.begin(), lengths.end());

		double refpos = (double)lengths.size() / 5.0; //(int)floor(sqrt((double)lengths.size())); // / 20;
		double reflen = sqrt(lengths[(int)refpos]);

		// scaling factor is 1 if reference distance is 5
		// (arbitrary thing, tuned for vlsi)
		double scale = 5.0 / reflen;

		fprintf (stderr, "Reference length is %.2f, using scaling factor %.5f\n", reflen, scale);

		return scale; 
	}

	template <class GRAPH> static void OutputVertices(FILE *file, GRAPH &g, vector<EdgeType> &etype, vector<VertexType> &vtype) {
		int n = g.VertexCount();
		double scale = FindScalingFactor(g,etype);
		for (int i=0; i<2; i++) {
			if (i==0) {
				fprintf (file, "use style vertex\n");
			} else {
				fprintf (file, "use style terminal\n");
			}
			for (int v=1; v<=n; v++) {
				if (g.IsTerminal(v) == (i!=0)) { 
					double x, y;
					g.coord.GetCoordinates(v,x,y);
					x *= scale;
					y *= scale;
					if (i==0) {
						if (vtype[v] == GDV_DELETED) fprintf (file, "use style delvertex\n");
						else if (vtype[v] == GDV_STANDARD) fprintf (file, "use style vertex\n");
					}
					fprintf (file, "add vertex %d %.0f %.0f\n", v, x, y);
				}
			}
		}
	}

	template <class GRAPH> static void OutputEdges(FILE *file, GRAPH &g, vector<EdgeType> &edgeinfo) {
		int m = g.EdgeCount();
		fprintf (file, "use style edge\n");
		for (int e=1; e<=m; e++) {
			if (edgeinfo[e]!=GDE_STANDARD) continue;
			int v,w;
			g.GetEndpoints(e,v,w);
			fprintf (file, "add edge %d %d %d \n", e, v, w);
		}
	}

	template <class GRAPH> void MarkAllEdges (GRAPH &g, vector<int> &etype, EdgeType t) {
		int m = g.EdgeCount();


	}

public:

	template <class GRAPH> static void DrawSolution (FILE *file, GRAPH &g, vector<EdgeType> &edgeinfo, vector<VertexType> &vinfo) {
		OutputHeader(file);
		OutputVertices(file,g,edgeinfo,vinfo);
		OutputEdges(file,g,edgeinfo);
	}

	template <class GRAPH> static void DrawGraph (char *filename, GRAPH &g) {
		FILE *file = fopen (filename, "w");
		if (!file) fatal ("could not open filename for drawing");
		int m = g.EdgeCount();
		int n = g.VertexCount();
		vector<EdgeType> etype(m+1,GDE_STANDARD);
		vector<VertexType> vtype(n+1,GDV_STANDARD);
		DrawSolution (file, g, etype, vtype);
		fclose(file);
	}

	template <class GRAPH> static void DrawSubgraph (const string &filename, GRAPH &g, vector<bool> &emarked, vector<bool> &vkept) {
		FILE *file = fopen (filename.c_str(), "w");
		if (!file) fatal ("could not open filename for drawing");
		int m = g.EdgeCount();
		int n = g.VertexCount();
		vector<EdgeType> etype(m+1,GDE_STANDARD);
		vector<VertexType> vtype(n+1,GDV_STANDARD);
		for (int e=1; e<=m; e++) {
			if (!emarked[e]) etype[e] = GDE_DELETED;
		}
		for (int v=1; v<=n; v++) {
			if (!vkept[v]) vtype[v] = GDV_DELETED;
		}

		DrawSolution (file, g, etype, vtype);
		fclose(file);
	}

};

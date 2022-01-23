/*
MIT License

Copyright (c) 2016, Hao Wei.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <queue>
#include <algorithm>
#include <utility>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <chrono>

#include "Util.h"
#include "UnitHeap.h"

namespace Gorder
{

using namespace std;

class Vertex{
public:
	int outstart;
	int outdegree;
	int instart;
	int indegree;

	Vertex(){
		outdegree=indegree=0;
		outstart=instart=-1;
	}
};

class Graph{
	public:
		int vsize;
    int vnum;
		long long edgenum;
		string name;
		
		vector<Vertex> graph;
		vector<int> outedge;
		vector<int> inedge;
	
		string getFilename();
		void setFilename(string name);

		Graph();
		~Graph();
		void clear();
		void readGraph(const string& fullname);
		void writeGraph(ostream&);
		void PrintReOrderedGraph(const vector<int>& order);
		void GraphAnalysis();
		void RemoveDuplicate(const string& fullname);
    void SubGraphTest(Graph& sub, const vector<int>& candidate, double p);
    void SubGraphTest2(Graph& sub, const vector<int>& candidate);
    void SubGraph(Graph& sub, const vector<int>& candidate);
    void WriteSampleOriginalGraph(set<int>& visited, ofstream& out);
    void WriteSampleGraph(set<int>& visited, vector<int>& retorder, ofstream& out);
    void WriteSampleRandomGraph(set<int>& visited, ofstream& out);
    //void SubGraph(Graph& sub, const set<int>& candidate);
		
		void strTrimRight(string& str);
		static vector<string> split(const string &s, char delim);
		static vector<string>& split(const string &s, char delim, vector<string> &elems);

		void GapCount();
		double GapCost(vector<int>& order);
		double GapCostV(vector<int>& order, set<int>& visited);
		//void Transform();
		vector<int> Transform();
		void NDGorderGreedy(vector<int>& order, int window, vector<int>& candidate);
		void GorderGreedy(vector<int>& order, int window, vector<int>& candidate);
		void GorderSubGreedy(vector<int>& order, int window, vector<int>& candidate);
		void GorderSubGreedy2(vector<int>& order, int window, vector<int>& candidate, float p);
		void GorderTestSubGreedy(vector<int>& order, int window, vector<int>& candidate, vector<int>& collected);

		void RCMOrder(vector<int>& retorder);
    void ReRCMOrder(vector<int>& reorder);
		unsigned long long LocalityScore(const int w);
};

}

#endif


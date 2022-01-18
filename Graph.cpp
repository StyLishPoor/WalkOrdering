#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <unordered_set>
#include <functional>
#include <climits>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <random>

#ifdef __GNUC__
#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#else
#define likely(cond) (!!(cond))
#define unlikely(cond) (!!(cond))
#endif // GNUC

namespace Gorder
{

string Graph::getFilename(){
	return name;
}

void Graph::setFilename(string name){
	this->name.swap(name);
}

Graph::Graph() {
	edgenum=vsize=0;
}

Graph::~Graph() {
}

void Graph::clear() {
	vsize = 0;
	edgenum=0;
	name.clear();
	graph.clear();
	outedge.clear();
	inedge.clear();
}

void Graph::SubGraphTest2(Graph& sub, const vector<int>& candidate, double p) {
  clock_t sub_start = clock();
  vector<bool> exist(vsize, false);
  for (const int v : candidate) {
    exist[v] = true;
  }
  //random_device seed_gen;
  //mt19937 engine {seed_gen()};
  //uniform_real_distribution<double> dist(0, 1);
  vector< pair<int, int> > edges;
  edges.reserve(10000000);
  for (size_t u = 0; u < vsize; u++) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[u]) {
        edges.push_back(make_pair(u, v));
      } else {
        if(exist[v]) {
          edges.push_back(make_pair(u, v));
        }
      }
    }
  }
  cout << "Edge Num: " << edges.size() << endl;
  sub.edgenum = edges.size();
  sub.vsize = vsize;
  sub.vnum = vsize;
  sub.graph.resize(vsize+1);
  for(long long i=0; i<edges.size(); i++){
    sub.graph[edges[i].first].outdegree++;
    sub.graph[edges[i].second].indegree++;
  }
  sub.graph[0].outstart=0;
  sub.graph[0].instart=0;
   for(int i=1; i<vsize; i++){
    sub.graph[i].outstart=sub.graph[i-1].outstart+sub.graph[i-1].outdegree;
    sub.graph[i].instart=sub.graph[i-1].instart+sub.graph[i-1].indegree;
  }
 
  sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
   if(a.first<b.first) return true;
   else if(a.first>b.first) return false;
   else {
     if(a.second<=b.second) return true;
     else return false;
   }
  });

  sub.outedge.resize(sub.edgenum);
  sub.inedge.resize(sub.edgenum);
  for(long long i=0; i<edges.size(); i++){
    sub.outedge[i]=edges[i].second;
    sub.inedge[i]=edges[i].second;
  }
  clock_t sub_end = clock();
  cout << "Sub Graph Construce: " << (double)(sub_end - sub_start)/CLOCKS_PER_SEC << endl;
}

void Graph::SubGraphTest(Graph& sub, const vector<int>& candidate, double p) {
  //clock_t sub_start = clock();
  vector<bool> exist(vsize, false);
  vector<bool> added_flag(vsize, false);
  //vector<int> added_candidate;
  vector<int> neighbors;
  //added_candidate.reserve(vsize);
  neighbors.reserve(vsize);

  random_device seed_gen;
  mt19937 engine {seed_gen()};
  uniform_real_distribution<double> dist(0, 1);

  vector< pair<int, int> > edges;
  edges.reserve(10000000);

  /* この中は正常に動く
  //for (const int u : candidate) {
  //  exist[u] = true;
  //  added_candidate.push_back(u);
  //}

  //for (const int u : candidate) {
  //  for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
  //    int v = outedge[i];
  //    if (exist[v] == false && added_flag[v] == false) {
  //      neighbors.push_back(v);
  //      added_flag[v] = true;
  //    }
  //  }
  //}
  //
  //for (const int n : neighbors) {
  //  if (dist(engine) <= p) {
  //    added_candidate.push_back(n);
  //    exist[n] = true;
  //  }
  //}
  */

  //for (const int u : added_candidate) {
  //  for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
  //    int v = outedge[i];
  //    if(exist[v]) edges.push_back(make_pair(u, v));
  //  }
  //}

  for (const int u : candidate) {
    exist[u] = true;
  }

  // ノード集合を選んでから全部にエッジを張るパターン
  for (const int u : candidate) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[v]) {
        edges.push_back(make_pair(u, v));
      } else {
        if (!added_flag[v]) {
          added_flag[v] = true;
          neighbors.push_back(v);
        }
      }
    }
  }

  for (const int n : neighbors) {
    if (dist(engine) <= p) {
      for (size_t i = graph[n].outstart; i < graph[n].outstart + graph[n].outdegree; i++) {
        int v = outedge[i];
        if (exist[v]) {
          edges.push_back(make_pair(n, v));
          edges.push_back(make_pair(v, n));
        }
      }
    }
  }
  

  // cand-candは必ずエッジを張り，cand-cand外は確率pでエッジを張るパターン
  /*
  for (const int u : candidate) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[v]) {
        edges.push_back(make_pair(u, v));
      } else {
        if (dist(engine) <= p) {
          edges.push_back(make_pair(u, v));
          edges.push_back(make_pair(v, u));
        }
      }
    }
  }
  */

  sub.edgenum = edges.size();
  sub.vsize = vsize;
  sub.vnum = vsize;
  //sub.vnum = added_candidate.size();
  sub.graph.resize(vsize+1);
  for(long long i=0; i<edges.size(); i++){
    sub.graph[edges[i].first].outdegree++;
    sub.graph[edges[i].second].indegree++;
  }
  sub.graph[0].outstart=0;
  sub.graph[0].instart=0;
  for(int i=1; i<vsize; i++){
    sub.graph[i].outstart=sub.graph[i-1].outstart+sub.graph[i-1].outdegree;
    sub.graph[i].instart=sub.graph[i-1].instart+sub.graph[i-1].indegree;
  }
 
  sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
   if(a.first<b.first) return true;
   else if(a.first>b.first) return false;
   else {
     if(a.second<=b.second) return true;
     else return false;
   }
  });

  sub.outedge.resize(sub.edgenum);
  sub.inedge.resize(sub.edgenum);
  for(long long i=0; i<edges.size(); i++){
    sub.outedge[i]=edges[i].second;
    sub.inedge[i]=edges[i].second;
  }
  //clock_t sub_end = clock();
  //cout << "Sub Graph Construce: " << (double)(sub_end - sub_start)/CLOCKS_PER_SEC << endl;
}

void Graph::SubGraph(Graph& sub, const vector<int>& candidate) {
  clock_t sub_start = clock();
  vector<bool> exist(vsize, false); // exist[i]==trueのときiはcandidateに含まれる
  vector< pair<int, int> > edges;
	edges.reserve(10000000);
  for (const int u : candidate) {
    exist[u] = true;
  }

  sub.edgenum=0;
  for (const int u : candidate) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
        edges.push_back(make_pair(u, v));
    }
  }
  sub.edgenum = edges.size();
  sub.vsize = vsize;
  sub.vnum = candidate.size();
  sub.graph.resize(vsize+1);
  for(long long i=0; i<edges.size(); i++){
    sub.graph[edges[i].first].outdegree++;
    sub.graph[edges[i].second].indegree++;
  }
  sub.graph[0].outstart=0;
  sub.graph[0].instart=0;
   for(int i=1; i<vsize; i++){
    sub.graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
    sub.graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
  }
 
  sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
   if(a.first<b.first) return true;
   else if(a.first>b.first) return false;
   else {
     if(a.second<=b.second) return true;
     else return false;
   }
  });

  sub.outedge.resize(edgenum);
  sub.inedge.resize(edgenum);
  for(long long i=0; i<edges.size(); i++){
    sub.outedge[i]=edges[i].second;
    sub.inedge[i]=edges[i].second;
  }
  clock_t sub_end = clock();
  cout << "Sub Graph Construce: " << (double)(sub_end - sub_start)/CLOCKS_PER_SEC << endl;
}

void Graph::readGraph(const string& fullname) {
	FILE* fp;
	fp=fopen(fullname.c_str(), "r");
	if(fp==NULL){
		cout << "Fail to open " << fullname << endl;
		quit();
	}

	char line[40];
	int u, v;
	const char* str=NULL;
	
	vsize=0;
	edgenum=0;
	vector< pair<int, int> > edges;
	edges.reserve(100000000);

	while(feof(fp)!=true){
		if(fgets(line, 40, fp)){
			u=v=0;
			str=line;
			while(isdigit(*str))
				u=u*10+(*str++ - '0');
			str++;
			while(isdigit(*str))
				v=v*10+(*str++ - '0');

			if(u==v)
				continue;
			edgenum++;
			if(u>vsize)
				vsize=u;
			if(v>vsize)
				vsize=v;

			edges.push_back(make_pair(u, v));
		}
	}
	vsize++;

	fclose(fp);
	graph.resize(vsize+1);
	for(long long i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}
	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	
	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}

	});
	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}

	vector< pair<int, int> >().swap(edges);
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;
}

//void Graph::Transform(){
vector<int> Graph::Transform(){
	vector<int> order;
	RCMOrder(order);
  vector<int> original_order(order.size());
  for (size_t i = 0; i < order.size(); i++) {
    original_order[order[i]] = i;
  }
	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}
	if(graph.size()!=(vsize+1)){
		cout << "graph.size()!=(vsize+1)" << endl;
		quit();
	}

	vector<int>().swap(inedge);
	vector< pair<int, int> > edges;
	edges.reserve(edgenum);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart, limit=graph[i+1].outstart; j<limit; j++)
			edges.push_back(make_pair(order[i], order[outedge[j]]));
	}
	if(edges.size()!=edgenum){
		cout << "edges.size()!=edgenum" << endl;
		quit();
	}

	for(int i=0; i<vsize; i++){
		graph[i].outdegree=graph[i].indegree=0;
	}
	for(int i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}

	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;

	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{
			if(a.second<=b.second)
				return true;
			else
				return false;
		}
	});

	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}
	vector< pair<int, int> >().swap(edges);
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
  return original_order;
}


void Graph::writeGraph(ostream& out){
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outdegree+graph[u].outstart; j++){
			int v=outedge[j];
			out << u << '\t' << v << endl;
		}
	}
}

void Graph::WriteSampleOriginalGraph(set<int>& visited, ofstream& out) {
  vector<bool> exist(vsize, false);
  //vector< pair<int, int> > edges;
	//edges.reserve(10000000);
  for (const int u : visited) {
    exist[u] = true;
  }

  for (const int u : visited) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[v]) {
        out << u << " " << v << endl;
      }
    }
  }
}

void Graph::WriteSampleGraph(set<int>& visited, vector<int>& retorder, ofstream& out) {
  vector<bool> exist(vsize, false);
  //vector< pair<int, int> > edges;
	//edges.reserve(10000000);
  for (const int u : visited) {
    exist[u] = true;
  }

  for (const int u : visited) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[v]) {
        out << retorder[u] << " " << retorder[v] << endl;
      }
    }
  }
}

void Graph::WriteSampleRandomGraph(set<int>& visited, ofstream& out) {
  vector<bool> exist(vsize, false);
  vector<int> random_retorder(vsize);
  random_device seed_gen;
  mt19937 engine(seed_gen());
  for (size_t i = 0; i <= vsize; i++) {
    random_retorder[i] = i;
  }

  shuffle(random_retorder.begin(), random_retorder.end(), engine);

  for (const int u : visited) {
    exist[u] = true;
  }

  for (const int u : visited) {
    for (size_t i = graph[u].outstart; i < graph[u].outstart + graph[u].outdegree; i++) {
      int v = outedge[i];
      if (exist[v]) {
        out << random_retorder[u] << " " << random_retorder[v] << endl;
      }
    }
  }
}


void Graph::PrintReOrderedGraph(const vector<int>& order){
	ofstream out((name+"_Gorder.txt").c_str());

	vector<int>().swap(inedge);

	vector< vector<int> > ReOrderedGraph(vsize);
	int u, v;
	for(int i=0; i<vsize; i++){
		u=order[i];
		ReOrderedGraph[u].reserve(graph[i+1].outstart-graph[i].outstart);
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			v=order[outedge[j]];
			ReOrderedGraph[u].push_back(v);
		}
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
/*
	for(int u=0; u<vsize; u++){
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
*/
	for(int u=0; u<vsize; u++){
		for(int j=0; j<ReOrderedGraph[u].size(); j++){
			out << u << '\t' << ReOrderedGraph[u][j] << endl;
		}
	}
	out.close();
}


void Graph::GraphAnalysis(){
	vector<int> tmp(vsize);
	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].outdegree;
	}
	sort(tmp.begin(), tmp.end());
	
	cout << "outdegree:" << endl;
	vector<int>::iterator tmpit1, tmpit2;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}


	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].indegree;
	}
	sort(tmp.begin(), tmp.end());
	
	cout << "indegree:" << endl;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}
}


void Graph::RemoveDuplicate(const string& fullname) {
	FILE* fp;
	fp=fopen(fullname.c_str(), "r");
	if(fp==NULL){
		cout << "Fail to open " << fullname << endl;
		quit();
	}

	char line[40];
	int u, v;
	const char* str=NULL;
	
	set< pair<int, int> > edges;

	while(feof(fp)!=true){
		if(fgets(line, 40, fp)){
			u=v=0;
			str=line;
			while(isdigit(*str))
				u=u*10+(*str++ - '0');
			str++;
			while(isdigit(*str))
				v=v*10+(*str++ - '0');

			if(u==v)
				continue;

			edges.insert(make_pair(u, v));
		}
	}
	
	fclose(fp);

	cout << "after remove, the size is " << edges.size() << endl;
	
	ofstream fout;
	fout.open("NoDuplicate.txt");
	for(set< pair<int, int> >::iterator it=edges.begin(); it!=edges.end(); it++){
		fout << it->first << '\t' << it->second << endl;
	}
	fout.close();
}


vector<string>& Graph::split(const string &s, char delim, vector<string> &elems) {
	int begin, end;

	begin=0;
	end=s.find(delim);
	while(end!=string::npos){
		elems.push_back(s.substr(begin, end-begin));
		begin=end+1;
		end=s.find(delim, begin);
	}
	if(begin!=s.size()){
		elems.push_back(s.substr(begin));
	}

	return elems;
}

vector<string> Graph::split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

void Graph::strTrimRight(string& str) {
	string whitespaces(" \t\r\n");
	int index = str.find_last_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(index+1);
	else
		str.clear();
}

void Graph::GapCount(){
	int* gap = new int[vsize];
	memset(gap, 0, sizeof(int)*vsize);

	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
				gap[outedge[j]-outedge[j-1]]++;
		}
		gap[outedge[graph[i].outstart]]++;
	}

	double entropy=0;
	for(int i=0; i<vsize; i++){
		if(gap[i]==0)
			continue;
		else{
			entropy+=(double)gap[i]/edgenum*log2((double)gap[i]/edgenum);
		}
	}
	cout << "shannon: " << entropy << endl;
	delete[] gap;
}

/*
double Graph::GapCost(vector<int>& order){
	double gaplog=0;
	double gaplog2=0;
	vector<int> edgelist;
	edgelist.reserve(100000);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
			if(outedge[j]-outedge[j-1])
				gaplog+=log(double(outedge[j]-outedge[j-1]))/log(double(2));
		}
		edgelist.clear();
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			edgelist.push_back(order[outedge[j]]);
		}
		sort(edgelist.begin(), edgelist.end());
		for(int j=1; j<edgelist.size(); j++){
			if(edgelist[j]-edgelist[j-1])
				gaplog2+=log(double(edgelist[j]-edgelist[j-1]))/log(double(2));
		}
	}
	//cout << "original average gap cost: " << gaplog/edgenum << endl;
	//cout << "new average gap cost: " << gaplog2/edgenum << endl;

	return gaplog2/edgenum;
}
*/

double Graph::GapCostV(vector<int>& order, set<int>& visited){
  double gaplog=0;
  double gaplog2=0;
  vector<bool> exist(vsize, false);
  for (const int v : visited) {
    exist[v] = true;
  }
  int edge_num=0; // 集めたグラフだけのエッジ数
  vector<int> edgelist;
  edgelist.reserve(100000);
  for(const int i : visited){
    edgelist.clear();
    for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
      if (exist[outedge[j]]) {
        edge_num++;
        edgelist.push_back(order[outedge[j]]);
      }
		}
		sort(edgelist.begin(), edgelist.end());
		for(int j=1; j<edgelist.size(); j++){
			if(edgelist[j]-edgelist[j-1])
				gaplog2+=log(double(edgelist[j]-edgelist[j-1]))/log(double(2));
		}
	}
	//cout << "original average gap cost: " << gaplog/edgenum << endl;
	//cout << "new average gap cost: " << gaplog2/edge_num << endl;

	return gaplog2/edge_num;
}

void Graph::GorderTestSubGreedy(vector<int>& order, int window, vector<int>& candidate, vector<int>& collected){
  clock_t start=clock();
  UnitHeap unitheap(vsize, candidate);
  vector<bool> popvexist(vsize, false);
  vector<bool> update_active(vsize, false);
  vector<bool> use(vsize, false);
  int count=0;
  int finish_num = candidate.size() - 1;
  int tmpindex, tmpweight;
  const int hugevertex=sqrt((double)vsize);
  //const int hugevertex=sqrt((double)vnum);
  tmpweight=-1;
  for(const auto i : candidate){
    update_active[i] = true;
    unitheap.LinkedList[i].key=graph[i].indegree;
    unitheap.update[i]=-graph[i].indegree;
    if(graph[i].indegree>tmpweight) {
      tmpweight=graph[i].indegree;
      tmpindex=i;
    }
  }
  cout << "Cand: " << candidate.size() << " Coll: " << collected.size() << endl;
  for (const auto i : collected) {
    use[i] = true;
  }
  unitheap.ReConstruct(candidate);
  order.push_back(tmpindex);
  unitheap.update[tmpindex]=INT_MAX/2;
  unitheap.DeleteElement(tmpindex);
  // Initial Increment
  for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
    int u=inedge[i];
    if (use[u]) {
      if(graph[u].outdegree<=hugevertex){
        if(unitheap.update[u]==0){
          if(update_active[u]) unitheap.IncrementKey(u);
        } else {
          if(update_active[u]) unitheap.update[u]++;
        }
        if(graph[u].outdegree>1)
        for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
          int w=outedge[j];
          if(unitheap.update[w]==0){
            if(update_active[w]) unitheap.IncrementKey(w);
          } else {
            if(update_active[w]) unitheap.update[w]++;
          }
        }
      }
    }
  }
  if(graph[tmpindex].outdegree<=hugevertex){
    for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
      int w=outedge[i];
      if(unitheap.update[w]==0){
        if(update_active[w]) unitheap.IncrementKey(w);
      }else{
        if(update_active[w]) unitheap.update[w]++;
      }
    }
  }
  while(count<finish_num){
    int v=unitheap.ExtractMax();
    count++;
    order.push_back(v);
    unitheap.update[v]=INT_MAX/2;
    int popv;
    if(count-window>=0)
      popv=order[count-window];
    else
      popv=-1;
    // Decrement
		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					if(update_active[w]) unitheap.update[w]--;
				}
			}
			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
        if (use[u]) {
          if(graph[u].outdegree<=hugevertex){
            if(update_active[u]) unitheap.update[u]--;
            if(graph[u].outdegree>1)
            if(binary_search(outedge.data() + graph[u].outstart, outedge.data() + graph[u+1].outstart, v)==false){
              for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
                int w=outedge[j];
                if(update_active[w]) unitheap.update[w]--;
              }
            } else {
              popvexist[u]=true;
            }
          }
        }
			}
		}
    // Increment
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(unlikely(unitheap.update[w]==0)){
					if(update_active[w]) unitheap.IncrementKey(w);
				} else {
					if(update_active[w]) unitheap.update[w]++;
				}
			}
		}
		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
      if(use[u]) {
        if(graph[u].outdegree<=hugevertex){
          if(unlikely(unitheap.update[u]==0)){
            if(update_active[u]) unitheap.IncrementKey(u);
          } else {
            if(update_active[u]) unitheap.update[u]++;
          }
          if(popvexist[u]==false){
            if(graph[u].outdegree>1)
            for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
              int w=outedge[j];
              if(unlikely(unitheap.update[w]==0)){
                if(update_active[w]) unitheap.IncrementKey(w);
              }else{
                if(update_active[w]) unitheap.update[w]++;
              }
            }
          } else {
            popvexist[u]=false;
          }
        }
      }
		}
	}
  clock_t end=clock();
  cout << "Greedy: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
}

void Graph::GorderSubGreedy(vector<int>& order, int window, vector<int>& candidate){
  //clock_t start=clock();
  UnitHeap unitheap(vsize, candidate);
  vector<bool> popvexist(vsize, false);
  int count=0;
  int finish_num = candidate.size() - 1;
  int tmpindex, tmpweight;
  //const int hugevertex=sqrt((double)vsize);
  const int hugevertex=sqrt((double)vnum);
  tmpweight=-1;
  for(const auto i : candidate){
    unitheap.LinkedList[i].key=graph[i].indegree;
    unitheap.update[i]=-graph[i].indegree;
    if(graph[i].indegree>tmpweight) {
      tmpweight=graph[i].indegree;
      tmpindex=i;
    }
  }
  unitheap.ReConstruct(candidate);
  order.push_back(tmpindex);
  unitheap.update[tmpindex]=INT_MAX/2;
  unitheap.DeleteElement(tmpindex);
  // Initial Increment
  for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
    int u=inedge[i];
    if(graph[u].outdegree<=hugevertex){
      if(unitheap.update[u]==0){
        unitheap.IncrementKey(u);
      } else {
        unitheap.update[u]++;
      }
      if(graph[u].outdegree>1)
      for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
        int w=outedge[j];
        if(unitheap.update[w]==0){
          unitheap.IncrementKey(w);
        } else {
          unitheap.update[w]++;
        }
      }
    }
  }
  if(graph[tmpindex].outdegree<=hugevertex){
    for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
      int w=outedge[i];
      if(unitheap.update[w]==0){
        unitheap.IncrementKey(w);
      }else{
        unitheap.update[w]++;
      }
    }
  }
  while(count<finish_num){
    int v=unitheap.ExtractMax();
    count++;
    order.push_back(v);
    unitheap.update[v]=INT_MAX/2;
    int popv;
    if(count-window>=0)
      popv=order[count-window];
    else
      popv=-1;
    // Decrement
		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					unitheap.update[w]--;
				}
			}
			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
				if(graph[u].outdegree<=hugevertex){
					unitheap.update[u]--;
					if(graph[u].outdegree>1)
					if(binary_search(outedge.data() + graph[u].outstart, outedge.data() + graph[u+1].outstart, v)==false){
						for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
							int w=outedge[j];
							unitheap.update[w]--;
						}
					} else {
						popvexist[u]=true;
					}
				}
			}
		}
    // Increment
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(unlikely(unitheap.update[w]==0)){
					unitheap.IncrementKey(w);
				} else {
					unitheap.update[w]++;
				}
			}
		}
		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
			if(graph[u].outdegree<=hugevertex){
				if(unlikely(unitheap.update[u]==0)){
					unitheap.IncrementKey(u);
				} else {
					unitheap.update[u]++;
				}
				if(popvexist[u]==false){
					if(graph[u].outdegree>1)
					for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
						int w=outedge[j];
						if(unlikely(unitheap.update[w]==0)){
							unitheap.IncrementKey(w);
						}else{
							unitheap.update[w]++;
						}
					}
				} else {
					popvexist[u]=false;
				}
			}
		}
	}
  //clock_t end=clock();
  //cout << "Greedy: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
}

void Graph::GorderGreedy(vector<int>& order, int window, vector<int>& candidate){
  UnitHeap unitheap(vsize, candidate);
  vector<bool> popvexist(vsize, false);
  int count=0;
  int finish_num;
  int tmpindex, tmpweight;
  const int hugevertex=sqrt((double)vsize);
  for(const auto i : candidate){
      unitheap.LinkedList[i].key=graph[i].indegree;
      unitheap.update[i]=-graph[i].indegree;
  }
  unitheap.ReConstruct(candidate);
  finish_num = candidate.size() - 1;
  tmpweight=-1;
  for (const int i : candidate) {
    if(graph[i].indegree>tmpweight) {
      tmpweight=graph[i].indegree;
      tmpindex=i;
    }
  }
  order.push_back(tmpindex);
  unitheap.update[tmpindex]=INT_MAX/2;
  unitheap.DeleteElement(tmpindex);
  for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
    int u=inedge[i];
    if(graph[u].outdegree<=hugevertex){
      if(unitheap.update[u]==0){
        unitheap.IncrementKey(u);
      } else {
        unitheap.update[u]++;
      }
      if(graph[u].outdegree>1)
      for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
        int w=outedge[j];
        if(unitheap.update[w]==0){
          unitheap.IncrementKey(w);
        } else {
          unitheap.update[w]++;
        }
      }
    }
  }
  if(graph[tmpindex].outdegree<=hugevertex){
    for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
      int w=outedge[i];
      if(unitheap.update[w]==0){
        unitheap.IncrementKey(w);
      }else{
        unitheap.update[w]++;
      }
    }
  }
	while(count<finish_num){
    
		int v=unitheap.ExtractMax();
		count++;
		order.push_back(v);
		unitheap.update[v]=INT_MAX/2;

		int popv;
    if(count-window>=0)
      popv=order[count-window];
    else
      popv=-1;

		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					unitheap.update[w]--;
				}
			}
			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
				if(graph[u].outdegree<=hugevertex){
					unitheap.update[u]--;
					if(graph[u].outdegree>1)
					if(binary_search(outedge.data() + graph[u].outstart, outedge.data() + graph[u+1].outstart, v)==false){
						for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
							int w=outedge[j];
							unitheap.update[w]--;
						}
					} else {
						popvexist[u]=true;
					}
				}
			}
		}
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(unlikely(unitheap.update[w]==0)){
					unitheap.IncrementKey(w);
				} else {
					unitheap.update[w]++;
				}
			}
		}
		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
			if(graph[u].outdegree<=hugevertex){
				if(unlikely(unitheap.update[u]==0)){
					unitheap.IncrementKey(u);
				} else {
					unitheap.update[u]++;
				}
				if(popvexist[u]==false){
					if(graph[u].outdegree>1)
					for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
						int w=outedge[j];
						if(unlikely(unitheap.update[w]==0)){
							unitheap.IncrementKey(w);
						}else{
							unitheap.update[w]++;
						}
					}
				} else {
					popvexist[u]=false;
				}
			}
		}
	}
}

/*
void Graph::GorderGreedy(vector<int>& order, int window, vector<int>& candidate, int start_index, int end_index){
  UnitHeap unitheap(vsize, candidate);
  vector<bool> popvexist(vsize, false);
  vector<int> init_nodes;
  int count=0;
  int finish_num;
  int tmpindex, tmpweight;
  //const int hugevertex=1*sqrt((double)vsize);
  const int hugevertex=0.25*sqrt((double)vsize);
  for(const auto i : candidate){
      unitheap.LinkedList[i].key=graph[i].indegree;
      unitheap.update[i]=-graph[i].indegree;
  }
  unitheap.ReConstruct(candidate);
  finish_num = candidate.size() - 1;
  tmpweight=-1;
  for (const int i : candidate) {
    if(graph[i].indegree>tmpweight) {
      tmpweight=graph[i].indegree;
      tmpindex=i;
    }
  }
  order.push_back(tmpindex);
  init_nodes.push_back(tmpindex);
  unitheap.update[tmpindex]=INT_MAX/2;
  unitheap.DeleteElement(tmpindex);
  for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
    int u=inedge[i];
    if(graph[u].outdegree<=hugevertex){
      if(unitheap.update[u]==0){
        unitheap.IncrementKey(u);
        //if(unitheap.LinkedList[u].active) unitheap.IncrementKey(u);
      } else {
        unitheap.update[u]++;
        //if(unitheap.LinkedList[u].active) unitheap.update[u]++;
      }
      if(graph[u].outdegree>1)
      for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
        int w=outedge[j];
        if(unitheap.update[w]==0){
          unitheap.IncrementKey(w);
          //if(unitheap.LinkedList[w].active) unitheap.IncrementKey(w);
        } else {
          unitheap.update[w]++;
          //if(unitheap.LinkedList[w].active) unitheap.update[w]++;
        }
      }
    }
  }
  if(graph[tmpindex].outdegree<=hugevertex){
    for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
      int w=outedge[i];
      if(unitheap.update[w]==0){
        unitheap.IncrementKey(w);
        //if(unitheap.LinkedList[w].active) unitheap.IncrementKey(w);
      }else{
        unitheap.update[w]++;
        //if(unitheap.LinkedList[w].active) unitheap.update[w]++;
      }
    }
  }
	while(count<finish_num){
    
		int v=unitheap.ExtractMax();
		count++;
		order.push_back(v);
		unitheap.update[v]=INT_MAX/2;

		int popv;
    if(count-window>=0)
      popv=order[count-window];
    else
      popv=-1;

		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					unitheap.update[w]--;
					//if(unitheap.LinkedList[w].active) unitheap.update[w]--;
				}
			}
			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
				if(graph[u].outdegree<=hugevertex){
					unitheap.update[u]--;
					//if(unitheap.LinkedList[u].active) unitheap.update[u]--;
					if(graph[u].outdegree>1)
					if(binary_search(outedge.data() + graph[u].outstart, outedge.data() + graph[u+1].outstart, v)==false){
						for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
							int w=outedge[j];
							unitheap.update[w]--;
							//if(unitheap.LinkedList[w].active) unitheap.update[w]--;
						}
					} else {
						popvexist[u]=true;
					}
				}
			}
		}
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(unlikely(unitheap.update[w]==0)){
					unitheap.IncrementKey(w);
					//if(unitheap.LinkedList[w].active) unitheap.IncrementKey(w);
				} else {
					unitheap.update[w]++;
					//if(unitheap.LinkedList[w].active) unitheap.update[w]++;
				}
				
			}
		}
		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
			if(graph[u].outdegree<=hugevertex){
				if(unlikely(unitheap.update[u]==0)){
					unitheap.IncrementKey(u);
					//if(unitheap.LinkedList[u].active) unitheap.IncrementKey(u);
				} else {
					unitheap.update[u]++;
					//if(unitheap.LinkedList[u].active) unitheap.update[u]++;
				}
				if(popvexist[u]==false){
					if(graph[u].outdegree>1)
					for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
						int w=outedge[j];
						if(unlikely(unitheap.update[w]==0)){
							unitheap.IncrementKey(w);
							//if(unitheap.LinkedList[w].active) unitheap.IncrementKey(w);
						}else{
							unitheap.update[w]++;
							//if(unitheap.LinkedList[w].active) unitheap.update[w]++;
						}
					}
				} else {
					popvexist[u]=false;
				}
			}
		}
	}
 //clock_t hage = clock();
 //cout << "Loop: " << (double)(hage - hoge)/CLOCKS_PER_SEC << endl;
}
*/

void Graph::ReRCMOrder(vector<int>& reorder) {
  vector<int> order;
  RCMOrder(order);
  if(order.size()!=vsize){
    cout << "order.size()!=vsize" << endl;
    quit();
  }

  for (size_t i = 0; i < vsize; i++) {
    reorder[order[i]] = i;
  }
}
void Graph::RCMOrder(vector<int>& retorder){
	queue<int> que;
	bool* BFSflag=new bool[vsize];
	bool* QueFlag=new bool[vsize];
	memset(BFSflag, 0, sizeof(bool)*vsize);
	memset(QueFlag, 0, sizeof(bool)*vsize);

	vector<int> tmp;
	vector<int> degreevertex(vsize);
	for(int i=0; i<vsize; i++){
		degreevertex[i]=i;
	}

	sort(degreevertex.begin(), degreevertex.end(), [&](const int& a, const int& b)->bool{
		if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
			return true;
		else
			return false;
	});

        int now;
	vector<int> order;

	for(int k=0; k<vsize; k++){
		int i=degreevertex[k];
		if(BFSflag[i]==false){
			que.push(i);
//			QueFlag[i]=true;
			BFSflag[i]=true;
			order.push_back(i);

			while(que.empty()==false){
				now=que.front();
				que.pop();

//				BFSflag[now]=true;
				tmp.clear();
				for(int it=graph[now].outstart, limit=graph[now+1].outstart; it<limit; it++){
					tmp.push_back(outedge[it]);
				}
				sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b)->bool{
					if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
						return true;
					else
						return false;
				});
				if(tmp.size()!=graph[now].outdegree)
					cout << "tmp.size()!=graph[now].outdegree" << endl;

				for(int i=0; i<tmp.size(); i++){
//					if((BFSflag[tmp[i]]==false)&&(QueFlag[tmp[i]]==false)){
					if(BFSflag[tmp[i]]==false){
						que.push(tmp[i]);
                        			BFSflag[tmp[i]]=true;
						order.push_back(tmp[i]);
                        		}
				}
        		}
		}
	}

        delete[] BFSflag;
        delete[] QueFlag;

	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}

	retorder.resize(vsize);
	for(int i=0; i<order.size(); i++){
		retorder[order[i]]=order.size()-1-i;
	}
}

unsigned long long Graph::LocalityScore(const int w){
	unsigned long long sum=0;
	for(int i=0; i<vsize; i++){
		for(int j=i-1; j>=i-w && j>=0; j--){
			sum+=IntersectionSize(inedge.data()+graph[i].instart, inedge.data()+graph[j].instart, graph[i].indegree, graph[j].indegree, -1);
			if(binary_search(inedge.data()+graph[i].instart, inedge.data()+graph[i].instart+graph[i].indegree, j))
				sum++;
			if(binary_search(inedge.data()+graph[j].instart, inedge.data()+graph[j].instart+graph[j].indegree, i))
				sum++;
		}
	}
	return sum;
}
}

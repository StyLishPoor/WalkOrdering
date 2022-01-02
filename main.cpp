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
#include <algorithm>
#include <mpi.h>

#include "Graph.h"
#include "Util.h"

using namespace std;
using namespace Gorder;

const int INPUTNUM=1;
const int MSIZE = 10;

int main(int argc, char* argv[]){
	ios::sync_with_stdio(false);
	int i;
	int W=5;
	clock_t start, end, gorder_start, gorder_end, gorder_total, rw_start, rw_end, rw_total;
  gorder_total = 0;
  rw_total = 0;
	string filename;



	if(argc==1){
		cout << "please provide parameter" << endl;
		quit();
	}

  filename=argv[1];
  //cout << filename << endl;
	srand(time(0));

	Graph g;
	string name;
	name=extractFilename(filename.c_str());
	g.setFilename(name);
	start=clock();
	g.readGraph(filename);
	g.Transform();
	//cout << name << " readGraph is complete." << endl;
	end=clock();
	//cout << "Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
  int max=-1;
  int start_node = -1;
  for (size_t i = 0; i < g.vsize; i++) {
    if(g.graph[i].indegree>max) {
      max=g.graph[i].indegree;
      start_node  = i;
    }
  }

  // MPI Init
  int pid, nprocs;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  vector< pair<int, int> > path;
  set<int> visited;
  int sample_size = g.vsize * stof(argv[2]);
  int candidate_size = sample_size / stoi(argv[3]);
  int execute_num = stoi(argv[4]);
  vector<int> candidate;
  vector<int> recv(candidate_size);
  vector<int> order, seqseq;
  vector<double> gapvector, timevector;
  order.reserve(g.vsize);
  random_device seed_gen;
  mt19937 engine(seed_gen());
  uniform_real_distribution<> dist(0.0, 1.0);

  //cout << "PID " << pid << endl;
  int current_node, next_node, init_node;

  // pid 0's job
  while(execute_num-- != 0){
  visited.clear();
  order.clear();
  order.reserve(g.vsize);
  start=clock();
  if (pid==0) {
    // initialize
    candidate.push_back(start_node);
    //seqseq.push_back(start_node);
    visited.insert(start_node);
    init_node = start_node;
    current_node = start_node;
    int check = 0;

    // RW
    while(check < stoi(argv[3])) {
      while(candidate.size() != candidate_size) {
        next_node = g.outedge[g.graph[current_node].outstart + (engine() % g.graph[current_node].outdegree)];
        if (dist(engine) <= static_cast<double>(g.graph[current_node].outdegree) / g.graph[next_node].outdegree) {
          path.push_back(make_pair(current_node, next_node));
          if (visited.count(next_node) == 0) {
            candidate.push_back(next_node);
            visited.insert(next_node);
            //seqseq.push_back(next_node);
          }
          if(g.graph[next_node].outdegree == 0) {
            current_node = start_node;
          } else {
            current_node = next_node;
          }
        } else {
          continue;
        }
      }
     check++;
      MPI_Send(&candidate[0], candidate_size, MPI_INT, check, 0, MPI_COMM_WORLD);
      candidate.clear();
    }
    for (size_t i = 1; i < nprocs; i++) {
      // Receive reordered vector from proc nproc
      MPI_Recv(&recv[0], candidate_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      //create true order
      for (const int v : recv) {
        order.push_back(v);
      }
    }
  } else { // pid > 0 : Gorder process
      MPI_Recv(&recv[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      vector<int> partial_order;
      partial_order.reserve(candidate_size);
      g.GorderGreedy(partial_order, W, recv);
      MPI_Send(&partial_order[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD);

  }

  //MPI_Finalize(); // MPI Finish

  /*
  double rw_ave, gorder_ave;
  rw_ave = (double)(rw_total/check);
  gorder_ave = (double)(gorder_total/check);
  double execute_time;
  cout << "---[Parallel]---" << endl;
  cout << "Execute Time: " << (execute_time)/CLOCKS_PER_SEC << endl;
  cout << "Ideal Time: " << (double)(rw_total + gorder_ave)/CLOCKS_PER_SEC << endl;
  cout << "Size: " << order.size() << endl;
  */
  if (pid==0) {
    //vector<int> retorder;
    end=clock();
    //cout << "Execute Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;;
    timevector.push_back((double)(end-start)/CLOCKS_PER_SEC);
    vector<int> retorder(g.vsize);
    int new_id = 0;
    //cout << "Size: " << order.size() << endl;
    for(const auto v : order) {
      retorder[v] = new_id++;
    }
    gapvector.push_back(g.GapCost(retorder));
  }
  //ofstream original_ofs(argv[4]);
  //ofstream gorder_ofs(argv[5]);
  ///*
  //for (size_t i = 1; i < path.size(); i++) {
  //  original_ofs << path[i-1] << " " << path[i] << endl;
  //  gorder_ofs << retorder[path[i-1]] << " " << retorder[path[i]] << endl;
  //}
  //*/
  //for (const auto edge : path) {
  //  original_ofs << edge.first << " " << edge.second << endl;
  //  gorder_ofs << retorder[edge.first] << " " << retorder[edge.second] << endl;
  //}
  //original_ofs.close();
  //gorder_ofs.close();
  }
  if (pid==0) {
    string tmp(argv[3]);
    string output_file=tmp + "-075.ans";
    double gap_average=0;
    double time_average=0;
    for (size_t i = 0; i < stoi(argv[4]); i++) {
      gap_average += gapvector[i];
      time_average += timevector[i];
    }
    gap_average = gap_average/gapvector.size();
    time_average = time_average/timevector.size();
    ofstream ans(output_file);
    ans << gap_average << " " << time_average;
  }
  /*
  if (pid==0) {
    cout << "---[Perfect]---" << endl;
    vector<int> true_order;
    vector<int> test;
    for (const int i : visited) {
      test.push_back(i);
    }
    true_order.reserve(g.vsize);
    start = clock();
    //g.GorderGreedy(true_order, W, visited);
    g.GorderGreedy(true_order, W, test);
    end = clock();
    cout << "Size: " << true_order.size() << endl;
    cout << "True Time: " << (double)(rw_total + (end-start))/CLOCKS_PER_SEC << endl;

    vector<int> true_retorder(g.vsize);
    int true_new_id = 0;
    for (const auto v : true_order) {
      true_retorder[v] = true_new_id++;
    }
    g.GapCost(true_retorder);
  }
  */
  MPI_Finalize();

  //ofstream gorder(argv[6]);
  //for (const auto edge : path) {
  //  gorder << true_retorder[edge.first] << " " << true_retorder[edge.second] << endl;
  //}
  ///*
  //for (size_t i = 1; i < path.size(); i++) {
  //  gorder << true_retorder[path[i-1]] << " " << true_retorder[path[i]] << endl;
  //}
  //*/
  //*/
  ////g.GapCost(seq_retorder);
  /*
  cout << "---[Seq]---" << endl;
  vector<int> seq_retorder(g.vsize);
  new_id = 0;
  cout << "Size: " << seqseq.size() << endl;
  for (const auto v : seqseq) {
    seq_retorder[v] = new_id++;
  }
  g.GapCost(seq_retorder);
  cout << "---[random]---" << endl;
  //g.GapCost(random_retorder);
  shuffle(seqseq.begin(), seqseq.end(), engine);
  cout << seqseq.size() << endl;
  vector<int> random_retorder(g.vsize);
  new_id = 0;
  cout << "Size: " << seqseq.size() << endl;
  for (const auto v : seqseq) {
    random_retorder[v] = new_id++;
  }
  g.GapCost(random_retorder);
  }
  //ofstream seq_order(argv[7]);
  //for (const auto edge : path) {
  //  seq_order << seq_retorder[edge.first] << " " << seq_retorder[edge.second] << endl;
  //}
  */
  ///*
  //for (size_t i = 1; i < path.size(); i++) {
  //  seq_order << seq_retorder[path[i-1]] << " " << seq_retorder[path[i]] << endl;
  //}

}

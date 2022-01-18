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

int main(int argc, char* argv[]){
	ios::sync_with_stdio(false);
	int i;
	int W=5;
	clock_t start, end, gorder_start, gorder_end, gorder_total, rw_start, rw_end, rw_total;
  gorder_total = 0;
  //rw_total = 0;
	string filename;

	if(argc==1){
		cout << "please provide parameter" << endl;
		quit();
	}

  filename=argv[1];
	srand(time(0));
	Graph g;
	string name;
	name=extractFilename(filename.c_str());
	g.setFilename(name);
	g.readGraph(filename);
  vector<int> original_order(g.vsize);
	original_order = g.Transform();
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

  bool rw_time_flag = true;
  float p = stof(argv[4]);
  int M = stoi(argv[3]);
  vector< pair<int, int> > path;
  set<int> visited;
  int sample_size = g.vsize * stof(argv[2]);
  int candidate_size = sample_size / M;
  int execute_num = stoi(argv[5]);
  vector<int> candidate, collected;
  vector<int> recv(candidate_size);
  vector<int> order;
  //vector<int> seqseq;
  vector<double> gapvector, timevector;
  order.reserve(g.vsize);
  //retorder.reserve(g.vsize);
  random_device seed_gen;
  mt19937 engine(seed_gen());
  uniform_real_distribution<> dist(0.0, 1.0);

  //cout << "PID " << pid << endl;
  int current_node, next_node, init_node;

  // pid 0's job
  while(execute_num-- != 0){
    vector<int> retorder(g.vsize, -1);
    visited.clear();
    order.clear();
    order.reserve(g.vsize);
    start=clock();
    if (pid==0) {
      if(rw_time_flag) {
        set<int> rw;
        rw.insert(start_node);
        int rw_current_node, rw_next_node;
        rw_current_node = start_node;
        // RW time eval
        rw_start = clock();
        while(rw.size() != sample_size) {
          rw_next_node = g.outedge[g.graph[rw_current_node].outstart + (engine() % g.graph[rw_current_node].outdegree)];
          if (dist(engine) <= static_cast<double>(g.graph[rw_current_node].outdegree) / g.graph[rw_next_node].outdegree) {
            if (rw.count(rw_next_node) == 0) {
              rw.insert(rw_next_node);
            }
            if(g.graph[rw_next_node].outdegree == 0) {
              rw_current_node = start_node;
            } else {
              rw_current_node = rw_next_node;
            }
          } else {
            continue;
          }
        }
        rw_end = clock();
        ofstream rw_time("rw-time.txt", std::ios::app);
        rw_time << (double)(rw_end - rw_start)/CLOCKS_PER_SEC << endl;
        rw_time_flag = false;
      }
      // initialize
      rw_total = 0;
      candidate.push_back(start_node);
      //seqseq.push_back(start_node);
      visited.insert(start_node);
      init_node = start_node;
      current_node = start_node;
      int check = 0;

      // RW
      while(check < M) {
        while(candidate.size() != candidate_size) {
          next_node = g.outedge[g.graph[current_node].outstart + (engine() % g.graph[current_node].outdegree)];
          if (dist(engine) <= static_cast<double>(g.graph[current_node].outdegree) / g.graph[next_node].outdegree) {
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
        //vector<int> collected(visited.begin(), visited.end());
        //MPI_Send(&collected[0], candidate_size * check, MPI_INT, check, 0, MPI_COMM_WORLD);
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
        //int collected_size = candidate_size * pid;
        MPI_Recv(&recv[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        //vector<int> recv_collected(collected_size);
        //MPI_Recv(&recv_collected[0], collected_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        Graph sub_g;
        //g.SubGraph(sub_g, recv);
        //g.SubGraphTest(sub_g, recv_collected);
        g.SubGraphTest(sub_g,recv,p);
        //g.SubGraphTest2(sub_g,recv,(double)stof(argv[4]));
        vector<int> partial_order;
        partial_order.reserve(candidate_size);
        //g.NDGorderGreedy(partial_order, W, recv);
        sub_g.GorderSubGreedy(partial_order, W, recv);
        //g.GorderGreedy(debug_order, W, recv);
        //set<int> pset(partial_order.begin(), partial_order.end());
        //set<int> recv_test(recv.begin(), recv.end());
        //set<int> dset(debug_order.begin(), debug_order.end());
        //set<int> result;
        //set_intersection(pset.begin(), pset.end(), recv_test.begin(), recv_test.end(), inserter(result, result.end()));
        //cout << "-----" << endl;
        //cout << pset.size() << " " << recv.size() << " " << result.size() << endl;
        //cout << "-----" << endl;
        //g.GorderTestSubGreedy(partial_order, W, recv, recv_collected);
        //g.GorderGreedy(partial_order, W, recv);
        //cout << "TEST: " << partial_order.size() << endl;
        MPI_Send(&partial_order[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD);

    }
    if (pid==0) {
      //vector<int> retorder;
      end=clock();
      //cout << "Execute Time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;;
      timevector.push_back((double)(end-start)/CLOCKS_PER_SEC);
      //vector<int> retorder(g.vsize);
      int new_id = 0;
      for(const auto v : order) {
        retorder[v] = new_id++;
      }
      gapvector.push_back(g.GapCostV(retorder, visited));
      string tmp1(argv[3]);
      string tmp2 = to_string(execute_num+1);
      string output_graph = tmp1 + "-" + tmp2 + "-sample.txt";
      ofstream outgraph(output_graph);
      g.WriteSampleGraph(visited, retorder, outgraph);
      retorder.clear();
    }
  } // end average loop

  if (pid==0) {
    string tmp(argv[3]);
    string output_file=tmp + "-.ans";
    //string output_graph = tmp + "-sample.txt";
    //string output_rwtime = "rw-time.txt";
    double gap_average=0;
    double time_average=0;
    for (size_t i = 0; i < stoi(argv[5]); i++) {
      gap_average += gapvector[i];
      time_average += timevector[i];
    }
    gap_average = gap_average/gapvector.size();
    time_average = time_average/timevector.size();
    //cout << "Gap: " << gap_average << " Time: " << time_average << endl;
    ofstream ans(output_file);
    ans << gap_average << " " << time_average;
    //ofstream outgraph(output_graph);
    //g.WriteSampleGraph(visited, retorder, outgraph);
    //ofstream rwtime(output_rwtime);
    //cout << rw_total << endl;
    //rwtime << rw_total/CLOCKS_PER_SEC;
    /*
    for (const auto edge : path) {
      outgraph << retorder[edge.first] << " " << retorder[edge.second] << endl;
    }
    */
  }
  if (pid==0 && M==1) {
    //vector<int> seq_retorder(g.vsize, -1);
    ////seq_retorder.reserve(g.vsize);
    //int new_id = 0;
    //for (const auto v : seqseq) {
    //  seq_retorder[v] = new_id++;
    //}
    ////cout << "SEQ GAP" << endl;
    //g.GapCostV(seq_retorder, visited);
    ////shuffle(seqseq.begin(), seqseq.end(), engine);
    //vector<int> random_retorder(g.vsize, -1);
    ////random_retorder.reserve(g.vsize);
    //new_id = 0;
    //for (const auto v : seqseq) {
    //  random_retorder[v] = new_id++;
    //}
    ////cout << "RANDOM GAP" << endl;
    //g.GapCostV(random_retorder, visited);
    //vector<int> reorder;
    //reorder.reserve(g.vsize);
    //g.ReRCMOrder(reorder);

    ofstream original_graph("original.txt");
    g.WriteSampleGraph(visited, original_order, original_graph);
    ofstream random_graph("random.txt");
    g.WriteSampleRandomGraph(visited, random_graph);
    /*
    for (const auto edge : path) {
      original_graph << edge.first << " " << edge.second << endl;
    }
    */
    //ofstream seq_graph("seq.txt");
    //g.WriteSampleGraph(visited, seq_retorder, seq_graph);
    ///*
    //for (const auto edge : path) {
    //  seq_graph << seq_retorder[edge.first] << " " << seq_retorder[edge.second] << endl;
    //}
    //*/
    //ofstream random_graph("random.txt");
    //g.WriteSampleGraph(visited, random_retorder, random_graph);
    /*
    for (const auto edge : path) {
      random_graph << random_retorder[edge.first] << " " << random_retorder[edge.second] << endl;
    }
    */
  }
  MPI_Finalize();
}

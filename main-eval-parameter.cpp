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

  // MPI Init
  int pid, nprocs;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  filename=argv[1];
	srand(time(0));
	Graph g;
	string name;
	name=extractFilename(filename.c_str());
	g.setFilename(name);
  //clock_t tmp_start = clock();
	g.readGraph(filename);
  //clock_t tmp_end = clock();
  //cout << "Read Time: " << (double)(tmp_end - tmp_start)/CLOCKS_PER_SEC << endl;
  vector<int> original_order;
	if (pid !=0) {
    g.Transform();
  } else {
    original_order = g.Transform();
  }

  int max=-1;
  int start_node = -1;
  for (size_t i = 0; i < g.vsize; i++) {
    if(g.graph[i].indegree>max) {
      max=g.graph[i].indegree;
      start_node  = i;
    }
  }

  bool rw_time_flag = true;
  float p = stof(argv[3]);
  int M = 5;
  //vector< pair<int, int> > path;
  int sample_size = g.vnum * stof(argv[2]);
  int candidate_size = sample_size / M;
  int execute_num = stoi(argv[4]);
  vector<int> recv(candidate_size);
  vector<double> gapvector, timevector;
  //order.reserve(g.vsize);
  random_device seed_gen;
  mt19937 engine(seed_gen());
  uniform_real_distribution<> dist(0.0, 1.0);

  //cout << "PID " << pid << endl;
  int current_node, next_node, init_node;

  // pid 0's job
  while(execute_num-- != 0){
    if (pid==0) {
      vector<int> candidate;
      set<int> visited;
      vector<int> order;
      order.reserve(g.vsize);
      vector<int> retorder(g.vsize, -1);
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
        string graphname(argv[1]);
        int split_start=graphname.find_last_of('/');
        string after_graphname = graphname.substr(split_start+1, graphname.size()-split_start-5);
        string tmp(argv[3]);
        string output_file=after_graphname+".rwtime";
        ofstream rw_time(output_file, std::ios::app);
        rw_time << (double)(rw_end - rw_start)/CLOCKS_PER_SEC << endl;
        rw_time_flag = false;
      }
      // initialize
      start=clock();
      candidate.push_back(start_node);
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
        //int collect_size = visited.size();
        //MPI_Send(&collected[0], collect_size, MPI_INT, check, 0, MPI_COMM_WORLD);
        candidate.clear();
      }
      // RW finish
      for (size_t i = 1; i < nprocs; i++) {
        // Receive reordered vector from proc nproc
        MPI_Recv(&recv[0], candidate_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        //create true order
        for (const int v : recv) {
          order.push_back(v);
        }
      }
      end=clock();
      timevector.push_back((double)(end-start)/CLOCKS_PER_SEC);
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
      //retorder.clear();
      if (execute_num == 0) {
        string graphname(argv[1]);
        int split_start=graphname.find_last_of('/');
        string after_graphname = graphname.substr(split_start+1, graphname.size()-split_start-5);
        string tmp(argv[3]);
        string output_file=after_graphname+"-"+tmp + ".ans";
        double gap_average=0;
        double time_average=0;
        sort(gapvector.begin(), gapvector.end());
        sort(timevector.begin(), timevector.end());
        int out_range = stoi(argv[4]) * 0.1;
        for (size_t i = out_range; i < stoi(argv[4])-out_range; i++) {
          gap_average += gapvector[i];
          time_average += timevector[i];
        }
        gap_average = gap_average/gapvector.size();
        time_average = time_average/timevector.size();
        ofstream ans(output_file);
        ans << gap_average << " " << time_average;

        if (p==0) {
          ofstream original_graph("original.txt");
          g.WriteSampleGraph(visited, original_order, original_graph);
          vector<int> random_retorder(g.vsize);
          for (size_t i = 0; i <= g.vsize; i++) {
            random_retorder[i] = i;
          }
          shuffle(random_retorder.begin(), random_retorder.end(), engine);
          ofstream random_graph("random.txt");
          g.WriteSampleGraph(visited, random_retorder, random_graph);
        }
      }
    } else { // pid > 0 : Gorder process
        MPI_Recv(&recv[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        //vector<int> recv_collected(candidate_size * pid);
        //MPI_Recv(&recv_collected[0], candidate_size * pid, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        Graph sub_g;
        //g.SubGraphTest(sub_g,recv,p);
        g.SubGraphTest2(sub_g, recv);
        //g.SubGraphTest2(sub_g,recv_collected);
        vector<int> partial_order;
        partial_order.reserve(candidate_size);
        //sub_g.GorderSubGreedy(partial_order, W, recv);
        //if (pid==5) { p = 0; }
        sub_g.GorderSubGreedy2(partial_order, W, recv, p);
        MPI_Send(&partial_order[0], candidate_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        //cout << "Finish: " << pid << endl;
    }
  } // end average loop

  // after process
  //if (pid==0) {
  //  string tmp(argv[3]);
  //  string output_file=tmp + "-.ans";
  //  double gap_average=0;
  //  double time_average=0;
  //  for (size_t i = 0; i < stoi(argv[5]); i++) {
  //    gap_average += gapvector[i];
  //    time_average += timevector[i];
  //  }
  //  cout << gap_average/gapvector.size() << endl;
  //  gap_average = gap_average/gapvector.size();
  //  time_average = time_average/timevector.size();
  //  ofstream ans(output_file);
  //  ans << gap_average << " " << time_average;

  //  if (M==1) {
  //    ofstream original_graph("original.txt");
  //    g.WriteSampleGraph(visited, original_order, original_graph);
  //    vector<int> random_retorder(g.vsize);
  //    for (size_t i = 0; i <= g.vsize; i++) {
  //      random_retorder[i] = i;
  //    }
  //    shuffle(random_retorder.begin(), random_retorder.end(), engine);
  //    ofstream random_graph("random.txt");
  //    g.WriteSampleGraph(visited, random_retorder, random_graph);
  //  }
  //}
  MPI_Finalize();
}

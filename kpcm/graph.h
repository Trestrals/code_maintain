#pragma once
#ifndef GRAPH_H
#define GRAPH_H
#include "utility.h"
class Graph {
public:

	/**
		p_t : std::pair<long long int, long long int>
		p_t.first : subgraph degree.
		p_t.second : degree of the whole graph.
	*/

	int num_node;
	int num_edge;
	vid_t max_id;
	std::vector<std::vector<vid_t>> neighbor; // neighbors of a node
	std::vector<int> degree; // node degree
	std::vector<std::vector<p_t>> node_kpcore_info; 
		// node_kpcore_info[a][b] : a's p value in b + 1 core.


	// k-order
	std::vector<vid_t> k_order;
	std::vector<int> k_rank;
	std::vector<int> k_deg_plus;
	std::vector<int> k_deg_star;
	std::vector<int> k_mcd;

	std::vector<std::vector<std::pair<p_t, std::vector<vid_t>>>> naive_index;

	//new neighbor
	std::vector<std::vector<vid_t>> order_neighbor; // neighbors of a node
	std::vector<std::vector<int>> begin_index; 

	// supporter part
	std::vector<std::vector<int>> not_supporter;
	std::vector<std::vector<int>> maybe_supporter;
	std::vector<std::vector<int>> supporter;

	//pOrder
	std::vector<std::vector<int>> rank; //rank[a][b] : vid = a, k = b;
	std::vector<std::vector<int>> new_order;
	bool* dirty;
	int* change_tag;
	std::vector<std::vector<std::pair<p_t, vid_t>>> p_order;
	std::vector<int> deg_star;
	// p_order for k-core : p_order[k]
	// rank from p_order[k].size()-1 to 0 (p from small to large)
	//end pOrder

	//affect subgraph
	int* include_in_subgraph;
	vid_t* queue_for_aff_graph;
	int maintain_cnt;
	//end affect subgraph


	Graph();
	vid_t get_node_id(int id);
	void insert_edge(vid_t u, vid_t v);
	void remove_edge(vid_t u, vid_t v);
	void loadGraph(std::string dir);
	double kpcore_decom();
	double global_insert_edge(vid_t u, vid_t v);
	double global_delete_edge(vid_t u, vid_t v);
	bool isEdge(vid_t u, vid_t v);
	void generate_random_edge(std::string stream_location);
	bool compare_kpcore(Graph* g);
	void local_insert_maintain_p(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);

	bool check_decrease(int k, int u, int new_link_to);

	long long get_decrease_set_vertex_statistics(int k, p_t p, std::unordered_set<vid_t> &decrease_set);
	std::vector<long long> pIncrease_vertex_statistics(int k, int begin_rank);
	std::vector<long long> local_insert_edge_vertex_statistics(vid_t u, vid_t v);
	std::vector<long long> local_insert_maintain_p_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);
	std::vector<long long> pDecrease_multiple_vertices_vertex_statistics(int k, std::vector<vid_t> othernodes);


	void pIncrease(int k, int begin_rank);
	void pDecrease_multiple_vertices(int k, std::vector<vid_t> othernodes);
	void get_decrease_set(int k, p_t p, std::unordered_set<vid_t> &decrease_set);

	void local_delete_maintain_p_for_changed_core_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);

	void print_new_order(int k);
	void print_rank(int k);

	p_t get_p_with_default(vid_t u, int k);

	std::vector<vid_t> get_bfs_order();
	void print_subgraph(std::vector<vid_t> nodes, std::string dir);


	void local_insert_edge(vid_t u, vid_t v);
	void local_insert_edge_woc(vid_t u, vid_t v);
	std::vector<double> local_insert_edge_test(vid_t u, vid_t v);
	long long global_insert_edge_vertex_statistics(vid_t u, vid_t v);

	void local_delete_edge(vid_t u, vid_t v);
	void local_delete_edge_woc(vid_t u, vid_t v);

	void local_delete_maintain_p(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v);
	void local_delete_maintain_p_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v);
	void local_delete_maintain_p_for_changed_core(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);
	void local_insert_maintain_p_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);
	bool check_increase_delete(int k, int u);

	std::vector<long long> compare_p_change_inc_dec(Graph* g, vid_t u, vid_t v);
	std::vector<long long> local_delete_edge_vertex_statistics(vid_t u, vid_t v);
	std::vector<long long> local_delete_maintain_p_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v);
	std::vector<long long> local_delete_maintain_p_for_changed_core_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes);
	long long global_delete_edge_vertex_statistics(vid_t u, vid_t v);
	void print_k_order();

	void pIncrease_woc(int k, int begin_rank);

	void find_smallest_k_p_set(int k, p_t p, int target);

};

#endif // !GRAPH_H

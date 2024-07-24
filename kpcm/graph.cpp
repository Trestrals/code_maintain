#include "graph.h"

using namespace std;


Graph::Graph() {
	num_node = 0;
	num_edge = 0;
	max_id = 0;
	maintain_cnt = 0;
}

vid_t Graph::get_node_id(int id) {
	return id;
}

void Graph::insert_edge(vid_t u, vid_t v) {
	neighbor[u].push_back(v);
	neighbor[v].push_back(u);
	degree[u]++;
	degree[v]++;
	num_edge++;
}

void Graph::remove_edge(vid_t u, vid_t v) {
	for (int i = 0; i < degree[u]; i++) {
		int vv = neighbor[u][i];
		if (vv == v) {
			swap(neighbor[u][i], neighbor[u][degree[u] - 1]);
			--degree[u];
			neighbor[u].pop_back();
			num_edge--; //only once!!!
			break;
		}
	}
	for (int i = 0; i < degree[v]; i++) {
		int vv = neighbor[v][i];
		if (vv == u) {
			swap(neighbor[v][i], neighbor[v][degree[v] - 1]);
			--degree[v];
			neighbor[v].pop_back();
			break;
		}
	}

}

void Graph::loadGraph(string dir) {
	ifstream input;
	input.open(dir);
	int u, v;
	input >> u >> v;
	cout << "num node : " << u << endl;
	cout << "num edge : " << v << endl;
	num_node = u; //num_edge = v;
	degree.resize(num_node); neighbor.resize(num_node);
	dirty = new bool[num_node];
	change_tag = new int[num_node];
	for(int i = 0; i < num_node; i++){
		dirty[i] = false;
		change_tag[i] = 0;
	}

	//add for affect graph
	include_in_subgraph = new int[num_node * 2];
	queue_for_aff_graph = new vid_t[num_node * 2];
	//end

	fill_n(degree.begin(), degree.size(), 0);
	while (input >> u >> v) {
		vid_t x = get_node_id(u);
		vid_t y = get_node_id(v);
		insert_edge(x, y);
	}
	input.close();
}

double Graph::kpcore_decom() {
	auto deg_partial = degree; // degree for partial graph? (vector<int>)
	// kpcore_decom_result.clear();
	naive_index.clear();
	node_kpcore_info.clear();
	
	//add
	rank.clear();
	// p_order.clear();
	new_order.clear();
	// change_tag.clear();

	k_order.clear();
	k_rank.clear();
	k_deg_plus.clear();
	deg_star.clear();
	k_mcd.clear();

	//end

	vector<bool> is_removed;
	is_removed.resize(num_node);
	fill_n(is_removed.begin(), is_removed.size(), false);
	int max_deg = 0;
	int min_deg = degree[0];
	for (int i = 0; i < num_node; i++) {
		max_deg = max_deg < degree[i] ? degree[i] : max_deg;
		min_deg = min_deg < degree[i] ? min_deg : degree[i];
	}
	// attention! mindeg is not 0, since there may be some same core!
	// printf("min_deg: %d\n", min_deg);
	// in Facebook, min_deg = 1, seems no use.
	naive_index.resize(max_deg + 1); //naive_index ?

	//add
	new_order.resize(max_deg + 1);
	// change_tag.resize(num_node + 1);
	rank.resize(num_node + 1);
	//for(auto x : rank) x.push_back(0);
	for(int i = 0; i <= num_node; i++) {
		rank[i].push_back(0);
		// change_tag[i] = 0;
	}
	int subgraph_num_node = num_node;

	k_rank.resize(num_node + 1); 
	k_deg_plus.resize(num_node + 1);
	k_deg_star.resize(num_node + 1);
	deg_star.resize(num_node + 1);
	k_mcd.resize(num_node + 1);
	
	//end

	int pre = 0;
	auto start = chrono::system_clock::now();
	for (int k = min_deg; k != MAX_DEG; k = min_deg) {
		// printf("BIGK - %d\n", k);
		// process pcore decomposition for some fixed k
		for (k = pre + 1; k <= min_deg; k++) {
			// printf("SMALL K - %d\n", k);
			//add
			new_order[k].resize(subgraph_num_node);
			int now_rank = subgraph_num_node - 1;//rank 从后往前排，增加节点时加在最后符合操作逻辑
			//end

			auto deg_in_pcore = deg_partial;
			auto is_removed_in_pcore = is_removed;
			// create a min heap
			min_heap mheap;
			mheap.node2index.clear();
			mheap.node2index.resize(num_node, -1);

			for (vid_t i = 0; i < num_node; i++) {
				if (is_removed_in_pcore[i]) continue;
				p_t t = make_pair(deg_in_pcore[i], degree[i]);
				mheap.insert({ t,i });
			}

			std::queue<std::pair<p_t, int> > not_satisfy_k;
			
			bool notSatisfy = false;
			p_t now_p = make_pair(0, 1);
			std::vector<vid_t> naive_index_vec;
			naive_index_vec.clear();
			while((!mheap.empty()) || (!not_satisfy_k.empty())) {
				
				std::pair<p_t, vid_t> top_ele;

				if(!not_satisfy_k.empty()) {
					top_ele = not_satisfy_k.front();
					top_ele.first = make_pair(deg_in_pcore[top_ele.second], degree[top_ele.second]);
					not_satisfy_k.pop();
					notSatisfy = true;
				}
				else {
					notSatisfy = false;
					top_ele = mheap.top();
					mheap.pop(); // ?
				}

				int now_id = top_ele.second;
				

				if(compare_pvalue(now_p, top_ele.first) == -1 && notSatisfy == false) {
					if(!naive_index_vec.empty()) {
						naive_index[k].push_back(make_pair(now_p, naive_index_vec)); 
						// is vector a pointer? naive_index_vec will change in the process?
						// clear, it's a copy, not a pointer.
						naive_index_vec.clear();
					}
					now_p = top_ele.first;
				}
				naive_index_vec.push_back(now_id);
				new_order[k][now_rank] = top_ele.second;
				rank[now_id].push_back(now_rank);
				now_rank--;

				is_removed_in_pcore[now_id] = true;

				for(auto x : neighbor[now_id]) {
					if(is_removed_in_pcore[x]) continue;
					else {
						deg_in_pcore[x]--;
						if(deg_in_pcore[x] == k - 1) {
							not_satisfy_k.push(make_pair(make_pair(0, 0), x));
							mheap.remove(mheap.node2index[x]);
						}
						else if(deg_in_pcore[x] >= k){
							mheap.decrease_key(make_pair(deg_in_pcore[x], degree[x]), mheap.node2index[x]);
						}
						
						
					}
				}
				
			}

			if(!naive_index_vec.empty()) naive_index[k].push_back(make_pair(now_p, naive_index_vec)); 

			//print_p_order(k);

		}

		

		k = min_deg;

		
		// peel graph for next k-core
		vector<vid_t> stack;
		for (vid_t i = 0; i < num_node; i++) {
			if (!is_removed[i] && deg_partial[i] <= k) {
				stack.push_back(i);
			}
		}
		while (!stack.empty()) {
			vid_t v = stack.back();
			stack.pop_back();
			if (is_removed[v]) continue;
			is_removed[v] = true;
			
			//add
			subgraph_num_node --;
			k_rank[v] = k_order.size();
			k_order.push_back(v);
			k_deg_plus[v] = deg_partial[v];
			//end

			for (int i = 0; i < neighbor[v].size(); i++) {
				vid_t neigh = neighbor[v][i];
				if (is_removed[neigh]) continue;
				deg_partial[neigh]--;
				if (deg_partial[neigh] == k) {
					stack.push_back(neigh);
				}
			}
		}
		pre = min_deg;
		min_deg = MAX_DEG;
		for (int i = 0; i < num_node; i++) {
			if (is_removed[i]) continue;
			min_deg = min_deg < deg_partial[i] ? min_deg : deg_partial[i];
		}
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	node_kpcore_info.resize(num_node);
	for (auto e : node_kpcore_info) {
		e.push_back(make_pair(0, 0));
	}
	for (int k = 1; k < naive_index.size(); k++) {
		for (int p = 0; p < naive_index[k].size(); p++) {
			for (auto node : naive_index[k][p].second) {
				node_kpcore_info[node].push_back(naive_index[k][p].first);
			}
		}
	}

	for (int i = 0; i < num_node; i++) {
		int temp = 0;
		int ik = node_kpcore_info[i].size();
		for (int j : neighbor[i]) {
			if (node_kpcore_info[j].size() >= ik) {
				temp++;
			}
		}
		k_mcd[i] = temp;
	}

	//print_p_order(2);

	return elapsed_seconds.count();
}

bool cmp_sort_p(p_t a, p_t b) {
	return compare_pvalue(a, b) < 0;
}

void Graph::local_delete_edge_woc(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}

	int pre_k_u = k_u;
	int pre_k_v = k_v;

	int k_min = MIN(k_u, k_v); 
	remove_edge(u, v);

	if (k_u <= k_v) {
		k_mcd[u] --;
	}
	if(k_v <= k_u) {
		k_mcd[v] --;
	}

	bool kcore_change = false;

	// calculate mcd brute force
	unordered_map<vid_t, int> cd;
	vector<vid_t> V_star;
	unordered_set<vid_t> in_V_star;

	int head = 0;
	int tail = 0;

	cd[u] = k_mcd[u];
	cd[v] = k_mcd[v];

	if (cd[u] < k_min) {
		V_star.push_back(u);
		in_V_star.insert(u);
		tail++;
	}
	if (cd[v] < k_min) {
		V_star.push_back(v);
		in_V_star.insert(v);
		tail++;
	}

	if(head < tail) {
		
		kcore_change = true;

		while (head < tail) {
			vid_t now = V_star[head++];
			for (int nei : neighbor[now]) {
				if (in_V_star.find(nei) != in_V_star.end()) {
					continue;
				}
				if (node_kpcore_info[nei].size() == k_min) {
					if (cd.find(nei) != cd.end()) {
						cd[nei] --;
					}
					else {
						cd[nei] = k_mcd[nei] - 1;
					}
					if(cd[nei] < k_min) {
						V_star.push_back(nei);
						in_V_star.insert(nei);
						tail++;
					}
				}
			}
		}

		int num = V_star.size();
		for (int i = 0; i < num; i++) {
			vid_t now = V_star[i];
			k_deg_plus[now] = 0;
			for (vid_t nei : neighbor[now]) {
				if(node_kpcore_info[nei].size() == k_min){
					if(in_V_star.find(nei) == in_V_star.end()) k_mcd[nei] --;
					if (k_rank[nei] < k_rank[now]) k_deg_plus[nei] --;
				} 
				if(node_kpcore_info[nei].size() >= k_min || in_V_star.find(nei) != in_V_star.end()) {
					k_deg_plus[now] ++;
				}
				if(node_kpcore_info[nei].size() == k_min - 1) {
					k_mcd[now]++;
				}
			}
			in_V_star.erase(now);
		}

		int begin_rank = k_rank[V_star[0]];
		int end_rank = k_rank[V_star[0]];

		while(begin_rank >= 0 && node_kpcore_info[k_order[begin_rank]].size() == k_min) begin_rank --;
		while(end_rank < k_order.size() && node_kpcore_info[k_order[begin_rank]].size() == k_min) end_rank ++;
		begin_rank ++;
		end_rank --;

		for (vid_t now : V_star) {
			in_V_star.insert(now);
		}

		int determined_rank = end_rank + 1;
		for (int i = end_rank; i >= begin_rank; i--) {
			if (in_V_star.find(k_order[i]) == in_V_star.end()) {
				k_order[--determined_rank] = k_order[i];
				k_rank[k_order[determined_rank]] = determined_rank;
			}
		}

		for (int i = 0; i < num; i++){
			k_order[begin_rank + i] = V_star[i];
			k_rank[V_star[i]] = begin_rank + i;
		}

		if (k_min == 1) {
			for (vid_t now : V_star) {
				node_kpcore_info[now].pop_back();
			}
		}
	}

	// puts("After k update");
	
	// need to deal with the condition that p_order[min_k]'s size equal to 0
	// after the deletion (there's no node in min_k core)
	
	int max_k =(k_u <  k_v) ? k_v : k_u;

	if(k_rank[v] < k_rank[u]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min && kcore_change) {
			// printf("new k = %d, size = %d\n", k_min, changed_nodes.size()); 
			// print_new_order(k);
			local_delete_maintain_p_for_changed_core_woc(k, u, v, pre_k_u, pre_k_v, V_star);
			// print_new_order(k);
		}
		else {
			// printf("else k = %d\n", k);
			// print_new_order(k);
			local_delete_maintain_p_woc(k, u, v, pre_k_u, pre_k_v); // revise
			// print_new_order(k);
		}
	}

	return;


}

void Graph::local_delete_edge(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}

	int pre_k_u = k_u;
	int pre_k_v = k_v;

	int k_min = MIN(k_u, k_v); 
	remove_edge(u, v);

	if (k_u <= k_v) {
		k_mcd[u] --;
	}
	if(k_v <= k_u) {
		k_mcd[v] --;
	}

	bool kcore_change = false;

	// calculate mcd brute force
	unordered_map<vid_t, int> cd;
	vector<vid_t> V_star;
	unordered_set<vid_t> in_V_star;

	int head = 0;
	int tail = 0;

	cd[u] = k_mcd[u];
	cd[v] = k_mcd[v];

	if (cd[u] < k_min) {
		V_star.push_back(u);
		in_V_star.insert(u);
		tail++;
	}
	if (cd[v] < k_min) {
		V_star.push_back(v);
		in_V_star.insert(v);
		tail++;
	}

	if(head < tail) {
		
		kcore_change = true;

		while (head < tail) {
			vid_t now = V_star[head++];
			for (int nei : neighbor[now]) {
				if (in_V_star.find(nei) != in_V_star.end()) {
					continue;
				}
				if (node_kpcore_info[nei].size() == k_min) {
					if (cd.find(nei) != cd.end()) {
						cd[nei] --;
					}
					else {
						cd[nei] = k_mcd[nei] - 1;
					}
					if(cd[nei] < k_min) {
						V_star.push_back(nei);
						in_V_star.insert(nei);
						tail++;
					}
				}
			}
		}

		int num = V_star.size();
		for (int i = 0; i < num; i++) {
			vid_t now = V_star[i];
			k_deg_plus[now] = 0;
			for (vid_t nei : neighbor[now]) {
				if(node_kpcore_info[nei].size() == k_min){
					if(in_V_star.find(nei) == in_V_star.end()) k_mcd[nei] --;
					if (k_rank[nei] < k_rank[now]) k_deg_plus[nei] --;
				} 
				if(node_kpcore_info[nei].size() >= k_min || in_V_star.find(nei) != in_V_star.end()) {
					k_deg_plus[now] ++;
				}
				if(node_kpcore_info[nei].size() == k_min - 1) {
					k_mcd[now]++;
				}
			}
			in_V_star.erase(now);
		}

		int begin_rank = k_rank[V_star[0]];
		int end_rank = k_rank[V_star[0]];

		while(begin_rank >= 0 && node_kpcore_info[k_order[begin_rank]].size() == k_min) begin_rank --;
		while(end_rank < k_order.size() && node_kpcore_info[k_order[begin_rank]].size() == k_min) end_rank ++;
		begin_rank ++;
		end_rank --;

		for (vid_t now : V_star) {
			in_V_star.insert(now);
		}

		int determined_rank = end_rank + 1;
		for (int i = end_rank; i >= begin_rank; i--) {
			if (in_V_star.find(k_order[i]) == in_V_star.end()) {
				k_order[--determined_rank] = k_order[i];
				k_rank[k_order[determined_rank]] = determined_rank;
			}
		}

		for (int i = 0; i < num; i++){
			k_order[begin_rank + i] = V_star[i];
			k_rank[V_star[i]] = begin_rank + i;
		}

		if (k_min == 1) {
			for (vid_t now : V_star) {
				node_kpcore_info[now].pop_back();
			}
		}
	}

	// puts("After k update");
	
	// need to deal with the condition that p_order[min_k]'s size equal to 0
	// after the deletion (there's no node in min_k core)
	
	int max_k =(k_u <  k_v) ? k_v : k_u;

	if(k_rank[v] < k_rank[u]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min && kcore_change) {
			// printf("new k = %d, size = %d\n", k_min, changed_nodes.size()); 
			// print_new_order(k);
			local_delete_maintain_p_for_changed_core(k, u, v, pre_k_u, pre_k_v, V_star);
			// print_new_order(k);
		}
		else {
			// printf("else k = %d\n", k);
			// print_new_order(k);
			local_delete_maintain_p(k, u, v, pre_k_u, pre_k_v); // revise
			// print_new_order(k);
		}
	}

	return;


}

void Graph::local_insert_edge(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	vector<vid_t> changed_nodes;
	vector<vid_t> new_k_order;

	// printf("before insert : ");
	// print_k_order();

	// add for k order insertion
	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}
	// end

	//inserted clip of my method -----------------
	int pre_k_u = k_u;
	int pre_k_v = k_v;
	//end of the clip ----------------------------

	int k_min = MIN(k_u, k_v); 
	insert_edge(u, v);


	// new method to maintain k value
	k_deg_plus[u] += 1;
	// printf("u %d, k_deg_plus = %d, k_mind = %d\n", u, k_deg_plus[u], k_min);

	if (k_deg_plus[u] <= k_min) {
		// puts("doneee");
		// done
	} else {
		// remember to reset k_deg_star value !
		// puts("IN else !");
		unordered_set<vid_t> vc;
		int now_rank = k_rank[u];
		int max_rank = k_order.size();
		while (now_rank < max_rank) {
			int now_id = k_order[now_rank];
			// printf("node_id : %d (+=%d, *=%d)\n", now_id, k_deg_plus[now_id], k_deg_star[now_id]);
			if (node_kpcore_info[k_order[now_rank]].size() > k_min) {
				break;
			} else {
				if (k_deg_star[now_id] + k_deg_plus[now_id] > k_min) {
					// puts("con1");
					vc.insert(now_id);
					for (vid_t nei : neighbor[now_id]) {
						if (node_kpcore_info[nei].size() == k_min && k_rank[nei] > now_rank) {
							k_deg_star[nei] ++;
						}
					}
					now_rank ++;
				} else if (k_deg_star[now_id] == 0) {
					// puts("con2");
					new_k_order.push_back(now_id);
					// k_rank[now_id] = determined_rank;
					now_rank ++;
				} else {
					// puts("con3");
					new_k_order.push_back(now_id);
					k_deg_plus[now_id] += k_deg_star[now_id];
					k_deg_star[now_id] = 0;

					vector<vid_t> que;
					unordered_set<vid_t> in_que;
					int head = 0;
					int tail = 0;
					for (vid_t nei : neighbor[now_id]) {
						// printf("nei %d, %d-%d\n", nei, k_deg_plus[nei], k_deg_star[nei]);
						if (vc.find(nei) != vc.end()) {
							k_deg_plus[nei] --;
							if (k_deg_plus[nei] + k_deg_star[nei] <= k_min) {
								que.push_back(nei);
								in_que.insert(nei);
								tail++;
							}
						}
					}

					while (head < tail) {
						// dequeue w
						int w = que[head++];
						// printf("%d out\n", w);

						// update k_deg_plus and k_deg_star
						k_deg_plus[w] += k_deg_star[w];
						k_deg_star[w] = 0;

						// delete w from Vc
						vc.erase(w);

						// append w to Ok
						new_k_order.push_back(w);
						// k_rank[w] = determined_rank;

						for (vid_t nei : neighbor[w]) {
							// printf("now's nei %d\n", nei);
							if (k_rank[nei] > k_rank[w]) {
								if(k_deg_star[nei] > 0) {
									k_deg_star[nei] --;
								}
								if (vc.find(nei) != vc.end()) {
									// printf("find! , +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
									// printf("%d\n", in_que.find(nei) == in_que.end() ? 0 : 1);
									if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
										que.push_back(nei);
										in_que.insert(nei);
										++ tail;
									}
								}
							}
							else if(vc.find(nei) != vc.end()) {
								k_deg_plus[nei] --;
								// printf("find2 ! +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
								if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
									que.push_back(nei);
									in_que.insert(nei);
									++ tail;
								}
							}
						}
					}

					now_rank ++;

				}
			}
			
		}

		vector<int> rank_sort;
	
		for (vid_t now : vc) {
			k_deg_star[now] = 0;
			node_kpcore_info[now].push_back({0,1});
			changed_nodes.push_back(now);
			rank[now].push_back(0);
			rank[now][k_min+1] = new_order[k_min + 1].size();
			new_order[k_min + 1].push_back(now);
			rank_sort.push_back(k_rank[now]);
		}

		sort(rank_sort.begin(), rank_sort.end());
		for (int temp : rank_sort) {
			new_k_order.push_back(k_order[temp]);
		}

		int min_rank = k_rank[u];
		for (vid_t now : new_k_order) {
			k_rank[now] = min_rank;
			k_order[min_rank++] = now;
		}

	}



	int now_k_u = node_kpcore_info[u].size();
	int now_k_v = node_kpcore_info[v].size();

	int max_k = now_k_u > now_k_v ? now_k_u : now_k_v;

	if(now_k_v < now_k_u) {
		std::swap(u, v);
		std::swap(now_k_u, now_k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min + 1 && changed_nodes.size() > 0) {
			local_insert_maintain_p(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
		}
			
		else {
			local_insert_maintain_p(k, u, v, pre_k_u, pre_k_v, empty_v);
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
		}
	}

	return;

}

void Graph::local_insert_edge_woc(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	vector<vid_t> changed_nodes;
	vector<vid_t> new_k_order;

	// printf("before insert : ");
	// print_k_order();

	// add for k order insertion
	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}
	// end

	//inserted clip of my method -----------------
	int pre_k_u = k_u;
	int pre_k_v = k_v;
	//end of the clip ----------------------------

	int k_min = MIN(k_u, k_v); 
	insert_edge(u, v);


	// new method to maintain k value
	k_deg_plus[u] += 1;
	// printf("u %d, k_deg_plus = %d, k_mind = %d\n", u, k_deg_plus[u], k_min);

	if (k_deg_plus[u] <= k_min) {
		// puts("doneee");
		// done
	} else {
		// remember to reset k_deg_star value !
		// puts("IN else !");
		unordered_set<vid_t> vc;
		int now_rank = k_rank[u];
		int max_rank = k_order.size();
		while (now_rank < max_rank) {
			int now_id = k_order[now_rank];
			// printf("node_id : %d (+=%d, *=%d)\n", now_id, k_deg_plus[now_id], k_deg_star[now_id]);
			if (node_kpcore_info[k_order[now_rank]].size() > k_min) {
				break;
			} else {
				if (k_deg_star[now_id] + k_deg_plus[now_id] > k_min) {
					// puts("con1");
					vc.insert(now_id);
					for (vid_t nei : neighbor[now_id]) {
						if (node_kpcore_info[nei].size() == k_min && k_rank[nei] > now_rank) {
							k_deg_star[nei] ++;
						}
					}
					now_rank ++;
				} else if (k_deg_star[now_id] == 0) {
					// puts("con2");
					new_k_order.push_back(now_id);
					// k_rank[now_id] = determined_rank;
					now_rank ++;
				} else {
					// puts("con3");
					new_k_order.push_back(now_id);
					k_deg_plus[now_id] += k_deg_star[now_id];
					k_deg_star[now_id] = 0;

					vector<vid_t> que;
					unordered_set<vid_t> in_que;
					int head = 0;
					int tail = 0;
					for (vid_t nei : neighbor[now_id]) {
						// printf("nei %d, %d-%d\n", nei, k_deg_plus[nei], k_deg_star[nei]);
						if (vc.find(nei) != vc.end()) {
							k_deg_plus[nei] --;
							if (k_deg_plus[nei] + k_deg_star[nei] <= k_min) {
								que.push_back(nei);
								in_que.insert(nei);
								tail++;
							}
						}
					}

					while (head < tail) {
						// dequeue w
						int w = que[head++];
						// printf("%d out\n", w);

						// update k_deg_plus and k_deg_star
						k_deg_plus[w] += k_deg_star[w];
						k_deg_star[w] = 0;

						// delete w from Vc
						vc.erase(w);

						// append w to Ok
						new_k_order.push_back(w);
						// k_rank[w] = determined_rank;

						for (vid_t nei : neighbor[w]) {
							// printf("now's nei %d\n", nei);
							if (k_rank[nei] > k_rank[w]) {
								if(k_deg_star[nei] > 0) {
									k_deg_star[nei] --;
								}
								if (vc.find(nei) != vc.end()) {
									// printf("find! , +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
									// printf("%d\n", in_que.find(nei) == in_que.end() ? 0 : 1);
									if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
										que.push_back(nei);
										in_que.insert(nei);
										++ tail;
									}
								}
							}
							else if(vc.find(nei) != vc.end()) {
								k_deg_plus[nei] --;
								// printf("find2 ! +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
								if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
									que.push_back(nei);
									in_que.insert(nei);
									++ tail;
								}
							}
						}
					}

					now_rank ++;

				}
			}
			
		}

		vector<int> rank_sort;
	
		for (vid_t now : vc) {
			k_deg_star[now] = 0;
			node_kpcore_info[now].push_back({0,1});
			changed_nodes.push_back(now);
			rank[now].push_back(0);
			rank[now][k_min+1] = new_order[k_min + 1].size();
			new_order[k_min + 1].push_back(now);
			rank_sort.push_back(k_rank[now]);
		}

		sort(rank_sort.begin(), rank_sort.end());
		for (int temp : rank_sort) {
			new_k_order.push_back(k_order[temp]);
		}

		int min_rank = k_rank[u];
		for (vid_t now : new_k_order) {
			k_rank[now] = min_rank;
			k_order[min_rank++] = now;
		}

	}



	int now_k_u = node_kpcore_info[u].size();
	int now_k_v = node_kpcore_info[v].size();

	int max_k = now_k_u > now_k_v ? now_k_u : now_k_v;

	if(now_k_v < now_k_u) {
		std::swap(u, v);
		std::swap(now_k_u, now_k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min + 1 && changed_nodes.size() > 0) {
			local_insert_maintain_p_woc(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
		}
			
		else {
			local_insert_maintain_p_woc(k, u, v, pre_k_u, pre_k_v, empty_v);
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
		}
	}

	return;

}


vector<long long> Graph::local_delete_edge_vertex_statistics(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	double avg_check_time_per_vertex = 0;
	std::vector<long long> tempvec;

	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}

	int pre_k_u = k_u;
	int pre_k_v = k_v;

	int k_min = MIN(k_u, k_v); 
	remove_edge(u, v);

	if (k_u <= k_v) {
		k_mcd[u] --;
	}
	if(k_v <= k_u) {
		k_mcd[v] --;
	}

	bool kcore_change = false;

	// calculate mcd brute force
	unordered_map<vid_t, int> cd;
	vector<vid_t> V_star;
	unordered_set<vid_t> in_V_star;

	int head = 0;
	int tail = 0;

	cd[u] = k_mcd[u];
	cd[v] = k_mcd[v];

	if (cd[u] < k_min) {
		V_star.push_back(u);
		in_V_star.insert(u);
		tail++;
	}
	if (cd[v] < k_min) {
		V_star.push_back(v);
		in_V_star.insert(v);
		tail++;
	}

	if(head < tail) {
		
		kcore_change = true;

		while (head < tail) {
			vid_t now = V_star[head++];
			for (int nei : neighbor[now]) {
				if (in_V_star.find(nei) != in_V_star.end()) {
					continue;
				}
				if (node_kpcore_info[nei].size() == k_min) {
					if (cd.find(nei) != cd.end()) {
						cd[nei] --;
					}
					else {
						cd[nei] = k_mcd[nei] - 1;
					}
					if(cd[nei] < k_min) {
						V_star.push_back(nei);
						in_V_star.insert(nei);
						tail++;
					}
				}
			}
		}

		int num = V_star.size();
		for (int i = 0; i < num; i++) {
			vid_t now = V_star[i];
			k_deg_plus[now] = 0;
			for (vid_t nei : neighbor[now]) {
				if(node_kpcore_info[nei].size() == k_min){
					if(in_V_star.find(nei) == in_V_star.end()) k_mcd[nei] --;
					if (k_rank[nei] < k_rank[now]) k_deg_plus[nei] --;
				} 
				if(node_kpcore_info[nei].size() >= k_min || in_V_star.find(nei) != in_V_star.end()) {
					k_deg_plus[now] ++;
				}
				if(node_kpcore_info[nei].size() == k_min - 1) {
					k_mcd[now]++;
				}
			}
			in_V_star.erase(now);
		}

		int begin_rank = k_rank[V_star[0]];
		int end_rank = k_rank[V_star[0]];

		while(begin_rank >= 0 && node_kpcore_info[k_order[begin_rank]].size() == k_min) begin_rank --;
		while(end_rank < k_order.size() && node_kpcore_info[k_order[begin_rank]].size() == k_min) end_rank ++;
		begin_rank ++;
		end_rank --;

		for (vid_t now : V_star) {
			in_V_star.insert(now);
		}

		int determined_rank = end_rank + 1;
		for (int i = end_rank; i >= begin_rank; i--) {
			if (in_V_star.find(k_order[i]) == in_V_star.end()) {
				k_order[--determined_rank] = k_order[i];
				k_rank[k_order[determined_rank]] = determined_rank;
			}
		}

		for (int i = 0; i < num; i++){
			k_order[begin_rank + i] = V_star[i];
			k_rank[V_star[i]] = begin_rank + i;
		}

		if (k_min == 1) {
			for (vid_t now : V_star) {
				node_kpcore_info[now].pop_back();
			}
		}
	}

	// puts("After k update");
	
	// need to deal with the condition that p_order[min_k]'s size equal to 0
	// after the deletion (there's no node in min_k core)
	
	int max_k =(k_u <  k_v) ? k_v : k_u;

	if(k_rank[v] < k_rank[u]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min && kcore_change) {
			// printf("new k = %d, size = %d\n", k_min, changed_nodes.size()); 
			// print_new_order(k);
			tempvec = local_delete_maintain_p_for_changed_core_vertex_statistics(k, u, v, pre_k_u, pre_k_v, V_star);
			// print_new_order(k);
		}
		else {
			// printf("else k = %d\n", k);
			// print_new_order(k);
			tempvec = local_delete_maintain_p_vertex_statistics(k, u, v, pre_k_u, pre_k_v); // revise
			// print_new_order(k);
		}

		dec_times += tempvec[0];
		checknum_nei_dec += tempvec[1];
		checknum_dec_dec += tempvec[2];
		total_check_number_dec += tempvec[3];
		distinct_check_num_dec += tempvec[4];
		inc_times += tempvec[5];
		checknum_nei_inc += tempvec[6];
		checknum_dec_inc += tempvec[7];

		if(tempvec[0] == 1 && tempvec[4] == 0) dec_times--;

		//{dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
		if (tempvec[0] == 1 && tempvec[4] != 0) {
			avg_check_time_per_vertex += 1.0 * tempvec[3] / tempvec[4];
			if(tempvec[4] == 0) {
				cout << " " << tempvec[3] << " " << tempvec[4] << endl;
				exit(0);
			}
			// cout << avg_check_time_per_vertex << endl;
		}
	}
	// cout << "add " << avg_check_time_per_vertex << endl;
	return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc, (long long)(avg_check_time_per_vertex * 1000)};


}

vector<long long> Graph::local_insert_edge_vertex_statistics(vid_t u, vid_t v) {
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size(); 

	vector<vid_t> changed_nodes;
	vector<vid_t> new_k_order;

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	double avg_check_time_per_vertex = 0;
	std::vector<long long> tempvec;

	// printf("before insert : ");
	// print_k_order();

	// add for k order insertion
	if(k_rank[u] > k_rank[v]) {
		std::swap(u, v);
		std::swap(k_u, k_v);
	}
	// end

	//inserted clip of my method -----------------
	int pre_k_u = k_u;
	int pre_k_v = k_v;
	//end of the clip ----------------------------

	int k_min = MIN(k_u, k_v); 
	insert_edge(u, v);


	// new method to maintain k value
	k_deg_plus[u] += 1;
	// printf("u %d, k_deg_plus = %d, k_mind = %d\n", u, k_deg_plus[u], k_min);

	if (k_deg_plus[u] <= k_min) {
		// puts("doneee");
		// done
	} else {
		// remember to reset k_deg_star value !
		// puts("IN else !");
		unordered_set<vid_t> vc;
		int now_rank = k_rank[u];
		int max_rank = k_order.size();
		while (now_rank < max_rank) {
			int now_id = k_order[now_rank];
			// printf("node_id : %d (+=%d, *=%d)\n", now_id, k_deg_plus[now_id], k_deg_star[now_id]);
			if (node_kpcore_info[k_order[now_rank]].size() > k_min) {
				break;
			} else {
				if (k_deg_star[now_id] + k_deg_plus[now_id] > k_min) {
					// puts("con1");
					vc.insert(now_id);
					for (vid_t nei : neighbor[now_id]) {
						if (node_kpcore_info[nei].size() == k_min && k_rank[nei] > now_rank) {
							k_deg_star[nei] ++;
						}
					}
					now_rank ++;
				} else if (k_deg_star[now_id] == 0) {
					// puts("con2");
					new_k_order.push_back(now_id);
					// k_rank[now_id] = determined_rank;
					now_rank ++;
				} else {
					// puts("con3");
					new_k_order.push_back(now_id);
					k_deg_plus[now_id] += k_deg_star[now_id];
					k_deg_star[now_id] = 0;

					vector<vid_t> que;
					unordered_set<vid_t> in_que;
					int head = 0;
					int tail = 0;
					for (vid_t nei : neighbor[now_id]) {
						// printf("nei %d, %d-%d\n", nei, k_deg_plus[nei], k_deg_star[nei]);
						if (vc.find(nei) != vc.end()) {
							k_deg_plus[nei] --;
							if (k_deg_plus[nei] + k_deg_star[nei] <= k_min) {
								que.push_back(nei);
								in_que.insert(nei);
								tail++;
							}
						}
					}

					while (head < tail) {
						// dequeue w
						int w = que[head++];
						// printf("%d out\n", w);

						// update k_deg_plus and k_deg_star
						k_deg_plus[w] += k_deg_star[w];
						k_deg_star[w] = 0;

						// delete w from Vc
						vc.erase(w);

						// append w to Ok
						new_k_order.push_back(w);
						// k_rank[w] = determined_rank;

						for (vid_t nei : neighbor[w]) {
							// printf("now's nei %d\n", nei);
							if (k_rank[nei] > k_rank[w]) {
								if(k_deg_star[nei] > 0) {
									k_deg_star[nei] --;
								}
								if (vc.find(nei) != vc.end()) {
									// printf("find! , +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
									// printf("%d\n", in_que.find(nei) == in_que.end() ? 0 : 1);
									if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
										que.push_back(nei);
										in_que.insert(nei);
										++ tail;
									}
								}
							}
							else if(vc.find(nei) != vc.end()) {
								k_deg_plus[nei] --;
								// printf("find2 ! +=%d, *=%d\n", k_deg_plus[nei], k_deg_star[nei]);
								if (k_deg_star[nei] + k_deg_plus[nei] <= k_min && in_que.find(nei) == in_que.end()) {
									que.push_back(nei);
									in_que.insert(nei);
									++ tail;
								}
							}
						}
					}

					now_rank ++;

				}
			}
			
		}

		vector<int> rank_sort;
	
		for (vid_t now : vc) {
			k_deg_star[now] = 0;
			node_kpcore_info[now].push_back({0,1});
			changed_nodes.push_back(now);
			rank[now].push_back(0);
			rank[now][k_min+1] = new_order[k_min + 1].size();
			new_order[k_min + 1].push_back(now);
			rank_sort.push_back(k_rank[now]);
		}

		sort(rank_sort.begin(), rank_sort.end());
		for (int temp : rank_sort) {
			new_k_order.push_back(k_order[temp]);
		}

		int min_rank = k_rank[u];
		for (vid_t now : new_k_order) {
			k_rank[now] = min_rank;
			k_order[min_rank++] = now;
		}

	}



	int now_k_u = node_kpcore_info[u].size();
	int now_k_v = node_kpcore_info[v].size();

	int max_k = now_k_u > now_k_v ? now_k_u : now_k_v;

	if(now_k_v < now_k_u) {
		std::swap(u, v);
		std::swap(now_k_u, now_k_v);
		std::swap(pre_k_u, pre_k_v);
	}

	vector<vid_t> empty_v;
	for(int k = 2; k <= max_k; k++) {
		// printf("%d\n", k);
		if(k == k_min + 1 && changed_nodes.size() > 0) {
			tempvec = local_insert_maintain_p_vertex_statistics(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			
			
			
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, changed_nodes);
		}
			
		else {
			tempvec = local_insert_maintain_p_vertex_statistics(k, u, v, pre_k_u, pre_k_v, empty_v);
			
			// if(k < 20) update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
			// else update_p_value_by_order_insertion(k, u, v, pre_k_u, pre_k_v, empty_v);
		}
		dec_times += tempvec[0];
		checknum_nei_dec += tempvec[1];
		checknum_dec_dec += tempvec[2];
		total_check_number_dec += tempvec[3];
		distinct_check_num_dec += tempvec[4];
		inc_times += tempvec[5];
		checknum_nei_inc += tempvec[6];
		checknum_dec_inc += tempvec[7];
		//{dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
		if (tempvec[0] == 1) {
			avg_check_time_per_vertex += 1.0 * tempvec[3] / tempvec[4];
		}
	}

	return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc, (long long)(avg_check_time_per_vertex * 1000)};
}


bool Graph::check_increase_delete(int k, int u) {
	int deg = degree[u];
	int supporter = 0;
	p_t now_p = node_kpcore_info[u][k - 1];
	for(int nei : neighbor[u]) {
		if(node_kpcore_info[nei].size() < k) {
			continue;
		}
		int compare = compare_pvalue(node_kpcore_info[nei][k - 1], now_p);
		if(compare > 0) {
			supporter ++;
		}
		else if(compare < 0) {
			continue;
		}
		else if(compare == 0) {
			supporter++;
			p_t p_nei = node_kpcore_info[nei][k - 1];
			int supporter_nei = 0;
			for(int neinei : neighbor[nei]) {
				if(node_kpcore_info[neinei].size() >= k && 
					compare_pvalue(node_kpcore_info[neinei][k - 1], p_nei) >= 0) {
						supporter_nei ++;
				}
			}
			if(compare_pvalue(now_p, (p_t){supporter *1ll, deg*1ll}) <= 0) {
				supporter++;
			}
		}
	}
	if(compare_pvalue({supporter, deg}, now_p) > 0) return true;
	else return false;
}

bool Graph::check_decrease(int k, int u, int new_link_to) {
	int deg = degree[u];
	int supporter = 0;
	p_t now_p = node_kpcore_info[u][k - 1];
	for(int nei : neighbor[u]) {
		if(node_kpcore_info[nei].size() >= k && compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0) {
			supporter ++;
		}
	}
	if(compare_pvalue((p_t){supporter *1ll, deg*1ll}, now_p) < 0 || supporter < k) return true;
	else return false;
}

void Graph::pIncrease(int k, int begin_rank) {
	
	// after checked, add a new tag to rank to show where p changes, instead using the compare_p() function
	// change order and rank to pointers
	// change supports and inqueue to int[], use maintain_cnt to represent inqueue

	// for change_tag : when a node is in queue, set its change_tag be zero

	// printf("k=%d\n", k);
	// print_new_order(k);

	std::unordered_map<vid_t, int> supports;
	std::unordered_map<vid_t, bool> dirty;	

	min_heap qu;
	qu.node2index.clear();
	qu.node2index.resize(num_node + 10, -1);
	int* que = (int*) malloc((num_node+10) * sizeof(int));
	// std::vector<int>que;
	// que.resize(num_node + 10);

	std::vector<int> *order_k = &(new_order[k]);

	int next_rank = begin_rank;
	int determined_rank = begin_rank + 1;

	

	// int undecided_cnt = 0;
	// puts("A");

	int begin_node = (*order_k)[begin_rank];
	// printf("begin %d %lld/%lld\n", begin_node, node_kpcore_info[begin_node][k - 1].first, node_kpcore_info[begin_node][k - 1].second);
	if(node_kpcore_info[begin_node][k - 1].first == 0) {
		for(int rk = next_rank; rk >= 0 && node_kpcore_info[(*order_k)[rk]][k - 1].first == 0; rk --) {
			// printf("rk %d\n", rk);
			change_tag[(*order_k)[rk]] = 1;
		}
	}
 
	// puts("B");

	int support_begin_node = 0;
	for(vid_t nei : neighbor[begin_node]) {
		if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
			support_begin_node ++;
		}
	}
	if(support_begin_node >= k && compare_pvalue((p_t){support_begin_node, degree[begin_node]}, node_kpcore_info[begin_node][k - 1]) > 0) {
		qu.insert(std::make_pair((p_t){support_begin_node, degree[begin_node]}, begin_node));
		for(vid_t nei : neighbor[begin_node]) {
			if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
				change_tag[nei] ++;
			}
		}
		// printf("in %d\n", begin_node);
		supports[begin_node] = support_begin_node;
		next_rank --;
	}

	// std::vector<int>::iterator order_k = new_order[k].begin();

	while(next_rank + 1 < determined_rank) {

		p_t next_p = (next_rank >= 0) ? node_kpcore_info[(*order_k)[next_rank]][k - 1] : (p_t){1, 1};

		//  printf("next_p %lld/%lld\n", next_p.first, next_p.second);

		if(next_p.first != 0 && compare_pvalue(qu.top().first, next_p) <= 0) {
			// case collapse
			 
			p_t now_p = qu.top().first;
			// printf("collapse now_p %lld/%lld\n", now_p.first, now_p.second);
			int head = 0;
			int tail = 0;
			que[tail++] = qu.top().second;
			supports.erase(que[head]);
			while(head < tail) {
				vid_t now = que[head++];
				// printf("out %d\n", now);
				(*order_k)[--determined_rank] = now;
				// printf("deter %d, next %d\n", determined_rank, next_rank);
				rank[now][k] = determined_rank;
				node_kpcore_info[now][k - 1] = now_p;
				// supports.erase(now);
				qu.remove(qu.node2index[now]);

				for(vid_t nei : neighbor[now]) {
					if(supports.find(nei) != supports.end()) {
						supports[nei] --;
						if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
							que[tail++] = nei;
							supports.erase(nei);
							dirty[nei] = false;
						}
						else {
							dirty[nei] = true;
						}
					}
					else if(change_tag[nei] > 0) {
						change_tag[nei] --;
					}
				}
			}

			for(int i = 0; i < tail; i++){
				vid_t now = que[i];
				for(vid_t nei : neighbor[now]) {
					if(dirty[nei] == true) {
						qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
						dirty[nei] = false;
					}
				}
			}
		}
		else {
			int next_node = (*order_k)[next_rank];
			// printf("next_node : %d, tag = %d\n", next_node, change_tag[next_node]);
			// printf("%lld/%lld\n", node_kpcore_info[next_node][k - 1].first, node_kpcore_info[next_node][k - 1].second);
			if(change_tag[next_node] == 0){
				// puts("con change_tag[next_node] == 0");
				(*order_k)[--determined_rank] = next_node;
				rank[next_node][k] = determined_rank;
				next_rank --;
			}
			else {
				change_tag[next_node] = 0;
				int temp = 0;
				for(vid_t nei : neighbor[next_node]) {
					// printf("nei %d\n", nei);
					if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
						change_tag[nei] ++;
						temp ++;
					}
					else if (supports.find(nei) != supports.end()) {
						temp ++;
					}
					
				}
				// printf("temp : %d\n", temp);
				if(temp >= k && compare_pvalue((p_t){temp, degree[next_node]}, node_kpcore_info[next_node][k - 1]) > 0) {
					qu.insert(std::make_pair((p_t){temp, degree[next_node]}, next_node));
					// for(vid_t nei : neighbor[next_node]) {
					// 	printf("checknei %d\n", nei);
					// 	if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
							
					// 		printf("nei : %d, tag = %d\n", nei, change_tag[nei]);
					// 	}
					// }
					supports[next_node] = temp;
					// printf("in %d\n", next_node);
				}
				else {
					// puts("con collapse");
					p_t now_p = node_kpcore_info[next_node][k - 1];
					// (*order_k)[--determined_rank] = next_node;
					// rank[next_node][k] = determined_rank;
					int head = 0, tail = 0;
					que[tail ++] = next_node;


					while(head < tail) {
						vid_t now = que[head++];
						(*order_k)[--determined_rank] = now;
						rank[now][k] = determined_rank;
						node_kpcore_info[now][k - 1] = now_p;
						
						qu.remove(qu.node2index[now]);

						for(vid_t nei : neighbor[now]) {
							if(supports.find(nei) != supports.end()) {
								supports[nei] --;
								if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
									que[tail++] = nei;
									supports.erase(nei);
									dirty[nei] = false;
								}
								else {
									dirty[nei] = true;
								}
							}
							else if(change_tag[nei] > 0) {
								change_tag[nei] --;
							}
						}
					}

					for(int i = 0; i < tail; i++){
						vid_t now = que[i];
						// printf("check now %d\n", now);
						for(vid_t nei : neighbor[now]) {
							if(dirty[nei] == true) {
								qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
								dirty[nei] = false;
							}
						}
					}
				}

				next_rank --;
			}
		}

	}

	free(que);

	return;

}

void Graph::pIncrease_woc(int k, int begin_rank) {
	
	// after checked, add a new tag to rank to show where p changes, instead using the compare_p() function
	// change order and rank to pointers
	// change supports and inqueue to int[], use maintain_cnt to represent inqueue

	// for change_tag : when a node is in queue, set its change_tag be zero

	// printf("k=%d\n", k);
	// print_new_order(k);

	std::unordered_map<vid_t, int> supports;
	std::unordered_map<vid_t, bool> dirty;	

	min_heap qu;
	qu.node2index.clear();
	qu.node2index.resize(num_node + 10, -1);
	int* que = (int*) malloc((num_node+10) * sizeof(int));
	// std::vector<int>que;
	// que.resize(num_node + 10);

	std::vector<int> *order_k = &(new_order[k]);

	int next_rank = begin_rank;
	int determined_rank = begin_rank + 1;

	

	// int undecided_cnt = 0;
	// puts("A");

	int begin_node = (*order_k)[begin_rank];
	// printf("begin %d %lld/%lld\n", begin_node, node_kpcore_info[begin_node][k - 1].first, node_kpcore_info[begin_node][k - 1].second);
	if(node_kpcore_info[begin_node][k - 1].first == 0) {
		for(int rk = next_rank; rk >= 0 && node_kpcore_info[(*order_k)[rk]][k - 1].first == 0; rk --) {
			// printf("rk %d\n", rk);
			change_tag[(*order_k)[rk]] = 1;
		}
	}
 
	// puts("B");

	int support_begin_node = 0;
	for(vid_t nei : neighbor[begin_node]) {
		if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
			support_begin_node ++;
		}
	}
	if(support_begin_node >= k && compare_pvalue((p_t){support_begin_node, degree[begin_node]}, node_kpcore_info[begin_node][k - 1]) > 0) {
		qu.insert(std::make_pair((p_t){support_begin_node, degree[begin_node]}, begin_node));
		for(vid_t nei : neighbor[begin_node]) {
			if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
				change_tag[nei] ++;
			}
		}
		// printf("in %d\n", begin_node);
		supports[begin_node] = support_begin_node;
		next_rank --;
	}

	// std::vector<int>::iterator order_k = new_order[k].begin();

	while(next_rank + 1 < determined_rank) {

		p_t next_p = (next_rank >= 0) ? node_kpcore_info[(*order_k)[next_rank]][k - 1] : (p_t){1, 1};

		//  printf("next_p %lld/%lld\n", next_p.first, next_p.second);

		if(next_p.first != 0 && compare_pvalue(qu.top().first, next_p) <= 0) {
			// case collapse
			 
			p_t now_p = qu.top().first;
			// printf("collapse now_p %lld/%lld\n", now_p.first, now_p.second);
			int head = 0;
			int tail = 0;
			que[tail++] = qu.top().second;
			supports.erase(que[head]);
			while(head < tail) {
				vid_t now = que[head++];
				// printf("out %d\n", now);
				(*order_k)[--determined_rank] = now;
				// printf("deter %d, next %d\n", determined_rank, next_rank);
				rank[now][k] = determined_rank;
				node_kpcore_info[now][k - 1] = now_p;
				// supports.erase(now);
				qu.remove(qu.node2index[now]);

				for(vid_t nei : neighbor[now]) {
					if(supports.find(nei) != supports.end()) {
						supports[nei] --;
						if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
							que[tail++] = nei;
							supports.erase(nei);
							dirty[nei] = false;
						}
						else {
							dirty[nei] = true;
						}
					}
					else if(change_tag[nei] > 0) {
						change_tag[nei] --;
					}
				}
			}

			for(int i = 0; i < tail; i++){
				vid_t now = que[i];
				for(vid_t nei : neighbor[now]) {
					if(dirty[nei] == true) {
						qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
						dirty[nei] = false;
					}
				}
			}
		}
		else {
			int next_node = (*order_k)[next_rank];
			// printf("next_node : %d, tag = %d\n", next_node, change_tag[next_node]);
			// printf("%lld/%lld\n", node_kpcore_info[next_node][k - 1].first, node_kpcore_info[next_node][k - 1].second);
			// if(change_tag[next_node] == 0){
			// 	// puts("con change_tag[next_node] == 0");
			// 	(*order_k)[--determined_rank] = next_node;
			// 	rank[next_node][k] = determined_rank;
			// 	next_rank --;
			// }
			// else {
				change_tag[next_node] = 0;
				int temp = 0;
				for(vid_t nei : neighbor[next_node]) {
					// printf("nei %d\n", nei);
					if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
						change_tag[nei] ++;
						temp ++;
					}
					else if (supports.find(nei) != supports.end()) {
						temp ++;
					}
					
				}
				// printf("temp : %d\n", temp);
				if(temp >= k && compare_pvalue((p_t){temp, degree[next_node]}, node_kpcore_info[next_node][k - 1]) > 0) {
					qu.insert(std::make_pair((p_t){temp, degree[next_node]}, next_node));
					// for(vid_t nei : neighbor[next_node]) {
					// 	printf("checknei %d\n", nei);
					// 	if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
							
					// 		printf("nei : %d, tag = %d\n", nei, change_tag[nei]);
					// 	}
					// }
					supports[next_node] = temp;
					// printf("in %d\n", next_node);
				}
				else {
					// puts("con collapse");
					p_t now_p = node_kpcore_info[next_node][k - 1];
					// (*order_k)[--determined_rank] = next_node;
					// rank[next_node][k] = determined_rank;
					int head = 0, tail = 0;
					que[tail ++] = next_node;

					while(head < tail) {
						vid_t now = que[head++];
						(*order_k)[--determined_rank] = now;
						rank[now][k] = determined_rank;
						node_kpcore_info[now][k - 1] = now_p;
						
						qu.remove(qu.node2index[now]);

						for(vid_t nei : neighbor[now]) {
							if(supports.find(nei) != supports.end()) {
								supports[nei] --;
								if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
									que[tail++] = nei;
									supports.erase(nei);
									dirty[nei] = false;
								}
								else {
									dirty[nei] = true;
								}
							}
							else if(change_tag[nei] > 0) {
								change_tag[nei] --;
							}
						}
					}

					for(int i = 0; i < tail; i++){
						vid_t now = que[i];
						// printf("check now %d\n", now);
						for(vid_t nei : neighbor[now]) {
							if(dirty[nei] == true) {
								qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
								dirty[nei] = false;
							}
						}
					}
				}

				next_rank --;
			// }
		}

	}

	free(que);

	return;

}


vector<long long> Graph::pIncrease_vertex_statistics(int k, int begin_rank) {
	
	// after checked, add a new tag to rank to show where p changes, instead using the compare_p() function
	// change order and rank to pointers
	// change supports and inqueue to int[], use maintain_cnt to represent inqueue

	// for change_tag : when a node is in queue, set its change_tag be zero

	// printf("k=%d\n", k);
	// print_new_order(k);

	std::unordered_map<vid_t, int> supports;
	std::unordered_map<vid_t, bool> dirty;	

	long long checknum_nei = 1;
	long long checknum_dec = 1;

	min_heap qu;
	qu.node2index.clear();
	qu.node2index.resize(num_node + 10, -1);
	int* que = (int*) malloc((num_node+10) * sizeof(int));
	// std::vector<int>que;
	// que.resize(num_node + 10);

	std::vector<int> *order_k = &(new_order[k]);

	int next_rank = begin_rank;
	int determined_rank = begin_rank + 1;

	

	// int undecided_cnt = 0;
	// puts("A");

	int begin_node = (*order_k)[begin_rank];
	// printf("begin %d %lld/%lld\n", begin_node, node_kpcore_info[begin_node][k - 1].first, node_kpcore_info[begin_node][k - 1].second);
	if(node_kpcore_info[begin_node][k - 1].first == 0) {
		for(int rk = next_rank; rk >= 0 && node_kpcore_info[(*order_k)[rk]][k - 1].first == 0; rk --) {
			// printf("rk %d\n", rk);
			change_tag[(*order_k)[rk]] = 1;
		}
	}
 
	// puts("B");

	int support_begin_node = 0;
	for(vid_t nei : neighbor[begin_node]) {
		if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
			support_begin_node ++;
		}
	}
	if(support_begin_node >= k && compare_pvalue((p_t){support_begin_node, degree[begin_node]}, node_kpcore_info[begin_node][k - 1]) > 0) {
		qu.insert(std::make_pair((p_t){support_begin_node, degree[begin_node]}, begin_node));
		for(vid_t nei : neighbor[begin_node]) {
			if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= begin_rank) {
				change_tag[nei] ++;
			}
		}
		// printf("in %d\n", begin_node);
		supports[begin_node] = support_begin_node;
		next_rank --;
	}

	// std::vector<int>::iterator order_k = new_order[k].begin();

	while(next_rank + 1 < determined_rank) {

		p_t next_p = (next_rank >= 0) ? node_kpcore_info[(*order_k)[next_rank]][k - 1] : (p_t){1, 1};

		//  printf("next_p %lld/%lld\n", next_p.first, next_p.second);

		if(next_p.first != 0 && compare_pvalue(qu.top().first, next_p) <= 0) {
			// case collapse
			 
			p_t now_p = qu.top().first;
			// printf("collapse now_p %lld/%lld\n", now_p.first, now_p.second);
			int head = 0;
			int tail = 0;
			que[tail++] = qu.top().second;
			supports.erase(que[head]);
			while(head < tail) {
				vid_t now = que[head++];
				// printf("out %d\n", now);
				(*order_k)[--determined_rank] = now;
				// printf("deter %d, next %d\n", determined_rank, next_rank);
				rank[now][k] = determined_rank;
				node_kpcore_info[now][k - 1] = now_p;
				// supports.erase(now);
				qu.remove(qu.node2index[now]);

				for(vid_t nei : neighbor[now]) {
					if(supports.find(nei) != supports.end()) {
						supports[nei] --;
						if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
							que[tail++] = nei;
							supports.erase(nei);
							dirty[nei] = false;
						}
						else {
							dirty[nei] = true;
						}
					}
					else if(change_tag[nei] > 0) {
						change_tag[nei] --;
					}
				}
			}

			for(int i = 0; i < tail; i++){
				vid_t now = que[i];
				for(vid_t nei : neighbor[now]) {
					if(dirty[nei] == true) {
						qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
						dirty[nei] = false;
					}
				}
			}
		}
		else {
			int next_node = (*order_k)[next_rank];
			// printf("next_node : %d, tag = %d\n", next_node, change_tag[next_node]);
			// printf("%lld/%lld\n", node_kpcore_info[next_node][k - 1].first, node_kpcore_info[next_node][k - 1].second);
			if(change_tag[next_node] == 0){
				// puts("con change_tag[next_node] == 0");
				(*order_k)[--determined_rank] = next_node;
				rank[next_node][k] = determined_rank;
				next_rank --;
			}
			else {
				checknum_nei++;
				change_tag[next_node] = 0;
				int temp = 0;
				for(vid_t nei : neighbor[next_node]) {
					// printf("nei %d\n", nei);
					if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
						change_tag[nei] ++;
						temp ++;
					}
					else if (supports.find(nei) != supports.end()) {
						temp ++;
					}
					
				}
				// printf("temp : %d\n", temp);
				if(temp >= k && compare_pvalue((p_t){temp, degree[next_node]}, node_kpcore_info[next_node][k - 1]) > 0) {
					qu.insert(std::make_pair((p_t){temp, degree[next_node]}, next_node));
					// for(vid_t nei : neighbor[next_node]) {
					// 	printf("checknei %d\n", nei);
					// 	if(node_kpcore_info[nei].size() >= k && rank[nei][k] <= next_rank) {
							
					// 		printf("nei : %d, tag = %d\n", nei, change_tag[nei]);
					// 	}
					// }
					supports[next_node] = temp;

					checknum_dec ++;
					// printf("in %d\n", next_node);
				}
				else {
					// puts("con collapse");
					p_t now_p = node_kpcore_info[next_node][k - 1];
					// (*order_k)[--determined_rank] = next_node;
					// rank[next_node][k] = determined_rank;
					int head = 0, tail = 0;
					que[tail ++] = next_node;
					// for(vid_t nei : neighbor[next_node]) {
					// 	if(supports.find(nei) != supports.end()) {
					// 		supports[nei] --;
					// 		if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
					// 			que[tail++] = nei;
					// 			supports.erase(nei);
					// 			dirty[nei] = false;
					// 		}
					// 		else {
					// 			dirty[nei] = true;
					// 			printf("dirty ! %d (%d)\n", nei, supports[nei]);
					// 		}
					// 	}
					// }
					while(head < tail) {
						vid_t now = que[head++];
						(*order_k)[--determined_rank] = now;
						rank[now][k] = determined_rank;
						node_kpcore_info[now][k - 1] = now_p;
						
						qu.remove(qu.node2index[now]);

						for(vid_t nei : neighbor[now]) {
							if(supports.find(nei) != supports.end()) {
								supports[nei] --;
								if(supports[nei] < k || compare_pvalue((p_t){supports[nei], degree[nei]}, now_p) <= 0) {
									que[tail++] = nei;
									supports.erase(nei);
									dirty[nei] = false;
								}
								else {
									dirty[nei] = true;
								}
							}
							else if(change_tag[nei] > 0) {
								change_tag[nei] --;
							}
						}
					}

					for(int i = 0; i < tail; i++){
						vid_t now = que[i];
						// printf("check now %d\n", now);
						for(vid_t nei : neighbor[now]) {
							if(dirty[nei] == true) {
								qu.decrease_key((p_t){supports[nei], degree[nei]}, qu.node2index[nei]);
								dirty[nei] = false;
							}
						}
					}
				}

				next_rank --;
			}
		}

	}

	free(que);

	return {checknum_dec, checknum_nei};

}

// traverse to find decrease set using some vertices that must decrease
void Graph::get_decrease_set(int k, p_t p, std::unordered_set<vid_t> &decrease_set) {
	vector<vid_t> que;
	int head = 0;
	int tail = 0;
	for (vid_t v : decrease_set) {
		que.push_back(v);
		++tail;
	}

	unordered_map<vid_t, int> deg_ge;

	while (head < tail) {
		int now = que[head++];
		for (vid_t nei : neighbor[now]) {
			if (node_kpcore_info[nei].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[nei][k - 1], p) == 0 && decrease_set.find(nei) == decrease_set.end()) {

				if (deg_ge.find(nei) != deg_ge.end()) {
					deg_ge[nei] --;
				} else {
					int deg_temp = 0;
					for (vid_t w : neighbor[nei]) {
						if (node_kpcore_info[w].size() >= k && compare_pvalue(node_kpcore_info[w][k - 1], p) >= 0) {
							deg_temp ++;
						}
					}
					deg_temp --;
					deg_ge[nei] = deg_temp;
				}
				
				if (compare_pvalue((p_t){deg_ge[nei], degree[nei]}, p) < 0 || deg_ge[nei] < k) {
					decrease_set.insert(nei);
					que.push_back(nei);
					tail++;
				}
			}
		}
	}
}

long long Graph::get_decrease_set_vertex_statistics(int k, p_t p, std::unordered_set<vid_t> &decrease_set) {
	vector<vid_t> que;
	int head = 0;
	int tail = 0;
	for (vid_t v : decrease_set) {
		que.push_back(v);
		++tail;
	}

	int checknum = 0;

	unordered_map<vid_t, int> deg_ge;

	while (head < tail) {
		int now = que[head++];
		for (vid_t nei : neighbor[now]) {
			if (node_kpcore_info[nei].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[nei][k - 1], p) == 0 && decrease_set.find(nei) == decrease_set.end()) {

				if (deg_ge.find(nei) != deg_ge.end()) {
					deg_ge[nei] --;
				} else {
					int deg_temp = 0;
					checknum++;
					for (vid_t w : neighbor[nei]) {
						if (node_kpcore_info[w].size() >= k && compare_pvalue(node_kpcore_info[w][k - 1], p) >= 0) {
							deg_temp ++;
						}
					}
					deg_temp --;
					deg_ge[nei] = deg_temp;
				}
				
				if (compare_pvalue((p_t){deg_ge[nei], degree[nei]}, p) < 0 || deg_ge[nei] < k) {
					decrease_set.insert(nei);
					que.push_back(nei);
					tail++;
				}
			}
		}
	}

	return checknum;
}

void Graph::pDecrease_multiple_vertices(int k, std::vector<vid_t> othernodes) {

	// please check v must decrease and do this 

	// printf("in decrease k=%d, othernodes.size=%d\n", k, othernodes.size());
	// print_new_order(k);

	int v = othernodes[0];

	for (vid_t u : othernodes) {
		// printf("%d ", u);
		change_tag[u] = 1;
		if(rank[u][k] < rank[v][k]) {
			v = u;
		}
	}
	// puts("");



	unordered_set<vid_t> decrease_set;
	unordered_set<vid_t> vc;
	// unordered_map<vid_t, int> deg_plus;
	unordered_map<vid_t, int> deg_plus_vc;

	// decrease_set.insert(v);
	
	// max_heap qu;
	// qu.node2index.clear();
	// qu.node2index.resize(num_node + 10, -1);

	min_heap min_qu;
	min_qu.node2index.clear();
	min_qu.node2index.resize(num_node + 10, -1);

	int begin_rank = rank[v][k];
	int end_rank = rank[v][k];
	int determined_rank;
	while (begin_rank > 0 && 
			compare_pvalue(node_kpcore_info[new_order[k][begin_rank]][k - 1], 
			node_kpcore_info[new_order[k][begin_rank - 1]][k - 1]) == 0 ) {
			
			begin_rank --;
	}
	end_rank = begin_rank;
	determined_rank = begin_rank - 1;
	int last_begin_rank = begin_rank;

	// bool first_loop = true;

	// [begin_rank, end_rank)

	int rank_size = new_order[k].size();

	int times = 0;

	vector<vid_t> decreased_order;
	vector<vid_t> remain_order;



	do {

		// printf("end_rank %d(node = %d)\n", end_rank, new_order[k][end_rank]);

		vector<vid_t> Oreverse;

		begin_rank = end_rank;
		p_t now_p = node_kpcore_info[new_order[k][begin_rank]][k - 1];

		while (end_rank < rank_size && 
		  compare_pvalue(node_kpcore_info[new_order[k][end_rank]][k - 1], now_p) == 0) {
			end_rank ++;
		}
 
		// printf("begin_rank = %d, end_rank = %d\n", begin_rank, end_rank);

		unordered_set<vid_t> determined_set;

		vector<vid_t> linked_node;


		for (int now_rank = begin_rank; now_rank < end_rank; now_rank ++) {
			int node = new_order[k][now_rank];

			if (change_tag[node] > 0) {
				change_tag[node] = 0;
				linked_node.push_back(node);
				for (vid_t w : neighbor[node]) {
					if (vc.find(w) != vc.end()) {
						// deg_plus[w]++;
						deg_plus_vc[w]++;
					}
				}
			}

		}
		// printf("linked_node.size = %d nowp %lld/%lld\n", linked_node.size(), now_p.first, now_p.second);

		if (linked_node.size() != 0) {

			for (vid_t node : vc) {
				// printf("node %d\n", node);
				if(deg_plus_vc[node] >= k && compare_pvalue({deg_plus_vc[node], degree[node]}, now_p) >= 0) {
					// puts("qu");
					// qu.remove(qu.node2index[node]);
					// puts("dp");
					// deg_plus.erase(node);
					// puts("dpv");
					deg_plus_vc.erase(node);
					// puts("ds");
					determined_set.insert(node);
					// puts("AA");
				}
			}

			// printf("determined : ");
			// for (vid_t node : determined_set) {
			// 	printf("%d ", node);
			// }
			// puts("");

			for (vid_t node : determined_set) {
				vc.erase(node);
				node_kpcore_info[node][k - 1] = now_p;
			}

			
			// printf("(%d-%d) %lld/%lld determined : %d, linked_node = %d\n", ++times, k, now_p.first, now_p.second, determined_set.size(), linked_node.size());
			
			
			decrease_set.clear();

			// if (first_loop) {
			// 	first_loop = false;
			// 	decrease_set.insert(v);
			// 	if(second_if_exist != -1) decrease_set.insert(second_if_exist);
			// } else {
				for (vid_t node : determined_set) {
					int deg_ge = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k && compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0) {
							deg_ge ++;
						}
					}
					if (compare_pvalue({deg_ge, degree[node]}, now_p) < 0 || deg_ge < k) {
						decrease_set.insert(node);
					}
				}
				for (vid_t node : linked_node) {
					int deg_ge = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k && compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0) {
							deg_ge ++;
						}
					}
					if (compare_pvalue({deg_ge, degree[node]}, now_p) < 0 || deg_ge < k) {
						decrease_set.insert(node);
					}
				}
			// }

			// printf("decrease_before %d\n", decrease_set.size());

			// puts("Step1 done");
			get_decrease_set(k, now_p, decrease_set);
			// printf("decrease : ");
			// for (vid_t node : decrease_set) {
			// 	printf("%d ", node);
			// }
			// puts("");
			// puts("Step2.1 done");

			for (vid_t node : decrease_set) {
				if (vc.find(node) == vc.end()) {
					node_kpcore_info[node][k - 1] = {-1, 1};
				}
			}

			int insert_num = 0;

			for (vid_t node : decrease_set) {
				if (vc.find(node) == vc.end()) {

					++insert_num;
					int deg_ge = 0;
					int deg_vc = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k) {
							// printf("%d->%d case1\n", node, nei);
							if ((compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0 || vc.find(nei) != vc.end() || decrease_set.find(nei) != decrease_set.end())) {
								deg_ge ++;
							} else {
								change_tag[nei] ++;
								// printf("change tag + : %d -> %d(%d)\n",node, nei, change_tag[nei]);
							}
							
						}
						// if (vc.find(nei) != vc.end()) {
						// 	// printf("%d->%d case2\n", node, nei);
						// 	deg_vc ++;
						// 	// deg_plus_vc[nei]++;
						// }
					}
					// deg_plus[node] = deg_ge;
					deg_plus_vc[node] = deg_ge;
					vc.insert(node);
					// qu.insert({{deg_plus[node], degree[node]}, node});

					if (determined_set.find(node) != determined_set.end()) {
						determined_set.erase(node);
					} 
				}
			}

			// printf("insert num %d\n", insert_num);
			

			if (!determined_set.empty()) {

				// maybe here too slow

				for (vid_t node : determined_set) {
					for (vid_t nei : neighbor[node]) {
						if (change_tag[nei] != 0) {
							change_tag[nei] --;
						}
					}
				}

				unordered_map<vid_t, p_t> remain;
				Oreverse.clear();
				p_t last_p = now_p;

				if (determined_set.size() <= 1000) {
					for (vid_t node : determined_set) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k && (compare_pvalue(node_kpcore_info[nei][k - 1], now_p) > 0 || determined_set.find(nei) != determined_set.end())) {
								deg_ge++;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
					}

					while (!remain.empty()) {
						bool p_lost = false;
						pair<vid_t, p_t> choose = {0, {1, 1}};

						for (auto entry : remain) {
							if (entry.second.first < k || compare_pvalue(entry.second, last_p) <= 0) {
								p_lost = true;
								choose = entry;
								break;
							}
						}

						if (!p_lost) {
							for (auto entry : remain) {
								if (compare_pvalue(entry.second, choose.second) < 0) {
									choose = entry;
								}
							}
						}

						for (vid_t nei : neighbor[choose.first]) {
							if (remain.find(nei) != remain.end()) {
								remain[nei].first --;
							}
						}

						if (p_lost) {
							node_kpcore_info[choose.first][k - 1] = last_p;
						} else {
							node_kpcore_info[choose.first][k - 1] = choose.second;
							last_p = choose.second;
						}

						remain.erase(choose.first);
						Oreverse.push_back(choose.first);
					}
				}

				else {
					min_heap mheap;
					mheap.node2index.clear();
					mheap.node2index.resize(num_node, -1);

					std::queue<std::pair<p_t, int> > not_satisfy_k;

					for (vid_t node : determined_set) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k && (compare_pvalue(node_kpcore_info[nei][k - 1], now_p) > 0 || determined_set.find(nei) != determined_set.end())) {
								deg_ge++;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
						if (deg_ge < k) not_satisfy_k.push({{deg_ge, degree[node]}, node});
						else mheap.insert({{deg_ge, degree[node]}, node});
					}

					// for (vid_t i = 0; i < num_node; i++) {
					// 	if (is_removed_in_pcore[i]) continue;
					// 	p_t t = make_pair(deg_in_pcore[i], degree[i]);
					// 	mheap.insert({ t,i });
					// }

					
					
					bool notSatisfy = false;
					p_t last_p = now_p;
					// std::vector<vid_t> naive_index_vec;
					// naive_index_vec.clear();
					while((!mheap.empty()) || (!not_satisfy_k.empty())) {
						
						std::pair<p_t, vid_t> top_ele;

						if(!not_satisfy_k.empty()) {
							top_ele = not_satisfy_k.front();
							top_ele.first = last_p;
							not_satisfy_k.pop();
							notSatisfy = true;
						}
						else {
							notSatisfy = false;
							top_ele = mheap.top();
							mheap.pop(); // ?
						}

						int now_id = top_ele.second;
						remain.erase(now_id);

						if(compare_pvalue(last_p, top_ele.first) == -1 && notSatisfy == false) {
							last_p = top_ele.first;
						}
						node_kpcore_info[now_id][k - 1] = last_p;
						Oreverse.push_back(now_id);

						for(auto x : neighbor[now_id]) {
							if(remain.find(x) == remain.end()) continue;
							else {
								remain[x].first --;
								if(remain[x].first == k - 1) {
									not_satisfy_k.push(make_pair(make_pair(0, 0), x));
									mheap.remove(mheap.node2index[x]);
								}
								else if(remain[x].first >= k){
									mheap.decrease_key(make_pair(remain[x].first, degree[x]), mheap.node2index[x]);
								}
								
								
							}
						}
						
					}
				}


				for (int i = Oreverse.size() - 1; i >= 0; i--) {
					decreased_order.push_back(Oreverse[i]);
				}
				

				// for (int i = Oreverse.size() - 1; i >= 0; i--){
				// 	new_order[k][++determined_rank] = Oreverse[i];
				// 	rank[Oreverse[i]][k] = determined_rank;
				// }
				
			}


		}

		for (int now_rank = begin_rank; now_rank < end_rank; now_rank ++) {
			if (vc.find(new_order[k][now_rank]) == vc.end()) {
				remain_order.push_back(new_order[k][now_rank]);
			}
		}

		// puts("HA");
		// print_new_order(k);

		if (end_rank == rank_size) {

			// puts("end rank = rank size");

			p_t next_p = {0, 1};

			// decompose
			if (!vc.empty()) {
				unordered_map<vid_t, p_t> remain;
				Oreverse.clear();
				p_t last_p = next_p;


				if (vc.size() <= 1000) {
					for (vid_t node : vc) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k) {
								deg_ge++;
								change_tag[nei] = 0;
							}
						}
						// printf("insert %d(%d/%d)\n", node, deg_ge, degree[node]);
						remain.insert({node, {deg_ge, degree[node]}});
					}

					while (!remain.empty()) {
						bool p_lost = false;
						pair<vid_t, p_t> choose = {0, {1, 1}};

						for (auto entry : remain) {
							if (entry.second.first < k || compare_pvalue(entry.second, last_p) <= 0) {
								p_lost = true;
								choose = entry;
								break;
							}
						}

						if (!p_lost) {
							for (auto entry : remain) {
								if (compare_pvalue(entry.second, choose.second) < 0) {
									choose = entry;
								}
							}
						}

						for (vid_t nei : neighbor[choose.first]) {
							if (remain.find(nei) != remain.end()) {
								remain[nei].first --;
							}
						}

						if (p_lost) {
							node_kpcore_info[choose.first][k - 1] = last_p;
						} else {
							node_kpcore_info[choose.first][k - 1] = choose.second;
							last_p = choose.second;
						}

						// printf("choose %d\n", choose.first);

						remain.erase(choose.first);
						Oreverse.push_back(choose.first);
					}
				}

				else {
					min_heap mheap;
					mheap.node2index.clear();
					mheap.node2index.resize(num_node, -1);

					std::queue<std::pair<p_t, int> > not_satisfy_k;

					for (vid_t node : vc) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k) {
								deg_ge++;
								change_tag[nei] = 0;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
						if (deg_ge < k) not_satisfy_k.push({{deg_ge, degree[node]}, node});
						else mheap.insert({{deg_ge, degree[node]}, node});
					}

					// for (vid_t i = 0; i < num_node; i++) {
					// 	if (is_removed_in_pcore[i]) continue;
					// 	p_t t = make_pair(deg_in_pcore[i], degree[i]);
					// 	mheap.insert({ t,i });
					// }

					
					
					bool notSatisfy = false;
					// std::vector<vid_t> naive_index_vec;
					// naive_index_vec.clear();
					while((!mheap.empty()) || (!not_satisfy_k.empty())) {
						
						std::pair<p_t, vid_t> top_ele;

						if(!not_satisfy_k.empty()) {
							top_ele = not_satisfy_k.front();
							top_ele.first = last_p;
							not_satisfy_k.pop();
							notSatisfy = true;
						}
						else {
							notSatisfy = false;
							top_ele = mheap.top();
							mheap.pop(); // ?
						}

						int now_id = top_ele.second;
						remain.erase(now_id);

						if(compare_pvalue(last_p, top_ele.first) == -1 && notSatisfy == false) {
							last_p = top_ele.first;
						}
						node_kpcore_info[now_id][k - 1] = last_p;
						Oreverse.push_back(now_id);

						for(auto x : neighbor[now_id]) {
							if(remain.find(x) == remain.end()) continue;
							else {
								remain[x].first --;
								if(remain[x].first == k - 1) {
									not_satisfy_k.push(make_pair(make_pair(0, 0), x));
									mheap.remove(mheap.node2index[x]);
								}
								else if(remain[x].first >= k){
									mheap.decrease_key(make_pair(remain[x].first, degree[x]), mheap.node2index[x]);
								}
								
								
							}
						}
						
					}
				}

				for (int i = Oreverse.size() - 1; i >= 0; i--) {
					decreased_order.push_back(Oreverse[i]);
				}

				
			}
			

			vc.clear();


		}

		

		// puts("Step3 done");
		

	} while (!vc.empty());



	int it_1 = 0;
	int it_2 = 0;
	int dec_size = decreased_order.size();
	int rem_size = remain_order.size();
	// int end_determined = determined_rank;
	// determined_rank = last_begin_rank
	// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
	// 	it_2 ++;
	// }

	while (it_1 < dec_size || it_2 < rem_size) {
		
		// printf("compare %d %d\n", Oreverse[it_1], new_order[k][it_2]);

		if (it_1 >= dec_size) {
			// if (vc.find(new_order[k][it_2]) == vc.end()) {
				new_order[k][++determined_rank] = remain_order[it_2];
				rank[remain_order[it_2]][k] = determined_rank;
			// }
			it_2 ++;
			// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
			// 	it_2 ++;
			// }
		} else if (it_2 >= rem_size) {
			// if (vc.find(Oreverse[it_1]) == vc.end()) {
				new_order[k][++determined_rank] = decreased_order[it_1];
				rank[decreased_order[it_1]][k] = determined_rank;
			// }
			it_1 ++;
		} else if (compare_pvalue(node_kpcore_info[decreased_order[it_1]][k - 1], node_kpcore_info[remain_order[it_2]][k - 1]) >= 0) {
			// if (vc.find(Oreverse[it_1]) == vc.end()) {
				new_order[k][++determined_rank] = decreased_order[it_1];
				rank[decreased_order[it_1]][k] = determined_rank;
			// }
			it_1 ++;
		} else {
			if (vc.find(new_order[k][it_2]) == vc.end()) {
				new_order[k][++determined_rank] = remain_order[it_2];
				rank[remain_order[it_2]][k] = determined_rank;
			}
			it_2 ++;
			// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
			// 	it_2 ++;
			// }
		}
	}


}

vector<long long> Graph::pDecrease_multiple_vertices_vertex_statistics(int k, std::vector<vid_t> othernodes) {

	// please check v must decrease and do this 

	// printf("in decrease k=%d, othernodes.size=%d\n", k, othernodes.size());
	// print_new_order(k);


	long long checknum_nei = 0;
	long long checknum_dec = 0;
	long long total_checknum_dec = 0;
	long long checknum_distinct_vertices = 0;


	int v = othernodes[0];

	for (vid_t u : othernodes) {
		// printf("%d ", u);
		change_tag[u] = 1;
		if(rank[u][k] < rank[v][k]) {
			v = u;
		}
	}
	// puts("");



	unordered_set<vid_t> decrease_set;
	unordered_set<vid_t> vc;
	// unordered_map<vid_t, int> deg_plus;
	unordered_map<vid_t, int> deg_plus_vc;

	// decrease_set.insert(v);
	
	// max_heap qu;
	// qu.node2index.clear();
	// qu.node2index.resize(num_node + 10, -1);

	min_heap min_qu;
	min_qu.node2index.clear();
	min_qu.node2index.resize(num_node + 10, -1);

	int begin_rank = rank[v][k];
	int end_rank = rank[v][k];
	int determined_rank;
	while (begin_rank > 0 && 
			compare_pvalue(node_kpcore_info[new_order[k][begin_rank]][k - 1], 
			node_kpcore_info[new_order[k][begin_rank - 1]][k - 1]) == 0 ) {
			
			begin_rank --;
	}
	end_rank = begin_rank;
	determined_rank = begin_rank - 1;
	int last_begin_rank = begin_rank;

	// bool first_loop = true;

	// [begin_rank, end_rank)

	int rank_size = new_order[k].size();

	int times = 0;

	vector<vid_t> decreased_order;
	vector<vid_t> remain_order;



	do {

		// printf("end_rank %d(node = %d)\n", end_rank, new_order[k][end_rank]);

		vector<vid_t> Oreverse;

		begin_rank = end_rank;
		p_t now_p = node_kpcore_info[new_order[k][begin_rank]][k - 1];

		while (end_rank < rank_size && 
		  compare_pvalue(node_kpcore_info[new_order[k][end_rank]][k - 1], now_p) == 0) {
			end_rank ++;
		}
 
		// printf("begin_rank = %d, end_rank = %d\n", begin_rank, end_rank);

		unordered_set<vid_t> determined_set;

		vector<vid_t> linked_node;


		for (int now_rank = begin_rank; now_rank < end_rank; now_rank ++) {
			int node = new_order[k][now_rank];

			if (change_tag[node] > 0) {
				change_tag[node] = 0;
				linked_node.push_back(node);
				for (vid_t w : neighbor[node]) {
					if (vc.find(w) != vc.end()) {
						// deg_plus[w]++;
						deg_plus_vc[w]++;
					}
				}
			}

		}
		// printf("linked_node.size = %d nowp %lld/%lld\n", linked_node.size(), now_p.first, now_p.second);

		if (linked_node.size() != 0) {

			checknum_nei += linked_node.size();

			for (vid_t node : vc) {
				// printf("node %d\n", node);
				if(deg_plus_vc[node] >= k && compare_pvalue({deg_plus_vc[node], degree[node]}, now_p) >= 0) {
					// puts("qu");
					// qu.remove(qu.node2index[node]);
					// puts("dp");
					// deg_plus.erase(node);
					// puts("dpv");
					deg_plus_vc.erase(node);
					// puts("ds");
					determined_set.insert(node);
					// puts("AA");
				}
			}

			for (vid_t node : determined_set) {
				vc.erase(node);
				node_kpcore_info[node][k - 1] = now_p;
			}

			decrease_set.clear();

				for (vid_t node : determined_set) {
					int deg_ge = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k && compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0) {
							deg_ge ++;
						}
					}
					if (compare_pvalue({deg_ge, degree[node]}, now_p) < 0 || deg_ge < k) {
						decrease_set.insert(node);
					}
				}
				for (vid_t node : linked_node) {
					int deg_ge = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k && compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0) {
							deg_ge ++;
						}
					}
					if (compare_pvalue({deg_ge, degree[node]}, now_p) < 0 || deg_ge < k) {
						decrease_set.insert(node);
					}
				}

			int temp_value = get_decrease_set_vertex_statistics(k, now_p, decrease_set);
			checknum_nei += temp_value;
			total_checknum_dec += decrease_set.size();


			for (int now_rank = begin_rank; now_rank < end_rank; now_rank++) {
				if (decrease_set.find(new_order[k][now_rank]) != decrease_set.end()) {
					checknum_distinct_vertices ++;
				}
			}			


			for (vid_t node : decrease_set) {
				if (vc.find(node) == vc.end()) {
					node_kpcore_info[node][k - 1] = {-1, 1};
				}
			}

			int insert_num = 0;

			for (vid_t node : decrease_set) {
				if (vc.find(node) == vc.end()) {

					++insert_num;
					int deg_ge = 0;
					int deg_vc = 0;
					for (vid_t nei : neighbor[node]) {
						if (node_kpcore_info[nei].size() >= k) {
							// printf("%d->%d case1\n", node, nei);
							if ((compare_pvalue(node_kpcore_info[nei][k - 1], now_p) >= 0 || vc.find(nei) != vc.end() || decrease_set.find(nei) != decrease_set.end())) {
								deg_ge ++;
							} else {
								change_tag[nei] ++;
								// printf("change tag + : %d -> %d(%d)\n",node, nei, change_tag[nei]);
							}
							
						}
					}
					// deg_plus[node] = deg_ge;
					deg_plus_vc[node] = deg_ge;
					vc.insert(node);
					// qu.insert({{deg_plus[node], degree[node]}, node});

					if (determined_set.find(node) != determined_set.end()) {
						determined_set.erase(node);
					} 
				}
			}

			// printf("insert num %d\n", insert_num);
			

			if (!determined_set.empty()) {

				checknum_dec += determined_set.size();

				// maybe here too slow

				for (vid_t node : determined_set) {
					for (vid_t nei : neighbor[node]) {
						if (change_tag[nei] != 0) {
							change_tag[nei] --;
						}
					}
				}

				unordered_map<vid_t, p_t> remain;
				Oreverse.clear();
				p_t last_p = now_p;

				if (determined_set.size() <= 1000) {
					for (vid_t node : determined_set) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k && (compare_pvalue(node_kpcore_info[nei][k - 1], now_p) > 0 || determined_set.find(nei) != determined_set.end())) {
								deg_ge++;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
					}

					while (!remain.empty()) {
						bool p_lost = false;
						pair<vid_t, p_t> choose = {0, {1, 1}};

						for (auto entry : remain) {
							if (entry.second.first < k || compare_pvalue(entry.second, last_p) <= 0) {
								p_lost = true;
								choose = entry;
								break;
							}
						}

						if (!p_lost) {
							for (auto entry : remain) {
								if (compare_pvalue(entry.second, choose.second) < 0) {
									choose = entry;
								}
							}
						}

						for (vid_t nei : neighbor[choose.first]) {
							if (remain.find(nei) != remain.end()) {
								remain[nei].first --;
							}
						}

						if (p_lost) {
							node_kpcore_info[choose.first][k - 1] = last_p;
						} else {
							node_kpcore_info[choose.first][k - 1] = choose.second;
							last_p = choose.second;
						}

						remain.erase(choose.first);
						Oreverse.push_back(choose.first);
					}
				}

				else {
					min_heap mheap;
					mheap.node2index.clear();
					mheap.node2index.resize(num_node, -1);

					std::queue<std::pair<p_t, int> > not_satisfy_k;

					for (vid_t node : determined_set) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k && (compare_pvalue(node_kpcore_info[nei][k - 1], now_p) > 0 || determined_set.find(nei) != determined_set.end())) {
								deg_ge++;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
						if (deg_ge < k) not_satisfy_k.push({{deg_ge, degree[node]}, node});
						else mheap.insert({{deg_ge, degree[node]}, node});
					}

					// for (vid_t i = 0; i < num_node; i++) {
					// 	if (is_removed_in_pcore[i]) continue;
					// 	p_t t = make_pair(deg_in_pcore[i], degree[i]);
					// 	mheap.insert({ t,i });
					// }

					
					
					bool notSatisfy = false;
					p_t last_p = now_p;
					// std::vector<vid_t> naive_index_vec;
					// naive_index_vec.clear();
					while((!mheap.empty()) || (!not_satisfy_k.empty())) {
						
						std::pair<p_t, vid_t> top_ele;

						if(!not_satisfy_k.empty()) {
							top_ele = not_satisfy_k.front();
							top_ele.first = last_p;
							not_satisfy_k.pop();
							notSatisfy = true;
						}
						else {
							notSatisfy = false;
							top_ele = mheap.top();
							mheap.pop(); // ?
						}

						int now_id = top_ele.second;
						remain.erase(now_id);

						if(compare_pvalue(last_p, top_ele.first) == -1 && notSatisfy == false) {
							last_p = top_ele.first;
						}
						node_kpcore_info[now_id][k - 1] = last_p;
						Oreverse.push_back(now_id);

						for(auto x : neighbor[now_id]) {
							if(remain.find(x) == remain.end()) continue;
							else {
								remain[x].first --;
								if(remain[x].first == k - 1) {
									not_satisfy_k.push(make_pair(make_pair(0, 0), x));
									mheap.remove(mheap.node2index[x]);
								}
								else if(remain[x].first >= k){
									mheap.decrease_key(make_pair(remain[x].first, degree[x]), mheap.node2index[x]);
								}
								
								
							}
						}
						
					}
				}


				for (int i = Oreverse.size() - 1; i >= 0; i--) {
					decreased_order.push_back(Oreverse[i]);
				}
				
				
			}


		}

		// add 0706
		for (int now_rank = begin_rank; now_rank < end_rank; now_rank ++) {
			if (vc.find(new_order[k][now_rank]) == vc.end()) {
				remain_order.push_back(new_order[k][now_rank]);
			}
		}

		if (end_rank == rank_size) {

			// puts("end rank = rank size");

			p_t next_p = {0, 1};

			checknum_dec += vc.size();
			// checknum_distinct_vertices += vc.size();

			// decompose
			if (!vc.empty()) {
				unordered_map<vid_t, p_t> remain;
				Oreverse.clear();
				p_t last_p = next_p;


				if (vc.size() <= 1000) {
					for (vid_t node : vc) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k) {
								deg_ge++;
								change_tag[nei] = 0;
							}
						}
						// printf("insert %d(%d/%d)\n", node, deg_ge, degree[node]);
						remain.insert({node, {deg_ge, degree[node]}});
					}

					while (!remain.empty()) {
						bool p_lost = false;
						pair<vid_t, p_t> choose = {0, {1, 1}};

						for (auto entry : remain) {
							if (entry.second.first < k || compare_pvalue(entry.second, last_p) <= 0) {
								p_lost = true;
								choose = entry;
								break;
							}
						}

						if (!p_lost) {
							for (auto entry : remain) {
								if (compare_pvalue(entry.second, choose.second) < 0) {
									choose = entry;
								}
							}
						}

						for (vid_t nei : neighbor[choose.first]) {
							if (remain.find(nei) != remain.end()) {
								remain[nei].first --;
							}
						}

						if (p_lost) {
							node_kpcore_info[choose.first][k - 1] = last_p;
						} else {
							node_kpcore_info[choose.first][k - 1] = choose.second;
							last_p = choose.second;
						}

						// printf("choose %d\n", choose.first);

						remain.erase(choose.first);
						Oreverse.push_back(choose.first);
					}
				}

				else {
					min_heap mheap;
					mheap.node2index.clear();
					mheap.node2index.resize(num_node, -1);

					std::queue<std::pair<p_t, int> > not_satisfy_k;

					for (vid_t node : vc) {
						int deg_ge = 0;
						for (vid_t nei : neighbor[node]) {
							if (node_kpcore_info[nei].size() >= k) {
								deg_ge++;
								change_tag[nei] = 0;
							}
						}
						remain.insert({node, {deg_ge, degree[node]}});
						if (deg_ge < k) not_satisfy_k.push({{deg_ge, degree[node]}, node});
						else mheap.insert({{deg_ge, degree[node]}, node});
					}

					
					bool notSatisfy = false;
					// std::vector<vid_t> naive_index_vec;
					// naive_index_vec.clear();
					while((!mheap.empty()) || (!not_satisfy_k.empty())) {
						
						std::pair<p_t, vid_t> top_ele;

						if(!not_satisfy_k.empty()) {
							top_ele = not_satisfy_k.front();
							top_ele.first = last_p;
							not_satisfy_k.pop();
							notSatisfy = true;
						}
						else {
							notSatisfy = false;
							top_ele = mheap.top();
							mheap.pop(); // ?
						}

						int now_id = top_ele.second;
						remain.erase(now_id);

						if(compare_pvalue(last_p, top_ele.first) == -1 && notSatisfy == false) {
							last_p = top_ele.first;
						}
						node_kpcore_info[now_id][k - 1] = last_p;
						Oreverse.push_back(now_id);

						for(auto x : neighbor[now_id]) {
							if(remain.find(x) == remain.end()) continue;
							else {
								remain[x].first --;
								if(remain[x].first == k - 1) {
									not_satisfy_k.push(make_pair(make_pair(0, 0), x));
									mheap.remove(mheap.node2index[x]);
								}
								else if(remain[x].first >= k){
									mheap.decrease_key(make_pair(remain[x].first, degree[x]), mheap.node2index[x]);
								}
								
								
							}
						}
						
					}
				}

				for (int i = Oreverse.size() - 1; i >= 0; i--) {
					decreased_order.push_back(Oreverse[i]);
				}

				
			}

			vc.clear();
		}

		// puts("Step3 done");
		

	} while (!vc.empty());


	int it_1 = 0;
	int it_2 = 0;
	int dec_size = decreased_order.size();
	int rem_size = remain_order.size();
	// int end_determined = determined_rank;
	// determined_rank = last_begin_rank
	// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
	// 	it_2 ++;
	// }

	while (it_1 < dec_size || it_2 < rem_size) {
		
		// printf("compare %d %d\n", Oreverse[it_1], new_order[k][it_2]);

		if (it_1 >= dec_size) {
			// if (vc.find(new_order[k][it_2]) == vc.end()) {
				new_order[k][++determined_rank] = remain_order[it_2];
				rank[remain_order[it_2]][k] = determined_rank;
			// }
			it_2 ++;
			// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
			// 	it_2 ++;
			// }
		} else if (it_2 >= rem_size) {
			// if (vc.find(Oreverse[it_1]) == vc.end()) {
				new_order[k][++determined_rank] = decreased_order[it_1];
				rank[decreased_order[it_1]][k] = determined_rank;
			// }
			it_1 ++;
		} else if (compare_pvalue(node_kpcore_info[decreased_order[it_1]][k - 1], node_kpcore_info[remain_order[it_2]][k - 1]) >= 0) {
			// if (vc.find(Oreverse[it_1]) == vc.end()) {
				new_order[k][++determined_rank] = decreased_order[it_1];
				rank[decreased_order[it_1]][k] = determined_rank;
			// }
			it_1 ++;
		} else {
			if (vc.find(new_order[k][it_2]) == vc.end()) {
				new_order[k][++determined_rank] = remain_order[it_2];
				rank[remain_order[it_2]][k] = determined_rank;
			}
			it_2 ++;
			// while(it_2 < end_rank && vc.find(new_order[k][it_2]) != vc.end()) {
			// 	it_2 ++;
			// }
		}
	}

	return {checknum_dec, checknum_nei, checknum_distinct_vertices, total_checknum_dec};


}

vector<long long> Graph::local_delete_maintain_p_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	std::vector<long long> tempvec;

	vid_t w_s;

	if(node_kpcore_info[v].size() < k) return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_increase_delete(k, v) == false) {
				// puts("HA Con1");
				return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
			}
			// puts("in");
			tempvec = pIncrease_vertex_statistics(k, rank[v][k]);
			checknum_nei_inc += tempvec[1];
			checknum_dec_inc += tempvec[0];
			inc_times ++;
			
			return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
		}
	}

	// puts("step1 done");

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}

	// p_t p_v = node_kpcore_info[v][k - 1];
	// if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	// else w_s = u;
	

	if(compare_pvalue(node_kpcore_info[v][k - 1], node_kpcore_info[u][k - 1]) == 0) {
		vector<vid_t> changed_node = {};
		if(check_decrease(k, u, v)) {
			changed_node.push_back(u);
		}
		if (check_decrease(k, v, u)) {
			changed_node.push_back(v);
		}
		if(changed_node.size() != 0) {
			tempvec = pDecrease_multiple_vertices_vertex_statistics(k, changed_node);
			checknum_dec_dec += tempvec[0];
			checknum_nei_dec += tempvec[1];
			dec_times ++;
			distinct_check_num_dec += tempvec[2];
			total_check_number_dec += tempvec[3];
		}
	}
	else if (check_decrease(k, u, v) == true){
		tempvec = pDecrease_multiple_vertices_vertex_statistics(k, {u});
		checknum_dec_dec += tempvec[0];
		checknum_nei_dec += tempvec[1];
		dec_times ++;
		distinct_check_num_dec += tempvec[2];
		total_check_number_dec += tempvec[3];

	}
	
	// puts("step2 done");

		//u v tongshi decrease ? 
	if(check_increase_delete(k, v) == true){
		tempvec = pIncrease_vertex_statistics(k, rank[v][k]);
		checknum_nei_inc += tempvec[1];
		checknum_dec_inc += tempvec[0];
		inc_times ++;
	}
	

	

	return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
}

void Graph::local_delete_maintain_p_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	vid_t w_s;

	if(node_kpcore_info[v].size() < k) return;
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_increase_delete(k, v) == false) {
				// puts("HA Con1");
				return;
			}
			// puts("in");
			pIncrease_woc(k, rank[v][k]);
			
			return;
		}
	}

	// puts("step1 done");

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}

	// p_t p_v = node_kpcore_info[v][k - 1];
	// if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	// else w_s = u;
	

	if(compare_pvalue(node_kpcore_info[v][k - 1], node_kpcore_info[u][k - 1]) == 0) {
		vector<vid_t> changed_node = {};
		if(check_decrease(k, u, v)) {
			changed_node.push_back(u);
		}
		if (check_decrease(k, v, u)) {
			changed_node.push_back(v);
		}
		if(changed_node.size() != 0) {
			pDecrease_multiple_vertices(k, changed_node);
		}
	}
	else if (check_decrease(k, u, v) == true)
		pDecrease_multiple_vertices(k, {u});
	
	// puts("step2 done");

		//u v tongshi decrease ? 
	if(check_increase_delete(k, v) == true)
		pIncrease_woc(k, rank[v][k]);
	


	

	return;
}


void Graph::local_delete_maintain_p(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	vid_t w_s;

	if(node_kpcore_info[v].size() < k) return;
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_increase_delete(k, v) == false) {
				// puts("HA Con1");
				return;
			}
			// puts("in");
			pIncrease(k, rank[v][k]);
			
			return;
		}
	}

	// puts("step1 done");

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}

	// p_t p_v = node_kpcore_info[v][k - 1];
	// if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	// else w_s = u;
	

	if(compare_pvalue(node_kpcore_info[v][k - 1], node_kpcore_info[u][k - 1]) == 0) {
		vector<vid_t> changed_node = {};
		if(check_decrease(k, u, v)) {
			changed_node.push_back(u);
		}
		if (check_decrease(k, v, u)) {
			changed_node.push_back(v);
		}
		if(changed_node.size() != 0) {
			pDecrease_multiple_vertices(k, changed_node);
		}
	}
	else if (check_decrease(k, u, v) == true)
		pDecrease_multiple_vertices(k, {u});
	
	// puts("step2 done");

		//u v tongshi decrease ? 
	if(check_increase_delete(k, v) == true)
		pIncrease(k, rank[v][k]);
	


	

	return;
}

void Graph::local_delete_maintain_p_for_changed_core_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	bool containv = false;
	bool containu = false;
	for (vid_t x : changed_nodes) {
		if(x == v) {
			containv = true;
		}
		if(x == u) {
			containu = true;
		}
	}

	if(containu == false && containv == true) {
		std::swap(u, v);
		std::swap(containu, containv);
	}

	vid_t w_s;
	bool add = false;
	if(containv == false && check_decrease(k, v, u) == true) {
		// puts("HA");
		add = true;
		changed_nodes.push_back(v);
	}
	pDecrease_multiple_vertices(k, changed_nodes);
	if(add) {
		changed_nodes.pop_back();
	}
	for (vid_t x : changed_nodes) {
		new_order[k].pop_back();
		node_kpcore_info[x].pop_back();
		rank[x].pop_back();
	}
	// print_new_order(k);
	// only u decreases
	if(containv == false) {
		if(check_increase_delete(k, v) == true) {
			pIncrease_woc(k, rank[v][k]);
		}
	}



	return;
}

void Graph::local_delete_maintain_p_for_changed_core(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	bool containv = false;
	bool containu = false;
	for (vid_t x : changed_nodes) {
		if(x == v) {
			containv = true;
		}
		if(x == u) {
			containu = true;
		}
	}

	if(containu == false && containv == true) {
		std::swap(u, v);
		std::swap(containu, containv);
	}

	vid_t w_s;
	bool add = false;
	if(containv == false && check_decrease(k, v, u) == true) {
		// puts("HA");
		add = true;
		changed_nodes.push_back(v);
	}
	pDecrease_multiple_vertices(k, changed_nodes);
	if(add) {
		changed_nodes.pop_back();
	}
	for (vid_t x : changed_nodes) {
		new_order[k].pop_back();
		node_kpcore_info[x].pop_back();
		rank[x].pop_back();
	}
	// print_new_order(k);
	// only u decreases
	if(containv == false) {
		if(check_increase_delete(k, v) == true) {
			pIncrease(k, rank[v][k]);
		}
	}



	return;
}

vector<long long> Graph::local_delete_maintain_p_for_changed_core_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, std::vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.
	// for the changed_nodes, just call it when k = min_k + 1
 
	// printf("k=%d, insert(%d,%d)\n",k ,u, v);
	// print_new_order(k);

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	std::vector<long long> tempvec;

	bool containv = false;
	bool containu = false;
	for (vid_t x : changed_nodes) {
		if(x == v) {
			containv = true;
		}
		if(x == u) {
			containu = true;
		}
	}

	if(containu == false && containv == true) {
		std::swap(u, v);
		std::swap(containu, containv);
	}

	vid_t w_s;
	bool add = false;
	if(containv == false && check_decrease(k, v, u) == true) {
		// puts("HA");
		add = true;
		changed_nodes.push_back(v);
	}
	tempvec = pDecrease_multiple_vertices_vertex_statistics(k, changed_nodes);
	checknum_dec_dec += tempvec[0];
	checknum_nei_dec += tempvec[1];
	dec_times ++;
	distinct_check_num_dec += tempvec[2];
	total_check_number_dec += tempvec[3];
	if(add) {
		changed_nodes.pop_back();
	}
	for (vid_t x : changed_nodes) {
		new_order[k].pop_back();
		node_kpcore_info[x].pop_back();
		rank[x].pop_back();
	}
	// print_new_order(k);
	// only u decreases
	if(containv == false) {
		if(check_increase_delete(k, v) == true) {
			tempvec = pIncrease_vertex_statistics(k, rank[v][k]);
			checknum_nei_inc += tempvec[1];
			checknum_dec_inc += tempvec[0];
			inc_times ++;
		}
	}



	return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
}

void Graph::local_insert_maintain_p(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.


	vid_t w_s;


	if(node_kpcore_info[v].size() < k) return;
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_decrease(k, v, u) == false) {
				// puts("HA Con1");
				return;
			}
			// puts("in");
			pDecrease_multiple_vertices(k, {v});
			return;
		}
	}

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}


	p_t p_v = node_kpcore_info[v][k - 1];
	if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	else w_s = u;
	pIncrease(k, rank[w_s][k]);
	if(compare_pvalue(p_v, node_kpcore_info[v][k - 1]) == 0) {
		if(check_decrease(k, v, u) == true) {
			// puts("in");
			pDecrease_multiple_vertices(k, {v});
			// puts("out");
			// w_s = get_w_start_for_V2_1_new_order(k, u, v, 1, node_kpcore_info[v][k - 1]);
			// update_p_value_baseline_given_w_start_V2_1_new_order(k, w_s, v);
		}
	}
	

	return;
}

void Graph::local_insert_maintain_p_woc(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.


	vid_t w_s;


	if(node_kpcore_info[v].size() < k) return;
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_decrease(k, v, u) == false) {
				// puts("HA Con1");
				return;
			}
			// puts("in");
			pDecrease_multiple_vertices(k, {v});
			return;
		}
	}

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}


	p_t p_v = node_kpcore_info[v][k - 1];
	if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	else w_s = u;
	pIncrease_woc(k, rank[w_s][k]);
	if(compare_pvalue(p_v, node_kpcore_info[v][k - 1]) == 0) {
		if(check_decrease(k, v, u) == true) {
			// puts("in");
			pDecrease_multiple_vertices(k, {v});
			// puts("out");
			// w_s = get_w_start_for_V2_1_new_order(k, u, v, 1, node_kpcore_info[v][k - 1]);
			// update_p_value_baseline_given_w_start_V2_1_new_order(k, w_s, v);
		}
	}
	

	return;
}

std::vector<long long> Graph::local_insert_maintain_p_vertex_statistics(int k, vid_t u, vid_t v, int pre_k_u, int pre_k_v, vector<vid_t> changed_nodes) {
	// check: now_k_v >= now_k_u must hold here.

	vid_t w_s;

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	std::vector<long long> tempvec;

	if(node_kpcore_info[v].size() < k) return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
	else {
		++ maintain_cnt;
		if(node_kpcore_info[u].size() < k) {
			// p_order[k][rank[v][k]].first.second ++;

			if(check_decrease(k, v, u) == false) {
				// puts("HA Con1");
				return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
			}
			// puts("in");
			tempvec = pDecrease_multiple_vertices_vertex_statistics(k, {v});
			checknum_dec_dec += tempvec[0];
			checknum_nei_dec += tempvec[1];
			dec_times ++;
			distinct_check_num_dec += tempvec[2];
			total_check_number_dec += tempvec[3];

			return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
		}
	}

	if(rank[u][k] < rank[v][k]) {
		std::swap(u, v);
		std::swap(pre_k_u, pre_k_v);
	}


	p_t p_v = node_kpcore_info[v][k - 1];
	if(node_kpcore_info[u][k - 1].first == 0) w_s = new_order[k][new_order[k].size() - 1];
	else w_s = u;
	tempvec = pIncrease_vertex_statistics(k, rank[w_s][k]);
	checknum_nei_inc += tempvec[1];
	checknum_dec_inc += tempvec[0];
	inc_times ++;
	if(compare_pvalue(p_v, node_kpcore_info[v][k - 1]) == 0) {
		if(check_decrease(k, v, u) == true) {
			// puts("in");
			tempvec = pDecrease_multiple_vertices_vertex_statistics(k, {v});
			checknum_dec_dec += tempvec[0];
			checknum_nei_dec += tempvec[1];
			dec_times ++;
			distinct_check_num_dec += tempvec[2];
			total_check_number_dec += tempvec[3];
			// puts("out");
			// w_s = get_w_start_for_V2_1_new_order(k, u, v, 1, node_kpcore_info[v][k - 1]);
			// update_p_value_baseline_given_w_start_V2_1_new_order(k, w_s, v);
		}
	}


	return {dec_times, checknum_nei_dec, checknum_dec_dec, total_check_number_dec, distinct_check_num_dec, inc_times, checknum_nei_inc, checknum_dec_inc};
}


double Graph::global_insert_edge(vid_t u, vid_t v) {

	auto start = chrono::system_clock::now();
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size();
	int k_min = MIN(k_u, k_v);
	insert_edge(u, v);


	// accelarate subgraph computation
	vector<bool> is_already_added_into_subgraph;
	is_already_added_into_subgraph.resize(num_node);
	// re-compute k_min+1-core
	// method : directly decomposition to the k_min+1 core.
	vector<vid_t> s;
	vector<bool> is_removed;
	is_removed.resize(num_node);
	fill_n(is_removed.begin(), is_removed.size(), false);
	auto tmp_deg = degree;
	for (int i = 0; i < num_node; i++) {
		if (tmp_deg[i] < k_min + 1) {
			s.push_back(i);
		}
	}
	while (!s.empty()) {
		vid_t cur = s.back();
		s.pop_back();
		is_removed[cur] = true;
		for (auto neigh : neighbor[cur]) {
			if (is_removed[neigh]) continue;
			tmp_deg[neigh]--;
			if (tmp_deg[neigh] == k_min) {
				s.push_back(neigh);
			}
		}
	}
	bool kcore_change = false;
	for (int i = 0; i < num_node; i++) {
		if (!is_removed[i] && node_kpcore_info[i].size() == k_min) {
			p_t tt = { 0,1 };
			node_kpcore_info[i].push_back(tt);
			kcore_change = true;
		}
	}
	//end of re-compute k_min+1-core
	

	int max_k = k_min; // kmin = min(ku, kv)
	if (kcore_change) {
		max_k = k_min + 1;
	}
	for (int k = 1; k <= max_k; k++) {


		// if(check_increase(k, u, v) == false && 
		// 	check_decrease(k, u, v) == false &&
		// 	check_increase(k, v, u) == false &&
		// 	check_decrease(k, v, u) == false) {
		// 		continue;
		// }
		
		// calculate lower bound
		p_t l_b;
		if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) < 0) {
			l_b = node_kpcore_info[u][k - 1];
		}
		else {
			l_b = node_kpcore_info[v][k - 1];
		}
		// node_kpcore_info[a][b] : pair<long long, long long> : p value of node a in b-core.
		// compute upper bound for u
		p_t u_b_u;
		vector<p_t> neighbor_upper_bound;
		int en = 0;
		for (auto n : neighbor[u]) {
			int t = 0;
			if (node_kpcore_info[n].size() < k) {
				continue;
			}
			else {
				en++;
			}
			for (auto nn : neighbor[n]) {
				if (node_kpcore_info[nn].size() >= k) {
					t++;
				}
			}
			neighbor_upper_bound.push_back(make_pair(t, neighbor[n].size()));	//O(du log du)

		}
		sort(neighbor_upper_bound.begin(), neighbor_upper_bound.end(), comp_p_desc);
		p_t naive_upper_bound_u = make_pair(en, neighbor[u].size());
		if (compare_pvalue(naive_upper_bound_u, neighbor_upper_bound[k - 1]) >= 0) {
			u_b_u = neighbor_upper_bound[k - 1];
		}
		else {
			u_b_u = naive_upper_bound_u;
			//cout << "useless top down u" << endl;
		}
		// compute upper bound for v
		p_t u_b_v;
		neighbor_upper_bound.clear();
		en = 0;
		for (auto n : neighbor[v]) {
			int t = 0;
			if (node_kpcore_info[n].size() < k) {
				continue;
			}
			else {
				en++;
			}
			for (auto nn : neighbor[n]) {
				if (node_kpcore_info[nn].size() >= k) {
					t++;
				}
			}
			neighbor_upper_bound.push_back(make_pair(t, neighbor[n].size()));

		}
		sort(neighbor_upper_bound.begin(), neighbor_upper_bound.end(), comp_p_desc);	//O(dv log dv)
		p_t naive_upper_bound_v = make_pair(en, neighbor[v].size());	
		if (compare_pvalue(naive_upper_bound_v, neighbor_upper_bound[k - 1]) >= 0) {
			u_b_v = neighbor_upper_bound[k - 1];
		}
		else {
			u_b_v = naive_upper_bound_v;
			//cout << "useless top down v" << endl;
		}
		p_t u_b;
		if (compare_pvalue(u_b_u, u_b_v) < 0) {
			u_b = u_b_u;
		}
		else {
			u_b = u_b_v;
		}
		if (compare_pvalue(u_b, node_kpcore_info[u][k - 1]) < 0) {
			u_b = node_kpcore_info[u][k - 1];
		}
		else if (compare_pvalue(u_b, node_kpcore_info[v][k - 1]) < 0) {
			u_b = node_kpcore_info[v][k - 1];
		}

		//end of compute lower bound and upper bound 

		

		// compute subgraph
		//auto degree_in_p = degree;
		//vector<bool> is_removed;
		//is_removed.resize(num_node);
		//fill_n(is_removed.begin(), is_removed.size(), false);
		//for (int i = 0; i < num_node; i++) {
		//	if (node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k], l_b) < 0) {
		//		is_removed[i] = true;
		//	}
		//	else {
		//		continue;
		//	}
		//	for (auto n : neighbor[i]) {
		//		degree_in_p[n]--;
		//	}
		//}

		// compute subgraph rooted at u and v
		//vector<bool> is_removed;
		//is_removed.resize(num_node);
		fill_n(is_removed.begin(), is_removed.size(), true); // is_removed : the nodes that are not in k,l_b-core but not in k,u_b core
		fill_n(is_already_added_into_subgraph.begin(), is_already_added_into_subgraph.size(), false);
		unordered_set<vid_t> subgraph_nodes;
		unordered_map<vid_t, int> subgraph_degree;
		unordered_set<vid_t> fixed; // the nodes that is in k,u_b-core (don't need to change thus fixed)
		vector<vid_t> stack;

		for(int i = 0; i < num_node; i++) {
			if (node_kpcore_info[i].size() >= k) {
				if( compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else if( compare_pvalue(node_kpcore_info[i][k - 1], l_b) >= 0) {
					subgraph_nodes.insert(i);
					is_removed[i] = false;
				}
			}
		}
	
		// compute degree
		for (auto sub_n : subgraph_nodes) {
			for (auto mn : neighbor[sub_n]) {
				if (!is_removed[mn]) {
					subgraph_degree[sub_n]++;
				}
				else if (fixed.find(mn) != fixed.end()) {
					subgraph_degree[sub_n]++;
				}
			}
		}

		//check  
		//printf("%d %d ", k, subgraph_nodes.size());
		int mheapSize = 0;

		// re-decomposition kpcore between l_b and u_b
		min_heap mheap; // min_heap write by the author
		mheap.node2index.resize(num_node, -1);
		for (auto i : subgraph_nodes) {
			if (is_removed[i]) continue;
			p_t t = make_pair(subgraph_degree[i], degree[i]); // p_t : pair<long long, long long>
			mheap.insert({ t,i }); // t : p_number(represent by two long int), i : node_id

			//  
			mheapSize++;

		}

		//printf("%d\n", mheapSize);

		while (!mheap.empty()) { 
			// peak at the current min p value and get the kp-core.
			auto top_ele = mheap.pop();
			p_t p_min = top_ele.first;
			vector<vid_t> stack_pcore = { top_ele.second };
			// first, we let those in the heap with the p value min_p get into stack_pcore,
			// while we doing this, we didn't change any p value in the heap.
			while (!mheap.empty()) {
				auto te = mheap.top();
				if (compare_pvalue(p_min, te.first) < 0) {
					break;
				}
				else {
					mheap.pop();
					stack_pcore.push_back(te.second);
				}
			}
			int d = mcd(p_min.first, p_min.second);
			p_min = make_pair(p_min.first / d, p_min.second / d);
			//if (d == 0) {
			//	p_min = { 0,1 };
			//}
			//else {
			//	p_min = make_pair(p_min.first / d, p_min.second / d);
			//}
			if (stack_pcore.empty()) {
				break; // terminate
			}
			// for those in stack_pcore, we remove it in the heap and 
			// update its neighbour's p number, and remove those unqualified
			// and update mheap (using decrease_key, faster?)
			while (!stack_pcore.empty()) {
				vid_t cur = stack_pcore.back();
				stack_pcore.pop_back();
				if (is_removed[cur]) continue;
				mheap.remove(mheap.node2index[cur]); // delete twice? already deleted in line 897?
				// modify index and node kpcore info
				auto old_p = node_kpcore_info[cur][k - 1];
				node_kpcore_info[cur][k - 1] = p_min;

				is_removed[cur] = true;
				for (int i = 0; i < neighbor[cur].size(); i++) {
					//
					vid_t neigh = neighbor[cur][i];
					if (is_removed[neigh]) {
						continue;
					}
					subgraph_degree[neigh]--;
					if (subgraph_degree[neigh] < k) {
						stack_pcore.push_back(neigh);
					}
					else {
						p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
						if (compare_pvalue(t, p_min) <= 0) {
							stack_pcore.push_back(neigh);
						}
						else {
							mheap.decrease_key(t, mheap.node2index[neigh]);
						}
					}
				}
			}
		}
	}
	//process large k
	vid_t root; int l_k = 0;
	k_v = node_kpcore_info[v].size();
	k_u = node_kpcore_info[u].size();
	vid_t other_node;
	if (k_v > k_u) {
		root = v;
		l_k = k_v;
		other_node = u;
	}
	else if (k_v < k_u) {
		root = u;
		l_k = k_u;
		other_node = v;
	}
	else {
		l_k = 0;
		root = 0;
	}
	for (int k = max_k + 1; k <= l_k; k++) {
		
		// if(check_decrease(k, root, other_node) == false) {
		// 	continue;
		// }

		// compute lower bound for root
		p_t p = node_kpcore_info[root][k - 1];
		int l = 0;
		for (auto neigh : neighbor[root]) {
			if (node_kpcore_info[neigh].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) >= 0) {
				l++;
			}
		}
		p_t l_b = { l,degree[root] };
		// need to recompute
		if (compare_pvalue(p, l_b) > 0) {
			//int st = 0;
			//for (auto neigh : neighbor[root]) {
			//	if (node_kpcore_info[neigh].size() < k) continue;
			//	if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) > 0) {
			//		st++;
			//	}
			//}
			//l_b = { st,degree[root] };
			p_t u_b = p;
			fill_n(is_removed.begin(), is_removed.size(), true);
			fill_n(is_already_added_into_subgraph.begin(), is_already_added_into_subgraph.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;

			for(int i = 0; i < num_node; i++) {
				if (node_kpcore_info[i].size() >= k) {
					if( compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
						fixed.insert(i);
						is_removed[i] = true;
					}
					else if( compare_pvalue(node_kpcore_info[i][k - 1], l_b) >= 0) {
						subgraph_nodes.insert(i);
						is_removed[i] = false;
					}
				}
			}

		
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);

			//check  
			//printf("%d %lu\n", k, subgraph_nodes.size());

			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				//if (d == 0) {
				//	p_min = { 0,1 };
				//}
				//else {
				//	p_min = make_pair(p_min.first / d, p_min.second / d);
				//}
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	return elapsed_seconds.count();
}

long long Graph::global_insert_edge_vertex_statistics(vid_t u, vid_t v) {

	long long checknum_nei = 0;
	long long checknum_dec = 0;

	auto start = chrono::system_clock::now();
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size();
	int k_min = MIN(k_u, k_v);
	insert_edge(u, v);


	// accelarate subgraph computation
	vector<bool> is_already_added_into_subgraph;
	is_already_added_into_subgraph.resize(num_node);
	// re-compute k_min+1-core
	// method : directly decomposition to the k_min+1 core.
	vector<vid_t> s;
	vector<bool> is_removed;
	is_removed.resize(num_node);
	fill_n(is_removed.begin(), is_removed.size(), false);
	auto tmp_deg = degree;
	for (int i = 0; i < num_node; i++) {
		if (tmp_deg[i] < k_min + 1) {
			s.push_back(i);
		}
	}
	while (!s.empty()) {
		vid_t cur = s.back();
		s.pop_back();
		is_removed[cur] = true;
		for (auto neigh : neighbor[cur]) {
			if (is_removed[neigh]) continue;
			tmp_deg[neigh]--;
			if (tmp_deg[neigh] == k_min) {
				s.push_back(neigh);
			}
		}
	}
	bool kcore_change = false;
	for (int i = 0; i < num_node; i++) {
		if (!is_removed[i] && node_kpcore_info[i].size() == k_min) {
			p_t tt = { 0,1 };
			node_kpcore_info[i].push_back(tt);
			kcore_change = true;
		}
	}
	//end of re-compute k_min+1-core
	

	int max_k = k_min; // kmin = min(ku, kv)
	if (kcore_change) {
		max_k = k_min + 1;
	}
	for (int k = 1; k <= max_k; k++) {


		// if(check_increase(k, u, v) == false && 
		// 	check_decrease(k, u, v) == false &&
		// 	check_increase(k, v, u) == false &&
		// 	check_decrease(k, v, u) == false) {
		// 		continue;
		// }
		
		// calculate lower bound
		p_t l_b;
		if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) < 0) {
			l_b = node_kpcore_info[u][k - 1];
		}
		else {
			l_b = node_kpcore_info[v][k - 1];
		}
		// node_kpcore_info[a][b] : pair<long long, long long> : p value of node a in b-core.
		// compute upper bound for u
		p_t u_b_u;
		vector<p_t> neighbor_upper_bound;
		int en = 0;
		for (auto n : neighbor[u]) {
			int t = 0;
			if (node_kpcore_info[n].size() < k) {
				continue;
			}
			else {
				en++;
			}
			for (auto nn : neighbor[n]) {
				if (node_kpcore_info[nn].size() >= k) {
					t++;
				}
			}
			neighbor_upper_bound.push_back(make_pair(t, neighbor[n].size()));	//O(du log du)

		}
		sort(neighbor_upper_bound.begin(), neighbor_upper_bound.end(), comp_p_desc);
		p_t naive_upper_bound_u = make_pair(en, neighbor[u].size());
		if (compare_pvalue(naive_upper_bound_u, neighbor_upper_bound[k - 1]) >= 0) {
			u_b_u = neighbor_upper_bound[k - 1];
		}
		else {
			u_b_u = naive_upper_bound_u;
			//cout << "useless top down u" << endl;
		}
		// compute upper bound for v
		p_t u_b_v;
		neighbor_upper_bound.clear();
		en = 0;
		for (auto n : neighbor[v]) {
			int t = 0;
			if (node_kpcore_info[n].size() < k) {
				continue;
			}
			else {
				en++;
			}
			for (auto nn : neighbor[n]) {
				if (node_kpcore_info[nn].size() >= k) {
					t++;
				}
			}
			neighbor_upper_bound.push_back(make_pair(t, neighbor[n].size()));

		}
		sort(neighbor_upper_bound.begin(), neighbor_upper_bound.end(), comp_p_desc);	//O(dv log dv)
		p_t naive_upper_bound_v = make_pair(en, neighbor[v].size());	
		if (compare_pvalue(naive_upper_bound_v, neighbor_upper_bound[k - 1]) >= 0) {
			u_b_v = neighbor_upper_bound[k - 1];
		}
		else {
			u_b_v = naive_upper_bound_v;
			//cout << "useless top down v" << endl;
		}
		p_t u_b;
		if (compare_pvalue(u_b_u, u_b_v) < 0) {
			u_b = u_b_u;
		}
		else {
			u_b = u_b_v;
		}
		if (compare_pvalue(u_b, node_kpcore_info[u][k - 1]) < 0) {
			u_b = node_kpcore_info[u][k - 1];
		}
		else if (compare_pvalue(u_b, node_kpcore_info[v][k - 1]) < 0) {
			u_b = node_kpcore_info[v][k - 1];
		}

		//end of compute lower bound and upper bound 

		// compute subgraph rooted at u and v
		//vector<bool> is_removed;
		//is_removed.resize(num_node);
		fill_n(is_removed.begin(), is_removed.size(), true); // is_removed : the nodes that are not in k,l_b-core but not in k,u_b core
		fill_n(is_already_added_into_subgraph.begin(), is_already_added_into_subgraph.size(), false);
		unordered_set<vid_t> subgraph_nodes;
		unordered_map<vid_t, int> subgraph_degree;
		unordered_set<vid_t> fixed; // the nodes that is in k,u_b-core (don't need to change thus fixed)
		vector<vid_t> stack;

		for(int i = 0; i < num_node; i++) {
			if (node_kpcore_info[i].size() >= k) {
				if( compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else if( compare_pvalue(node_kpcore_info[i][k - 1], l_b) >= 0) {
					subgraph_nodes.insert(i);
					is_removed[i] = false;
				}
			}
		}
		
	
		// compute degree
		for (auto sub_n : subgraph_nodes) {
			for (auto mn : neighbor[sub_n]) {
				if (!is_removed[mn]) {
					subgraph_degree[sub_n]++;
				}
				else if (fixed.find(mn) != fixed.end()) {
					subgraph_degree[sub_n]++;
				}
			}
		}

		//check  
		//printf("%d %d ", k, subgraph_nodes.size());
		int mheapSize = 0;

		// re-decomposition kpcore between l_b and u_b
		min_heap mheap; // min_heap write by the author
		mheap.node2index.resize(num_node, -1);
		for (auto i : subgraph_nodes) {
			if (is_removed[i]) continue;
			p_t t = make_pair(subgraph_degree[i], degree[i]); // p_t : pair<long long, long long>
			mheap.insert({ t,i }); // t : p_number(represent by two long int), i : node_id

			//  
			mheapSize++;

		}

		checknum_dec += subgraph_nodes.size();
		checknum_nei += subgraph_nodes.size() + fixed.size();

		//printf("%d\n", mheapSize);

		while (!mheap.empty()) { 
			// peak at the current min p value and get the kp-core.
			auto top_ele = mheap.pop();
			p_t p_min = top_ele.first;
			vector<vid_t> stack_pcore = { top_ele.second };
			// first, we let those in the heap with the p value min_p get into stack_pcore,
			// while we doing this, we didn't change any p value in the heap.
			while (!mheap.empty()) {
				auto te = mheap.top();
				if (compare_pvalue(p_min, te.first) < 0) {
					break;
				}
				else {
					mheap.pop();
					stack_pcore.push_back(te.second);
				}
			}
			int d = mcd(p_min.first, p_min.second);
			p_min = make_pair(p_min.first / d, p_min.second / d);
			//if (d == 0) {
			//	p_min = { 0,1 };
			//}
			//else {
			//	p_min = make_pair(p_min.first / d, p_min.second / d);
			//}
			if (stack_pcore.empty()) {
				break; // terminate
			}
			// for those in stack_pcore, we remove it in the heap and 
			// update its neighbour's p number, and remove those unqualified
			// and update mheap (using decrease_key, faster?)
			while (!stack_pcore.empty()) {
				vid_t cur = stack_pcore.back();
				stack_pcore.pop_back();
				if (is_removed[cur]) continue;
				mheap.remove(mheap.node2index[cur]); // delete twice? already deleted in line 897?
				// modify index and node kpcore info
				auto old_p = node_kpcore_info[cur][k - 1];
				node_kpcore_info[cur][k - 1] = p_min;

				is_removed[cur] = true;
				for (int i = 0; i < neighbor[cur].size(); i++) {
					//
					vid_t neigh = neighbor[cur][i];
					if (is_removed[neigh]) {
						continue;
					}
					subgraph_degree[neigh]--;
					if (subgraph_degree[neigh] < k) {
						stack_pcore.push_back(neigh);
					}
					else {
						p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
						if (compare_pvalue(t, p_min) <= 0) {
							stack_pcore.push_back(neigh);
						}
						else {
							mheap.decrease_key(t, mheap.node2index[neigh]);
						}
					}
				}
			}
		}
	}
	//process large k
	vid_t root; int l_k = 0;
	k_v = node_kpcore_info[v].size();
	k_u = node_kpcore_info[u].size();
	vid_t other_node;
	if (k_v > k_u) {
		root = v;
		l_k = k_v;
		other_node = u;
	}
	else if (k_v < k_u) {
		root = u;
		l_k = k_u;
		other_node = v;
	}
	else {
		l_k = 0;
		root = 0;
	}
	for (int k = max_k + 1; k <= l_k; k++) {
		
		// if(check_decrease(k, root, other_node) == false) {
		// 	continue;
		// }

		// compute lower bound for root
		p_t p = node_kpcore_info[root][k - 1];
		int l = 0;
		for (auto neigh : neighbor[root]) {
			if (node_kpcore_info[neigh].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) >= 0) {
				l++;
			}
		}
		p_t l_b = { l,degree[root] };
		// need to recompute
		if (compare_pvalue(p, l_b) > 0) {
			//int st = 0;
			//for (auto neigh : neighbor[root]) {
			//	if (node_kpcore_info[neigh].size() < k) continue;
			//	if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) > 0) {
			//		st++;
			//	}
			//}
			//l_b = { st,degree[root] };
			p_t u_b = p;
			fill_n(is_removed.begin(), is_removed.size(), true);
			fill_n(is_already_added_into_subgraph.begin(), is_already_added_into_subgraph.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;

			for(int i = 0; i < num_node; i++) {
				if (node_kpcore_info[i].size() >= k) {
					if( compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
						fixed.insert(i);
						is_removed[i] = true;
					}
					else if( compare_pvalue(node_kpcore_info[i][k - 1], l_b) >= 0) {
						subgraph_nodes.insert(i);
						is_removed[i] = false;
					}
				}
			}

			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);

			//check  
			//printf("%d %lu\n", k, subgraph_nodes.size());

			checknum_dec += subgraph_nodes.size();
			checknum_nei += subgraph_nodes.size() + fixed.size();

			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				//if (d == 0) {
				//	p_min = { 0,1 };
				//}
				//else {
				//	p_min = make_pair(p_min.first / d, p_min.second / d);
				//}
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	return checknum_dec;
}

long long Graph::global_delete_edge_vertex_statistics(vid_t u, vid_t v) {

	long long checknum_dec = 0;

	auto start = chrono::system_clock::now();
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size();
	int k_min = MIN(k_u, k_v);
	remove_edge(u, v);
	// re-compute k_min-core
	vector<vid_t> s;
	vector<bool> is_removed;
	is_removed.resize(num_node);
	fill_n(is_removed.begin(), is_removed.size(), false);
	auto tmp_deg = degree;
	for (int i = 0; i < num_node; i++) {
		if (tmp_deg[i] < k_min) {
			s.push_back(i);
		}
	}
	while (!s.empty()) {
		vid_t cur = s.back();
		s.pop_back();
		is_removed[cur] = true;
		for (auto neigh : neighbor[cur]) {
			if (is_removed[neigh]) continue;
			tmp_deg[neigh]--;
			if (tmp_deg[neigh] == k_min - 1) {
				s.push_back(neigh);
			}
		}
	}
	bool kcore_change = false;
	for (int i = 0; i < num_node; i++) {
		if (is_removed[i] && node_kpcore_info[i].size() == k_min) {
			//p_t tt = { 0,1 };
			node_kpcore_info[i].pop_back();
			kcore_change = true;
		}
	}
	int max_k = k_min;

	// puts("POINT Z");

	if (kcore_change) {
		max_k = k_min - 1;
		//recompute the entire max_k + 1 kpcore decomposition
		vector<int>& deg_in_pcore = tmp_deg;
		vector<bool>& is_removed_in_pcore = is_removed;
		min_heap mheap;
		mheap.node2index.resize(num_node, -1);
		for (vid_t i = 0; i < num_node; i++) {
			if (is_removed_in_pcore[i]) continue;
			p_t t = make_pair(deg_in_pcore[i], degree[i]);
			mheap.insert({ t,i });
		}
		while (!mheap.empty()) {
			auto top_ele = mheap.pop();
			p_t p_min = top_ele.first;
			vector<vid_t> stack_pcore = { top_ele.second };
			while (!mheap.empty()) {
				auto te = mheap.top();
				if (compare_pvalue(p_min, te.first) < 0) {
					break;
				}
				else {
					mheap.pop();
					stack_pcore.push_back(te.second);
				}
			}
			int d = mcd(p_min.first, p_min.second);
			p_min = make_pair(p_min.first / d, p_min.second / d);
			/*if (d == 0) {
				p_min = { 0,1 };
			}
			else {
				p_min = make_pair(p_min.first / d, p_min.second / d);
			}*/
			if (stack_pcore.empty()) {
				break; // terminate
			}
			vector<vid_t> vec;
			//naive_index[k].push_back(make_pair(p_min, vec));
			while (!stack_pcore.empty()) {
				vid_t v = stack_pcore.back();
				stack_pcore.pop_back();
				if (is_removed_in_pcore[v]) continue;
				mheap.remove(mheap.node2index[v]);
				node_kpcore_info[v][max_k] = p_min;
				is_removed_in_pcore[v] = true;
				for (int i = 0; i < neighbor[v].size(); i++) {
					vid_t neigh = neighbor[v][i];
					if (is_removed_in_pcore[neigh]) {
						continue;
					}
					deg_in_pcore[neigh]--;
					if (deg_in_pcore[neigh] < max_k + 1) {
						stack_pcore.push_back(neigh);
					}
					else {
						p_t t = make_pair(deg_in_pcore[neigh], degree[neigh]);
						if (compare_pvalue(t, p_min) <= 0) {
							stack_pcore.push_back(neigh);
						}
						else {
							mheap.decrease_key(t, mheap.node2index[neigh]);
						}
					}
				}
			}
		}
	}

	// puts("POINT A");

	for (int k = 2; k <= max_k; k++) {
		//  printf("now k %d\n", k);
		if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) > 0) {
			// puts("CHOICEA");
			// rooted at v
			p_t u_b = node_kpcore_info[v][k - 1]; // revised by  
			//p_t l_b = { 0, 1};
			int st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			int st2 = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st2++;
				}
			}
			p_t l_b;

			p_t l_b_v = {0, 1};
			vector<p_t> neisort;
			for (auto neigh : neighbor[v]) {
				if(node_kpcore_info[neigh].size() >= k) 
					neisort.push_back(node_kpcore_info[neigh][k - 1]);
			}
			std::sort(neisort.begin(), neisort.end(), cmp_sort_p);
			for (int i = 0; i < neisort.size(); i++) {
				if(i + 1 >= k && compare_pvalue({i + 1, degree[v]}, neisort[i]) >= 0  && compare_pvalue(node_kpcore_info[v][k - 1], neisort[i]) >= 0) {
					l_b_v = neisort[i];
					break;
				}
			}

			// if(compare_pvalue({st, degree[v]}, {st2, degree[u]}) < 0) l_b = {st, degree[v]};
			// else l_b = {st2, degree[u]};
			if(st < k) l_b = {0, 1};
			else l_b = {st, degree[v]};
			l_b = l_b_v;
			if(k == max_k) l_b = {0, 1};
			// l_b = {0, 1};
			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			//if(k == max_k) l_b = {0,1};

			//p_t l_b = { st,degree[v] }; // revised by  
			//if(k == max_k) l_b = {0, 1};
			//p_t u_b = {1,1};
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;
			stack.push_back(v);

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			

			checknum_dec += subgraph_nodes.size();
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("hA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
			// rooted at u
			l_b = node_kpcore_info[u][k - 1];
			st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) >= 0) {
					st++;
				}
			}
			u_b = { st,degree[u] };
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			subgraph_nodes.clear();
			subgraph_degree.clear();
			fixed.clear();

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
		
			// compute degree

			checknum_dec += subgraph_nodes.size();

			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			mheap.container.resize(1);
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("Hs2");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
		else if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) < 0) {
			// puts("CHOICEB");
			// puts("ELSE LOOP");
			p_t u_b = node_kpcore_info[u][k - 1]; // revised by  
			//p_t l_b = { 0, 1};
			int st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			int st2 = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st2++;
				}
			}
			// printf("st : %d, st2 : %d\n", st, st2);
			p_t l_b;

			p_t l_b_u = {0, 1};
			vector<p_t> neisort;
			for (auto neigh : neighbor[u]) {
				if(node_kpcore_info[neigh].size() >= k) 
					neisort.push_back(node_kpcore_info[neigh][k - 1]);
			}
			std::sort(neisort.begin(), neisort.end(), cmp_sort_p);
			for (int i = 0; i < neisort.size(); i++) {
				if(i + 1 >= k && compare_pvalue({i + 1, degree[u]}, neisort[i]) >= 0  && compare_pvalue(node_kpcore_info[u][k - 1], neisort[i]) >= 0) {
					l_b_u = neisort[i];
					break;
				}
			}

			// if(compare_pvalue({st, degree[v]}, {st2, degree[u]}) < 0) l_b = {st, degree[v]};
			// else l_b = {st2, degree[u]};

			if(st2 < k) l_b = {0, 1};
			else l_b = {st2, degree[u]};
			l_b = l_b_u;
			if(k == max_k) l_b = {0, 1};
			
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
		

			checknum_dec += subgraph_nodes.size();

			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d %d %d\n", d, p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("SADAAA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
			// rooted at v
			l_b = node_kpcore_info[v][k - 1];
			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) >= 0) {
					st++;
				}
			}
			u_b = { st,degree[v] };
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			subgraph_nodes.clear();
			subgraph_degree.clear();
			fixed.clear();

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			

			checknum_dec += subgraph_nodes.size();
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			//min_heap mheap;
			mheap.container.resize(1);
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("HSA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
		else { //compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) == 0
			   // rooted at u and v the most complicate part
			   //p_t u_b = node_kpcore_info[u][k - 1];
			//    puts("CHOICEC");


			int st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st++;
				}
			}
			
			p_t l_b_u = { st,degree[u] };
			
			

			st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) >= 0) {
					st++;
				}
			}
			p_t u_b_u = { st,degree[u] };
			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			p_t l_b_v = { st,degree[v] };

			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) >= 0) {
					st++;
				}
			}
			p_t u_b_v = { st,degree[v] };
			p_t l_b = { 0,1 };
			if (compare_pvalue(l_b_u, l_b_v) < 0) {
				l_b = l_b_u;
			}
			else {
				l_b = l_b_v;
			}
			if(l_b_u.first < k || l_b_v.first < k) l_b = {0, 1};
			if(k == max_k) l_b = {0,1};

			// l_b = (l_b_u < l_b_v) ? l_b_u : l_b_v;

			p_t u_b = { 1,1 };
			if (compare_pvalue(u_b_u, u_b_v) > 0) {
				u_b = u_b_u;
			}
			else {
				u_b = u_b_v;
			}
			if (compare_pvalue(u_b, node_kpcore_info[u][k - 1]) < 0) {
				u_b = node_kpcore_info[u][k - 1];
			}
			// l_b = {0, 1};
			// u_b = {node_kpcore_info[u][k - 1]};

			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}


			checknum_dec += subgraph_nodes.size();

			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("SADA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
	}
	//process large k
	
	// puts("POINT B");

	vid_t root; int l_k = 0;
	k_v = node_kpcore_info[v].size();
	k_u = node_kpcore_info[u].size();
	if (k_v > k_u) {
		root = v;
		l_k = k_v;
	}
	else if (k_v < k_u) {
		root = u;
		l_k = k_u;
	}
	else {
		l_k = 0;
		root = 0;
	}
	int start_kval = 0;
	if (kcore_change) {
		start_kval = max_k + 2;
	}
	else {
		start_kval = max_k + 1;
	}
	for (int k = start_kval; k <= l_k; k++) {
		// compute upper bound for root
		p_t p = node_kpcore_info[root][k - 1];
		int l = 0;
		for (auto neigh : neighbor[root]) {
			if (node_kpcore_info[neigh].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) >= 0) {
				l++;
			}
		}
		p_t u_b = { l,degree[root] };
		p_t l_b = p;
		// need to recompute
		//if (compare_pvalue(p, l_b) > 0) {
		if (true) {
			
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}

			checknum_dec += subgraph_nodes.size();

			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			while (true) {
				p_t p_min = make_pair(1, 1);
				for (auto i : subgraph_nodes) {
					if (is_removed[i]) continue;
					p_t t = make_pair(subgraph_degree[i], degree[i]);
					if (compare_pvalue(p_min, t) > 0) {
						p_min = t;
					}
				}
				int d = mcd(p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				vector<vid_t> stack_pcore;
				for (auto i : subgraph_nodes) {
					if (is_removed[i]) continue;
					p_t t = make_pair(subgraph_degree[i], degree[i]);
					if (compare_pvalue(p_min, t) == 0) {
						stack_pcore.push_back(i);
					}
				}
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;
					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
						}
					}
				}
			}
		}
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	return checknum_dec;
}

double Graph::global_delete_edge(vid_t u, vid_t v) {
	auto start = chrono::system_clock::now();
	int k_u = node_kpcore_info[u].size();
	int k_v = node_kpcore_info[v].size();
	int k_min = MIN(k_u, k_v);
	remove_edge(u, v);
	// re-compute k_min-core
	vector<vid_t> s;
	vector<bool> is_removed;
	is_removed.resize(num_node);
	fill_n(is_removed.begin(), is_removed.size(), false);
	auto tmp_deg = degree;
	for (int i = 0; i < num_node; i++) {
		if (tmp_deg[i] < k_min) {
			s.push_back(i);
		}
	}
	while (!s.empty()) {
		vid_t cur = s.back();
		s.pop_back();
		is_removed[cur] = true;
		for (auto neigh : neighbor[cur]) {
			if (is_removed[neigh]) continue;
			tmp_deg[neigh]--;
			if (tmp_deg[neigh] == k_min - 1) {
				s.push_back(neigh);
			}
		}
	}
	bool kcore_change = false;
	for (int i = 0; i < num_node; i++) {
		if (is_removed[i] && node_kpcore_info[i].size() == k_min) {
			//p_t tt = { 0,1 };
			node_kpcore_info[i].pop_back();
			kcore_change = true;
		}
	}
	int max_k = k_min;

	// puts("POINT Z");

	if (kcore_change) {
		max_k = k_min - 1;
		//recompute the entire max_k + 1 kpcore decomposition
		vector<int>& deg_in_pcore = tmp_deg;
		vector<bool>& is_removed_in_pcore = is_removed;
		min_heap mheap;
		mheap.node2index.resize(num_node, -1);
		for (vid_t i = 0; i < num_node; i++) {
			if (is_removed_in_pcore[i]) continue;
			p_t t = make_pair(deg_in_pcore[i], degree[i]);
			mheap.insert({ t,i });
		}
		while (!mheap.empty()) {
			auto top_ele = mheap.pop();
			p_t p_min = top_ele.first;
			vector<vid_t> stack_pcore = { top_ele.second };
			while (!mheap.empty()) {
				auto te = mheap.top();
				if (compare_pvalue(p_min, te.first) < 0) {
					break;
				}
				else {
					mheap.pop();
					stack_pcore.push_back(te.second);
				}
			}
			int d = mcd(p_min.first, p_min.second);
			p_min = make_pair(p_min.first / d, p_min.second / d);
			/*if (d == 0) {
				p_min = { 0,1 };
			}
			else {
				p_min = make_pair(p_min.first / d, p_min.second / d);
			}*/
			if (stack_pcore.empty()) {
				break; // terminate
			}
			vector<vid_t> vec;
			//naive_index[k].push_back(make_pair(p_min, vec));
			while (!stack_pcore.empty()) {
				vid_t v = stack_pcore.back();
				stack_pcore.pop_back();
				if (is_removed_in_pcore[v]) continue;
				mheap.remove(mheap.node2index[v]);
				node_kpcore_info[v][max_k] = p_min;
				is_removed_in_pcore[v] = true;
				for (int i = 0; i < neighbor[v].size(); i++) {
					vid_t neigh = neighbor[v][i];
					if (is_removed_in_pcore[neigh]) {
						continue;
					}
					deg_in_pcore[neigh]--;
					if (deg_in_pcore[neigh] < max_k + 1) {
						stack_pcore.push_back(neigh);
					}
					else {
						p_t t = make_pair(deg_in_pcore[neigh], degree[neigh]);
						if (compare_pvalue(t, p_min) <= 0) {
							stack_pcore.push_back(neigh);
						}
						else {
							mheap.decrease_key(t, mheap.node2index[neigh]);
						}
					}
				}
			}
		}
	}

	// puts("POINT A");

	for (int k = 2; k <= max_k; k++) {
		//  printf("now k %d\n", k);
		if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) > 0) {
			// puts("CHOICEA");
			// rooted at v
			p_t u_b = node_kpcore_info[v][k - 1]; // revised by  
			//p_t l_b = { 0, 1};
			int st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			int st2 = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st2++;
				}
			}
			p_t l_b;

			p_t l_b_v = {0, 1};
			vector<p_t> neisort;
			for (auto neigh : neighbor[v]) {
				if(node_kpcore_info[neigh].size() >= k) 
					neisort.push_back(node_kpcore_info[neigh][k - 1]);
			}
			std::sort(neisort.begin(), neisort.end(), cmp_sort_p);
			for (int i = 0; i < neisort.size(); i++) {
				if(i + 1 >= k && compare_pvalue({i + 1, degree[v]}, neisort[i]) >= 0  && compare_pvalue(node_kpcore_info[v][k - 1], neisort[i]) >= 0) {
					l_b_v = neisort[i];
					break;
				}
			}

			// if(compare_pvalue({st, degree[v]}, {st2, degree[u]}) < 0) l_b = {st, degree[v]};
			// else l_b = {st2, degree[u]};
			if(st < k) l_b = {0, 1};
			else l_b = {st, degree[v]};
			l_b = l_b_v;
			if(k == max_k) l_b = {0, 1};
			// l_b = {0, 1};
			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			//if(k == max_k) l_b = {0,1};

			//p_t l_b = { st,degree[v] }; // revised by  
			//if(k == max_k) l_b = {0, 1};
			//p_t u_b = {1,1};
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;
			stack.push_back(v);

			// puts("INWHILE");

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("hA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
			// rooted at u
			l_b = node_kpcore_info[u][k - 1];
			st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) >= 0) {
					st++;
				}
			}
			u_b = { st,degree[u] };
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			subgraph_nodes.clear();
			subgraph_degree.clear();
			fixed.clear();

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			mheap.container.resize(1);
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("Hs2");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
		else if (compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) < 0) {
			// puts("CHOICEB");
			// puts("ELSE LOOP");
			p_t u_b = node_kpcore_info[u][k - 1]; // revised by  
			//p_t l_b = { 0, 1};
			int st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			int st2 = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st2++;
				}
			}
			// printf("st : %d, st2 : %d\n", st, st2);
			p_t l_b;

			p_t l_b_u = {0, 1};
			vector<p_t> neisort;
			for (auto neigh : neighbor[u]) {
				if(node_kpcore_info[neigh].size() >= k) 
					neisort.push_back(node_kpcore_info[neigh][k - 1]);
			}
			std::sort(neisort.begin(), neisort.end(), cmp_sort_p);
			for (int i = 0; i < neisort.size(); i++) {
				if(i + 1 >= k && compare_pvalue({i + 1, degree[u]}, neisort[i]) >= 0  && compare_pvalue(node_kpcore_info[u][k - 1], neisort[i]) >= 0) {
					l_b_u = neisort[i];
					break;
				}
			}

			// if(compare_pvalue({st, degree[v]}, {st2, degree[u]}) < 0) l_b = {st, degree[v]};
			// else l_b = {st2, degree[u]};

			if(st2 < k) l_b = {0, 1};
			else l_b = {st2, degree[u]};
			l_b = l_b_u;
			if(k == max_k) l_b = {0, 1};
			// l_b = {0, 1};
			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			// rooted at u
			// p_t u_b = node_kpcore_info[u][k - 1]; // revised by  , wrong !
			// //p_t l_b = {0, 1};
			
			// int st = 0;
			// for (auto neigh : neighbor[u]) {
			// 	if (node_kpcore_info[neigh].size() < k) continue;
			// 	if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
			// 		st++;
			// 	}
			// }
			// //p_t u_b = {1, 1};
			// //p_t u_b = { st,degree[u] }; // revised by  , wrong !
			// p_t l_b = { st,degree[u] }; 
			// if(k == max_k) {l_b = {0, 1}; puts("MAX_K !");}
			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			

			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d %d %d\n", d, p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("SADAAA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
			// rooted at v
			l_b = node_kpcore_info[v][k - 1];
			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) >= 0) {
					st++;
				}
			}
			u_b = { st,degree[v] };
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			subgraph_nodes.clear();
			subgraph_degree.clear();
			fixed.clear();

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			//min_heap mheap;
			mheap.container.resize(1);
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("HSA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
		else { //compare_pvalue(node_kpcore_info[u][k - 1], node_kpcore_info[v][k - 1]) == 0
			   // rooted at u and v the most complicate part
			   //p_t u_b = node_kpcore_info[u][k - 1];
			//    puts("CHOICEC");


			int st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) > 0) {
					st++;
				}
			}
			
			p_t l_b_u = { st,degree[u] };
			
			

			st = 0;
			for (auto neigh : neighbor[u]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[u][k - 1]) >= 0) {
					st++;
				}
			}
			p_t u_b_u = { st,degree[u] };
			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) > 0) {
					st++;
				}
			}
			p_t l_b_v = { st,degree[v] };


			st = 0;
			for (auto neigh : neighbor[v]) {
				if (node_kpcore_info[neigh].size() < k) continue;
				if (compare_pvalue(node_kpcore_info[neigh][k - 1], node_kpcore_info[v][k - 1]) >= 0) {
					st++;
				}
			}
			p_t u_b_v = { st,degree[v] };
			p_t l_b = { 0,1 };
			if (compare_pvalue(l_b_u, l_b_v) < 0) {
				l_b = l_b_u;
			}
			else {
				l_b = l_b_v;
			}
			if(l_b_u.first < k || l_b_v.first < k) l_b = {0, 1};
			if(k == max_k) l_b = {0,1};

			// l_b = (l_b_u < l_b_v) ? l_b_u : l_b_v;

			p_t u_b = { 1,1 };
			if (compare_pvalue(u_b_u, u_b_v) > 0) {
				u_b = u_b_u;
			}
			else {
				u_b = u_b_v;
			}
			if (compare_pvalue(u_b, node_kpcore_info[u][k - 1]) < 0) {
				u_b = node_kpcore_info[u][k - 1];
			}
			// l_b = {0, 1};
			// u_b = {node_kpcore_info[u][k - 1]};

			// printf("%lf, %lf\n", 1.0 * l_b.first / l_b.second, 1.0 * u_b.first / u_b.second);
			// re-compute kcore between l_b and u_b
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;
			vector<vid_t> stack;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			min_heap mheap;
			mheap.node2index.resize(num_node, -1);
			for (auto i : subgraph_nodes) {
				if (is_removed[i]) continue;
				p_t t = make_pair(subgraph_degree[i], degree[i]);
				mheap.insert({ t,i });
			}
			while (!mheap.empty()) {
				auto top_ele = mheap.pop();
				p_t p_min = top_ele.first;
				vector<vid_t> stack_pcore = { top_ele.second };
				while (!mheap.empty()) {
					auto te = mheap.top();
					if (compare_pvalue(p_min, te.first) < 0) {
						break;
					}
					else {
						mheap.pop();
						stack_pcore.push_back(te.second);
					}
				}
				int d = mcd(p_min.first, p_min.second);
				// printf("%d\n", d);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				// puts("SADA");
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					mheap.remove(mheap.node2index[cur]);
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;

					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
							else {
								mheap.decrease_key(t, mheap.node2index[neigh]);
							}
						}
					}
				}
			}
		}
	}
	//process large k
	
	// puts("POINT B");

	vid_t root; int l_k = 0;
	k_v = node_kpcore_info[v].size();
	k_u = node_kpcore_info[u].size();
	if (k_v > k_u) {
		root = v;
		l_k = k_v;
	}
	else if (k_v < k_u) {
		root = u;
		l_k = k_u;
	}
	else {
		l_k = 0;
		root = 0;
	}
	int start_kval = 0;
	if (kcore_change) {
		start_kval = max_k + 2;
	}
	else {
		start_kval = max_k + 1;
	}
	for (int k = start_kval; k <= l_k; k++) {
		// compute upper bound for root
		p_t p = node_kpcore_info[root][k - 1];
		int l = 0;
		for (auto neigh : neighbor[root]) {
			if (node_kpcore_info[neigh].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) >= 0) {
				l++;
			}
		}
		p_t u_b = { l,degree[root] };
		p_t l_b = p;
		// need to recompute
		//if (compare_pvalue(p, l_b) > 0) {
		if (true) {
			/*int st = 0;
			for (auto neigh : neighbor[root]) {
			if (node_kpcore_info[neigh].size() < k) continue;
			if (compare_pvalue(node_kpcore_info[neigh][k - 1], p) > 0) {
			st++;
			}
			}
			l_b = { st,degree[root] };
			p_t u_b = p;*/
			fill_n(is_removed.begin(), is_removed.size(), false);
			unordered_set<vid_t> subgraph_nodes;
			unordered_map<vid_t, int> subgraph_degree;
			unordered_set<vid_t> fixed;

			for(int i = 0; i < num_node; i++) {
				
				if(node_kpcore_info[i].size() < k || compare_pvalue(node_kpcore_info[i][k - 1], l_b) < 0) {
					is_removed[i] = true;
				}
				else if (compare_pvalue(node_kpcore_info[i][k - 1], u_b) > 0) {
					fixed.insert(i);
					is_removed[i] = true;
				}
				else {
					subgraph_nodes.insert(i);
					subgraph_degree[i] = 0;
				}
			}
			
			// compute degree
			for (auto sub_n : subgraph_nodes) {
				for (auto mn : neighbor[sub_n]) {
					if (!is_removed[mn]) {
						subgraph_degree[sub_n]++;
					}
					else if (fixed.find(mn) != fixed.end()) {
						subgraph_degree[sub_n]++;
					}
				}
			}
			// re-decomposition kpcore between l_b and u_b
			while (true) {
				p_t p_min = make_pair(1, 1);
				for (auto i : subgraph_nodes) {
					if (is_removed[i]) continue;
					p_t t = make_pair(subgraph_degree[i], degree[i]);
					if (compare_pvalue(p_min, t) > 0) {
						p_min = t;
					}
				}
				int d = mcd(p_min.first, p_min.second);
				p_min = make_pair(p_min.first / d, p_min.second / d);
				/*if (d == 0) {
					p_min = { 0,1 };
				}
				else {
					p_min = make_pair(p_min.first / d, p_min.second / d);
				}*/
				vector<vid_t> stack_pcore;
				for (auto i : subgraph_nodes) {
					if (is_removed[i]) continue;
					p_t t = make_pair(subgraph_degree[i], degree[i]);
					if (compare_pvalue(p_min, t) == 0) {
						stack_pcore.push_back(i);
					}
				}
				if (stack_pcore.empty()) {
					break; // terminate
				}
				while (!stack_pcore.empty()) {
					vid_t cur = stack_pcore.back();
					stack_pcore.pop_back();
					if (is_removed[cur]) continue;
					// modify index and node kpcore info
					auto old_p = node_kpcore_info[cur][k - 1];
					node_kpcore_info[cur][k - 1] = p_min;
					is_removed[cur] = true;
					for (int i = 0; i < neighbor[cur].size(); i++) {
						vid_t neigh = neighbor[cur][i];
						if (is_removed[neigh]) {
							continue;
						}
						subgraph_degree[neigh]--;
						if (subgraph_degree[neigh] < k) {
							stack_pcore.push_back(neigh);
						}
						else {
							p_t t = make_pair(subgraph_degree[neigh], degree[neigh]);
							if (compare_pvalue(t, p_min) <= 0) {
								stack_pcore.push_back(neigh);
							}
						}
					}
				}
			}
		}
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	return elapsed_seconds.count();
}

std::vector<vid_t> Graph::get_bfs_order() {

	bool* visit = new bool[num_node + 10];

	visit[1] = true;

	std::vector<vid_t> order;

	order.push_back(1);

	int head = 0;
	int tail = 1;

	while(head < tail) {
		
		int now = order[head++];

		for(vid_t nei : neighbor[now]) {
			if(visit[nei]) {
				continue;
			}
			else {
				order.push_back(nei);
				tail++;
				visit[nei] = true;
			}
		}
	}

	delete[] visit;

	return order;

}

void Graph::print_subgraph(std::vector<vid_t> nodes, std::string dir) {

	std::sort(nodes.begin(), nodes.end());
	std::unordered_map<vid_t, int> map;
	int n = nodes.size();
	std::vector<std::pair<vid_t, vid_t>> edges;

	for(int i = 0; i < n; i++){
		map[nodes[i]] = i;
	} 

	for(int i = 0; i < n; i++){
		vid_t now = nodes[i];
		for(vid_t nei : neighbor[now]) {
			if(map.find(nei) != map.end() and map[nei] > map[now]) {
				edges.push_back({map[now], map[nei]});
			}
		}
	}

	std::sort(edges.begin(), edges.end());

	printf("%d %lu\n", n, edges.size()); 
	
	FILE* f = freopen(dir.c_str(), "w", stdout);
	printf("%d %lu\n", n, edges.size());
	for(int i = 0; i < n; i++) printf("%d, ", i);
	printf("\n");
	for(auto x : edges) {
		printf("(%d, %d),", x.first, x.second);
	}
	printf("\n");
	fclose(f);
	return;
}

bool Graph::isEdge(vid_t u, vid_t v) {
	for (auto it = neighbor[u].begin(); it != neighbor[u].end(); it++) {
		if (*it == v) return true;
	}
	return false;
}

void Graph::generate_random_edge(string stream_location) {
	int num_edge_generate = 190;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<vid_t> dis(0, num_node - 1);
	ofstream edge_stream_file;
	edge_stream_file.open(stream_location);
	for (int i = 0; i < num_edge_generate;) {
		int u = dis(gen);
		int v = dis(gen);
		if (!isEdge(u, v) && u != v) {
			insert_edge(u, v);
			edge_stream_file << u << "	" << v << endl;
			i++;
		}
	}
}

bool Graph::compare_kpcore(Graph* g) {

	puts("Compare !");

	
	bool ok = true;
	int node = 0;
	int lastK = 0;
	for (int i = 0; i < num_node; i++) {
		if (node_kpcore_info[i].size() != g->node_kpcore_info[i].size()) {
			printf("K DIFFERs FOR NODE %d : FIRST : %ld / SECOND : %ld\n", i, node_kpcore_info[i].size(), g->node_kpcore_info[i].size());
			return false;
		}
		for (int j = node_kpcore_info[i].size() - 1; j > 0; j--) {
			if (compare_pvalue(node_kpcore_info[i][j], g->node_kpcore_info[i][j]) != 0) {

				if(lastK < j + 1){
					lastK = j + 1;
					node = i;
				}
				
				ok = false;

			}
		}
	}

	if(ok == false) {
		printf("KP VALUE DIFFERS FOR NODE %d, K = %d : FIRST : %lld/%lld, SECOND : %lld/%lld\n", node, lastK, node_kpcore_info[node][lastK - 1].first, node_kpcore_info[node][lastK - 1].second, g->node_kpcore_info[node][lastK - 1].first, g->node_kpcore_info[node][lastK - 1].second);
		return false;
	}

	puts("Compare ok !");

	return true;
}

void Graph::print_new_order(int k) {
	printf("NEW ORDER FOR K = %d -----------------\n", k);
	// for(auto x : p_order[k]) {
	// 	printf("%u(%llu/%llu)<-", x.second, x.first.first, x.first.second);
	// }
	for(auto x : new_order[k]) {
		printf("%u(%llu/%llu)<-", x, node_kpcore_info[x][k-1].first, node_kpcore_info[x][k-1].second);
	}
	printf("\n");

}

void Graph::print_rank(int k) {
	printf("RANK INFO FOR K = %d ---------------\n", k);
	for(int i = 0; i < num_node; i++) {
		if(rank[i].size() > k) {
			printf("rank[%d]=%d, ", i, rank[i][k]);
		}
	}
	printf("\n");
}

p_t Graph::get_p_with_default(vid_t u, int k) {
	if(node_kpcore_info[u].size() < k) return {0, 0};
	else return node_kpcore_info[u][k - 1];
}


vector<long long> Graph::compare_p_change_inc_dec(Graph* g, vid_t u, vid_t v) {

	puts("in");
	long long ans = 0;
	int now_node = 0;

	long long inc = 0;
	long long dec = 0;
	for(int j = 1; ; j++) {
		if(get_p_with_default(u, j + 1).second == 0 && get_p_with_default(v, j + 1).second == 0) break;
		// printf("k = %d\n", j + 1);
		// printf("u = %d, p value from %lld/%lld to %lld/%lld\n", u, g->get_p_with_default(u, j + 1).first, g->get_p_with_default(u, j + 1).second, get_p_with_default(u, j + 1).first, get_p_with_default(u, j + 1).second);
		// printf("v = %d, p value from %lld/%lld to %lld/%lld\n", v, g->get_p_with_default(v, j + 1).first, g->get_p_with_default(v, j + 1).second, get_p_with_default(v, j + 1).first, get_p_with_default(v, j + 1).second);
		now_node = 0;
		int now_change = 0;
		int up = 0; 
		for(int i = 0; i < num_node; i++) {
			int k_a = node_kpcore_info[i].size();
			int k_b = g->node_kpcore_info[i].size();

			if(k_a > j && k_b > j) {
				now_node++;
				if(compare_pvalue(node_kpcore_info[i][j], g->node_kpcore_info[i][j]) != 0) {
					// printf("k = %d, %d change from %lld/%lld to %lld/%lld %d\n", j + 1, i,
					// 	g->node_kpcore_info[i][j].first, g->node_kpcore_info[i][j].second, 
					// 	node_kpcore_info[i][j].first, node_kpcore_info[i][j].second,
					// 	compare_pvalue(node_kpcore_info[i][j], g->node_kpcore_info[i][j]) > 0 ? 1 : -1);

					if (compare_pvalue(node_kpcore_info[i][j], g->node_kpcore_info[i][j]) > 0){
						up++;
					}
					ans++;
					now_change++;
				}
			}
			else if(k_a <= j && k_b <= j) {
				// printf("k = %d, changed_nodes = %d\n", j + 1, now_change);
				continue;
			}
			else if(k_a > j){
				ans ++;
				now_node++;
				now_change++;
				up++;
				// printf("k = %d, %d change from %lld/%lld to %lld/%lld\n", j + 1, i,
				// 		0ll, 0ll, 
				// 		node_kpcore_info[i][j].first, node_kpcore_info[i][j].second);
			}
			else if(k_b > j) {
				ans++;
				now_node++;
				now_change++;
				// printf("k = %d, %d change from %lld/%lld to %lld/%lld\n", j + 1, i,
				// 		g->node_kpcore_info[i][j].first, g->node_kpcore_info[i][j].second, 
				// 		0ll, 0ll);
			}
		}
		inc += up;
		// printf("k = %d, changed_nodes = %d(up=%d, down=%d)\n", j + 1, now_change, up, now_change - up);
		if(now_node == 0) break;
	}

	// for(int i = 0; i < num_node; i++) {
	// 	int k_a = node_kpcore_info[i].size();
	// 	int k_b = g->node_kpcore_info[i].size();
	// 	int x = (k_a > k_b) ? k_b : k_a;
	// 	ans += k_a - x + k_b - x;
	// 	for(int j = 0; j < x; j++){
	// 		// printf("%d %d\n", i, j + 1);
	// 		if(compare_pvalue(node_kpcore_info[i][j], g->node_kpcore_info[i][j]) != 0) {
	// 			printf("%d %d\n", i, j + 1);
	// 			ans++;
	// 		}
	// 	}
	// }
	puts("out");
	return {inc, ans -  inc};
}


void Graph::find_smallest_k_p_set(int k, p_t p, int target) {
	vector<int> visited;
	visited.resize(num_node, 0);
	vector<int> size;
	vector<int> nump;
	vector<int> tags;
	size.push_back(0);
	nump.push_back(0);
	tags.push_back(0);
	int tag = 1;
	int mxsize = num_node;
	int mxtag = -1;
	int mx_p = 0;
	for (int i = 0; i < num_node; i++) {
		if(visited[i] == 0 && node_kpcore_info[i].size() >= k) {
			visited[i] = tag;
			int has_p = 0;
			vector<int> que;
			int head = 0;
			int tail = 0;
			que.push_back(i);
			tail++;
			while(head < tail) {
				int now = que[head++];
				for(int x : neighbor[now]) {
					if(visited[x] == 0 && node_kpcore_info[x].size() >= k) {
						que.push_back(x);
						++tail;
						visited[x] = tag;
						if(compare_pvalue(node_kpcore_info[x][k - 1], p) >= 0) {
							has_p ++;
						}
					}
				}
			}

			// if(has_p != 0 && que.size() - has_p > 3) {
				size.push_back(que.size());
				nump.push_back(has_p);
				tags.push_back(tag);
			// }

			

			// if(has_p > 0 && que.size() < mxsize) {
			// 	mxtag = tag;
			// 	mxsize = que.size();
			// 	mx_p = has_p;
			// }
			++ tag;
		}
	}

	// if(mxtag != -1) {
	// 	printf("%d %d %d\n", mxtag, mxsize, mx_p);
	// }

	for(int i = 1; i < size.size(); i++){
		printf("%d %d %d\n", tags[i], size[i], nump[i]);
	} 

	if(target == -1) return;

	unordered_set<int> node_kp;
	unordered_set<int> node_k;
	unordered_set<int> node_link;

	vector<pair<int, int>> edges;

	for (int i = 0; i < num_node; i++) {
		if (visited[i] == target) {
			if(compare_pvalue(node_kpcore_info[i][k - 1], p) >= 0) {
				node_kp.insert(i);
			} else node_k.insert(i);
			for (auto y : neighbor[i]) {
				if(node_kp.find(y) == node_kp.end() && node_k.find(y) == node_k.end()) {
					edges.push_back({i, y});
				}
			}
		}
		
	}

	for (auto e : edges) {
		if(node_kp.find(e.first) == node_kp.end() && node_k.find(e.first) == node_k.end() && node_link.find(e.first) == node_link.end()) {
			node_link.insert(e.first);
		}
		if(node_kp.find(e.second) == node_kp.end() && node_k.find(e.second) == node_k.end() && node_link.find(e.second) == node_link.end()) {
			node_link.insert(e.second);
		}
	}


	FILE * f = freopen("case_edge.csv", "w", stdout);
	printf("source,target,weight\n");
	for (auto e : edges) {
		printf("%d,%d,%d\n", e.first, e.second, 1);
	}
	fclose(f);

	f = freopen("case_node.csv", "w", stdout);
	printf("id,module\n");
	for (auto x : node_kp) {
		printf("%d,blue\n", x);
	}
	for (auto x : node_k) {
		printf("%d,green\n", x);
	}
	for (auto x : node_link) {
		printf("%d,yellow\n", x);
	}
	fclose(f);


	return;
}

void Graph::print_k_order() {
	puts("K-Order :");
	for (int i = 0; i < k_order.size(); i++) {
		printf("%d(k = %lu, deg+ = %d)->", k_order[i], node_kpcore_info[k_order[i]].size(), k_deg_plus[k_order[i]]);
	}
	puts("");
}


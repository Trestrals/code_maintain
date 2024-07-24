#include "experiment.h"

using namespace std;

void correctness_test_local_insertion(std::string dataset) {
	cout << "Test Local Insertion Correctness, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";

	g->loadGraph(dir);
	g->kpcore_decom();

	Graph* g2 = new Graph();
	g2->loadGraph(dir);
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	vid_t u, v;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->local_insert_edge(u, v);

		g2->global_insert_edge(u, v);
		cout << "total number of edges : " << g->num_edge << endl;
		if (!g->compare_kpcore(g2)) {
			cout << "error in local insertion" << endl;
			cout << "FALSE" << endl;
			return;
		}
	}
	cout << "test passed" << endl;
	return;
}

void correctness_test_local_deletion(std::string dataset) {
	cout << "Test Local Deletion Correctness, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";

	g->loadGraph(dir);
	g->kpcore_decom();

	Graph* g2 = new Graph();
	g2->loadGraph(dir);
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	vid_t u, v;
	
	while (input_edge_stream >> u >> v) {
		printf("deleting (%d, %d)\n", u, v);
		g->local_delete_edge(u, v);

		g2->global_delete_edge(u, v);
		cout << "total number of edges : " << g->num_edge << endl;
		if (!g->compare_kpcore(g2)) {
			cout << "error in local deletion" << endl;
			cout << "FALSE" << endl;
			return;
		}
	}
	cout << "test passed" << endl;
	return;
}

void correctness_test_global_insertion(std::string dataset) {
	cout << "Test Global Insertion Correctness, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";

	g->loadGraph(dir);
	g->kpcore_decom();

	Graph* g2 = new Graph();
	g2->loadGraph(dir);
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	vid_t u, v;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->global_insert_edge(u, v);

		g2->insert_edge(u, v);
		g2->kpcore_decom();
		cout << "total number of edges : " << g->num_edge << endl;
		if (!g->compare_kpcore(g2)) {
			cout << "error in global insertion" << endl;
			cout << "FALSE" << endl;
			return;
		}
	}
	cout << "test passed" << endl;
	return;
}

void correctness_test_global_deletion(std::string dataset) {
	cout << "Test Global Deletion Correctness, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";

	g->loadGraph(dir);
	g->kpcore_decom();

	Graph* g2 = new Graph();
	g2->loadGraph(dir);
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	vid_t u, v;
	
	while (input_edge_stream >> u >> v) {
		printf("deleting (%d, %d)\n", u, v);
		g->global_delete_edge(u, v);

		g2->remove_edge(u, v);
		g2->kpcore_decom();
		cout << "total number of edges : " << g->num_edge << endl;
		if (!g->compare_kpcore(g2)) {
			cout << "error in global deletion" << endl;
			cout << "FALSE" << endl;
			return;
		}
	}
	cout << "test passed" << endl;
	return;
}

void runtime_test_local_insertion(std::string dataset) {
	cout << "Test Local Insertion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->local_insert_edge(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges inserted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge insertion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void runtime_test_local_deletion(std::string dataset) {
	cout << "Test Local Deletion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("deleting (%d, %d)\n", u, v);
		g->local_delete_edge(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges deleted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge deletion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void runtime_test_global_insertion(std::string dataset) {
	cout << "Test Global Insertion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->global_insert_edge(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges inserted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge insertion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void runtime_test_global_deletion(std::string dataset) {
	cout << "Test Global Deletion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("deleting (%d, %d)\n", u, v);
		g->global_delete_edge(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges deleted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge deletion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void runtime_test_local_woc_insertion(std::string dataset) {
	cout << "Test Local-woc Insertion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->local_insert_edge_woc(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges inserted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge insertion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void runtime_test_local_woc_deletion(std::string dataset) {
	cout << "Test Local-woc Deletion Runtime, Dataset : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";
	g->loadGraph(dir);

	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	vid_t u, v;
	int cnt = 0;
	 
	auto start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;
	
	while (input_edge_stream >> u >> v) {
		printf("inserting (%d, %d)\n", u, v);
		g->local_delete_edge_woc(u, v);
		elapsed_seconds = chrono::system_clock::now() - start;
		cout << ++cnt << " edges deleted, total time : " << elapsed_seconds.count() << endl;\
	}
	
	elapsed_seconds = chrono::system_clock::now() - start;
	cout << "average time for edge deletion: " << elapsed_seconds.count() / cnt << endl;
	 
	return;
}

void vertices_statistics_local_insertion(std::string dataset) {
	cout << "Vertex Statistics for Local algorithm during edge insertion : " << dataset << endl;
	Graph* g = new Graph();
	Graph* g2 = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";
	g->loadGraph(dir);
	g2->loadGraph(dir);

	g->kpcore_decom();
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	
	vid_t u, v;
	int cnt = 0;

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	double avg_check_time_per_vertex = 0;
	vector<long long> tempvec;

	vector<long long> tmp;
	long long inc = 0;
	long long dec = 0;

	
	while (input_edge_stream >> u >> v) {

		tempvec= g->local_insert_edge_vertex_statistics(u, v);
		dec_times += tempvec[0];
		checknum_nei_dec += tempvec[1];
		checknum_dec_dec += tempvec[2];
		total_check_number_dec += tempvec[3];
		distinct_check_num_dec += tempvec[4];
		inc_times += tempvec[5];
		checknum_nei_inc += tempvec[6];
		checknum_dec_inc += tempvec[7];
		avg_check_time_per_vertex += tempvec[8] * 1.0 / 1000;

		tmp = g->compare_p_change_inc_dec(g2, u, v);
		inc += tmp[0];
		dec += tmp[1];

		g2->local_insert_edge(u, v);

		++cnt;
		cout << cnt << " edges inserted." << endl;
		cout << "total distinct vertices maintained for pDecrease this time : " << tempvec[4] << endl;
		cout << "total vertices maintained for pDecrease this time : " << tempvec[3] << endl;
		cout << "total vertices maintained for pIncrease this time : " << tempvec[7] << endl;
		cout << "total vertices with pn changed this time : " << tmp[0] + tmp[1] << endl;
		
	}
	cout << "------------------------------------------------------" << endl;
	cout << "Process end, " << cnt << " edges inserted." << endl;
	cout << "total distinct vertices maintained for pDecrease: " << distinct_check_num_dec << endl;
	cout << "total vertices maintained for pDecrease: " << total_check_number_dec << endl;
	cout << "total vertices maintained for pIncrease: " << checknum_dec_inc << endl;
	cout << "total vertices maintained: " << checknum_dec_inc + total_check_number_dec << endl;
	cout << "total vertices with pn changed: " << inc + dec << endl;
	cout << "average check time for each distinct vertices during all pDecrease process: " << 1.0 * total_check_number_dec / distinct_check_num_dec << endl;
	
	return;
}

void vertices_statistics_local_deletion(std::string dataset) {
	cout << "Vertex Statistics for Local algorithm during edge deletion : " << dataset << endl;
	Graph* g = new Graph();
	Graph* g2 = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";
	g->loadGraph(dir);
	g2->loadGraph(dir);

	g->kpcore_decom();
	g2->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	
	vid_t u, v;
	int cnt = 0;

	long long checknum_nei_dec = 0;
	long long checknum_dec_dec = 0;
	long long checknum_nei_inc = 0;
	long long checknum_dec_inc = 0;
	long long inc_times = 0;
	long long dec_times = 0;
	long long total_check_number_dec = 0;
	long long distinct_check_num_dec = 0;
	double avg_check_time_per_vertex = 0;
	vector<long long> tempvec;

	vector<long long> tmp;
	long long inc = 0;
	long long dec = 0;

	
	while (input_edge_stream >> u >> v) {

		tempvec= g->local_delete_edge_vertex_statistics(u, v);
		dec_times += tempvec[0];
		checknum_nei_dec += tempvec[1];
		checknum_dec_dec += tempvec[2];
		total_check_number_dec += tempvec[3];
		distinct_check_num_dec += tempvec[4];
		inc_times += tempvec[5];
		checknum_nei_inc += tempvec[6];
		checknum_dec_inc += tempvec[7];
		avg_check_time_per_vertex += tempvec[8] * 1.0 / 1000;

		tmp = g2->compare_p_change_inc_dec(g, u, v);
		inc += tmp[1];
		dec += tmp[0];

		g2->local_delete_edge(u, v);

		++cnt;
		cout << cnt << " edges deleted." << endl;
		cout << "total distinct vertices maintained for pDecrease this time : " << tempvec[4] << endl;
		cout << "total vertices maintained for pDecrease this time : " << tempvec[3] << endl;
		cout << "total vertices maintained for pIncrease this time : " << tempvec[7] << endl;
		cout << "total vertices with pn changed this time : " << tmp[0] + tmp[1] << endl;
		
	}
	cout << "------------------------------------------------------" << endl;
	cout << "Process end, " << cnt << " edges deleted." << endl;
	cout << "total distinct vertices maintained for pDecrease: " << distinct_check_num_dec << endl;
	cout << "total vertices maintained for pDecrease: " << total_check_number_dec << endl;
	cout << "total vertices maintained for pIncrease: " << checknum_dec_inc << endl;
	cout << "total vertices maintained: " << checknum_dec_inc + total_check_number_dec << endl;
	cout << "total vertices with pn changed: " << inc + dec << endl;
	cout << "average check time for each distinct vertices during all pDecrease process: " << 1.0 * total_check_number_dec / distinct_check_num_dec << endl;
	
	return;
}

void vertices_statistics_global_insertion(std::string dataset) {
	cout << "Vertex Statistics for Global algorithm during edge insertion : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "2.txt";
	string dir2 = "../datasets/_" + dataset + "2_addGraph.txt";
	g->loadGraph(dir);
	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	
	vid_t u, v;
	int cnt = 0;

	long long total_check_number = 0;
	long long temp;

	while (input_edge_stream >> u >> v) {

		temp = g->global_insert_edge_vertex_statistics(u, v);
		total_check_number += temp;

		++cnt;
		cout << cnt << " edges inserted." << endl;
		cout << "total vertices maintained this time : " << temp << endl;
		
	}
	cout << "------------------------------------------------------" << endl;
	cout << "Process end, " << cnt << " edges inserted." << endl;
	cout << "total vertices maintained: " << total_check_number << endl;
	
	return;
}

void vertices_statistics_global_deletion(std::string dataset) {
	cout << "Vertex Statistics for Global algorithm during edge deletion : " << dataset << endl;
	Graph* g = new Graph();
	string dir = "../datasets/_" + dataset + "full.txt";
	string dir2 = "../datasets/_" + dataset + "full_delGraph.txt";
	g->loadGraph(dir);
	g->kpcore_decom();

	ifstream input_edge_stream;
	input_edge_stream.open(dir2);
	
	
	vid_t u, v;
	int cnt = 0;

	long long total_check_number = 0;
	long long temp;

	while (input_edge_stream >> u >> v) {

		temp = g->global_delete_edge_vertex_statistics(u, v);
		total_check_number += temp;

		++cnt;
		cout << cnt << " edges deleted." << endl;
		cout << "total vertices maintained this time : " << temp << endl;
		
	}
	cout << "------------------------------------------------------" << endl;
	cout << "Process end, " << cnt << " edges deleted." << endl;
	cout << "total vertices maintained: " << total_check_number << endl;
	
	return;
}

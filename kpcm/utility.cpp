#include "utility.h"

using namespace std;

double insert_total_time;
double decrease_total_time;

//std::vector<std::string> Datasets = { "Facebook", "Brightkite", "Gowalla", "YouTube",
//"Pokec", "DBLP", "LiveJournal", "Orkut" };

int mcd(int m, int n)
{
	int r = m;
	while (r != 0)
	{
		r = m % n;
		m = n;
		n = r;
	}
	return m;
}

// return -1 if p1 < p2; 0 if p1 == p2; 1 if p1 > p2
int compare_pvalue(const p_t &p1, const p_t &p2) {
	//if (p1.first == 0 || p2.first == 0 || p1.second == 0 || p2.second == 0) return 999;
	if (p1.first * p2.second == p1.second * p2.first) {
		return 0;
	}
	else if (p1.first * p2.second > p1.second * p2.first) {
		return 1;
	}
	else {
		return -1;
	}
}

bool comp_p_desc(const p_t &p1, const p_t &p2) {
	if (compare_pvalue(p1, p2) > 0) {
		return true;
	}
	return false;
}


p_t p_minus(const p_t &p_1, const p_t &p_2) {
	p_t res = make_pair(p_1.first * p_2.second - p_2.first * p_1.second, p_1.second * p_2.second);
	return res;
}

void min_heap::insert(std::pair<p_t, vid_t> e) {
	// auto bg = chrono::system_clock::now();
	container.push_back(e);
	int pos = container.size() - 1;
	node2index[e.second] = pos;
	while (pos != 1) {
		int parent = pos / 2;
		if (compare_pvalue(container[parent].first, e.first) > 0) {
			swap(node2index[e.second], node2index[container[parent].second]);
			swap(container[parent], container[pos]);
			pos /= 2;
		}
		else {
			break;
		}
	}
	// chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - bg;
	// insert_total_time += elapsed_seconds.count();
}

void min_heap::decrease_key(p_t newkey, int index) {
	// auto bg = chrono::system_clock::now();
	container[index].first = newkey;
	int pos = index;
	while (pos != 1) {
		int parent = pos / 2;
		if (compare_pvalue(container[parent].first, container[pos].first) > 0) {
			swap(node2index[container[pos].second], node2index[container[parent].second]);
			swap(container[parent], container[pos]);
			pos /= 2;
		}
		else {
			break;
		}
	}
	// chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - bg;
	// decrease_total_time += elapsed_seconds.count();
}

void min_heap::increase_key(p_t newkey, int index) {
	int max_po = container.size();
	container[index].first = newkey;
	int pos = index;

	int left = pos << 1;

	while(left < max_po) {
		if((left | 1) < max_po && compare_pvalue(container[left | 1].first, container[left].first) < 0) {
			left |= 1;
		}
		if(compare_pvalue(container[left].first, container[pos].first) < 0) {
			std::swap(node2index[container[pos].second], node2index[container[left].second]);
	 		std::swap(container[left], container[pos]);
			pos = left;
			left = pos << 1;
		}
		else break;
	}

}

pair<p_t, vid_t> min_heap::pop() {
	auto top_element = container[1];
	swap(node2index[container[1].second], node2index[container[container.size() - 1].second]);
	swap(container[1], container[container.size() - 1]);
	container.pop_back();
	node2index[top_element.second] = -1;
	int pos = 1;
	while (pos * 2 < container.size()){
		// choose the smallest child
		int child = -1;
		if (container.size() <= pos * 2 + 1) {
			child = pos * 2;
		}
		else {
			if (compare_pvalue(container[pos * 2].first, container[pos * 2 + 1].first) > 0) {
				child = pos * 2 + 1;
			}
			else {
				child = pos * 2;
			}
		}
		if (compare_pvalue(container[pos].first, container[child].first) > 0) {
			swap(node2index[container[pos].second], node2index[container[child].second]);
			swap(container[pos], container[child]);
			pos = child;
		}
		else {
			break;
		}
	}
	return top_element;
}

pair<p_t, vid_t> min_heap::top() {
	if (container.size() == 1) {
		return container[0];
	}
	else {
		return container[1];
	}
}

bool min_heap::empty() {
	if (container.size() == 1) {
		return true;
	}
	else {
		return false;
	}
}
void min_heap::remove(int index) {
	if (index == -1) return; // already being removed
	decrease_key({ -1,1 }, index);
	pop();
}


min_heap::min_heap() {
	vid_t k = MAX_DEG;
	container.push_back({ { 0,0 },k });
}

void max_heap::insert(std::pair<p_t, vid_t> e) {
	// auto bg = chrono::system_clock::now();
	container.push_back(e);
	int pos = container.size() - 1;
	node2index[e.second] = pos;
	while (pos != 1) {
		int parent = pos / 2;
		if (compare_pvalue(container[parent].first, e.first) < 0) {
			swap(node2index[e.second], node2index[container[parent].second]);
			swap(container[parent], container[pos]);
			pos /= 2;
		}
		else {
			break;
		}
	}
	// chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - bg;
	// insert_total_time += elapsed_seconds.count();
}

void max_heap::increase_key(p_t newkey, int index) {
	// auto bg = chrono::system_clock::now();
	container[index].first = newkey;
	int pos = index;
	while (pos != 1) {
		int parent = pos / 2;
		if (compare_pvalue(container[parent].first, container[pos].first) < 0) {
			swap(node2index[container[pos].second], node2index[container[parent].second]);
			swap(container[parent], container[pos]);
			pos /= 2;
		}
		else {
			break;
		}
	}
	// chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - bg;
	// decrease_total_time += elapsed_seconds.count();
}

void max_heap::decrease_key(p_t newkey, int index) {
	int max_po = container.size();
	container[index].first = newkey;
	int pos = index;

	int left = pos << 1;

	while(left < max_po) {
		if((left | 1) < max_po && compare_pvalue(container[left | 1].first, container[left].first) > 0) {
			left |= 1;
		}
		if(compare_pvalue(container[left].first, container[pos].first) > 0) {
			std::swap(node2index[container[pos].second], node2index[container[left].second]);
	 		std::swap(container[left], container[pos]);
			pos = left;
			left = pos << 1;
		}
		else break;
	}

}

pair<p_t, vid_t> max_heap::pop() {
	auto top_element = container[1];
	swap(node2index[container[1].second], node2index[container[container.size() - 1].second]);
	swap(container[1], container[container.size() - 1]);
	container.pop_back();
	node2index[top_element.second] = -1;
	int pos = 1;
	while (pos * 2 < container.size()){
		// choose the smallest child
		int child = -1;
		if (container.size() <= pos * 2 + 1) {
			child = pos * 2;
		}
		else {
			if (compare_pvalue(container[pos * 2].first, container[pos * 2 + 1].first) < 0) {
				child = pos * 2 + 1;
			}
			else {
				child = pos * 2;
			}
		}
		if (compare_pvalue(container[pos].first, container[child].first) < 0) {
			swap(node2index[container[pos].second], node2index[container[child].second]);
			swap(container[pos], container[child]);
			pos = child;
		}
		else {
			break;
		}
	}
	return top_element;
}

pair<p_t, vid_t> max_heap::top() {
	if (container.size() == 1) {
		return container[0];
	}
	else {
		return container[1];
	}
}

bool max_heap::empty() {
	if (container.size() == 1) {
		return true;
	}
	else {
		return false;
	}
}
void max_heap::remove(int index) {
	if (index == -1) return; // already being removed
	increase_key({ 2,1 }, index); 
	pop();
}


max_heap::max_heap() {
	vid_t k = MAX_DEG;
	container.push_back({ { 0,0 },k });
}
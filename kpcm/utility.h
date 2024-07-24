#pragma once
#ifndef UTILITY_H
#define UTILITY_H
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <random>
#include <chrono>
#include <queue>
#include <sstream>
#include <deque>

typedef unsigned int vid_t;

typedef std::pair<long long int, long long int> p_t;

#define DOUBLE(a) (double(a.first) / a.second)

#define MIN(a, b) a < b ? a : b

#define MAX(a, b) a > b ? a : b

#define MAX_DEG 1e9

#define BINARY true

extern double insert_total_time;
extern double decrease_total_time;

extern int mcd(int m, int n);

// return -1 if p1 < p2; 0 if p1 == p2; 1 if p1 > p2
extern int compare_pvalue(const p_t &p1, const p_t &p2);

extern bool comp_p_desc(const p_t &p1, const p_t &p2);

extern p_t p_minus(const p_t &p_1, const p_t &p_2);


// index starting from 1 not 0
class min_heap{
public:
	std::vector<int> node2index;
	std::vector<std::pair<p_t, vid_t>> container;
	void insert(std::pair<p_t, vid_t> e); // return the index of insert element
	void decrease_key(p_t newkey, int index);
	void increase_key(p_t newkey, int index);
	void remove(int index);
	std::pair<p_t, vid_t> pop();
	std::pair<p_t, vid_t> top();
	bool empty();
	min_heap();
};

class max_heap{
public:
	std::vector<int> node2index;
	std::vector<std::pair<p_t, vid_t>> container;
	void insert(std::pair<p_t, vid_t> e); // return the index of insert element
	void decrease_key(p_t newkey, int index);
	void increase_key(p_t newkey, int index);
	void remove(int index);
	std::pair<p_t, vid_t> pop();
	std::pair<p_t, vid_t> top();
	bool empty();
	max_heap();
};

#endif // !UTILITY_H

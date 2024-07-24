#pragma once
#ifndef EXPERIMENT_H
#define EXPERIMENT_H
#include "graph.h"

void correctness_test_local_insertion(std::string dataset);

void correctness_test_local_deletion(std::string dataset);

void runtime_test_local_insertion(std::string dataset);

void runtime_test_local_deletion(std::string dataset);

void runtime_test_local_woc_insertion(std::string dataset);

void runtime_test_local_woc_deletion(std::string dataset);

void vertices_statistics_local_insertion(std::string dataset);

void vertices_statistics_local_deletion(std::string dataset);



#endif // !EXPERIMENT_H

#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <climits>
#include <iostream>
#include <memory>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <fstream>
#include <math.h>
#include <cmath>
#include <map>
#include "boost/multi_array.hpp"
#include <cassert>
#include <time.h>
float** perform_DP(std::vector<int>, std::vector<int>, std::vector<std::vector<int>>, float, std::unordered_set<int> );
//std::tuple<float**, std::vector<int>> perform_DP(std::vector<int>, std::vector<int>, std::vector<std::vector<int>>, float, std::unordered_set<int> );
//std::tuple<float**,int> DP_partial(std::vector<int> , std::vector<int> , std::vector<std::vector<int>> , float, int, int, float**, std::unordered_set<int> , std::map<int,int> );
//std::tuple<std::vector<int>, std::vector<int>> backtrace_partial(float** , int, int, float**, std::unordered_set<int>, std::map<int,int>, int);
float mut(int, int, int, int, float);
void compute_scoring(char*, char*, char*, float, char*, char*);
std::tuple<float, int> minimum(float*, int, int, bool, int);

#endif
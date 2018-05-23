#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <climits>
#include <iostream>
#include <memory>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <fstream>
#include <math.h>
#include "boost/multi_array.hpp"
#include <cassert>
#include <time.h>
float** printlist(std::vector<int>, std::vector<int>, std::vector<std::vector<int>>, float, std::set<int> );
float mut(int, int, int, int, float);
void compute_scoring(char*, char*, char*, float, char*, char*);
std::tuple<float, int> minimum(float*, int, int, bool, int);
//void compute_score(std::vector<int>, std::vector<int>, std::set<int>, std::vector<std::vector<int> > , std::string, float, float);

#endif
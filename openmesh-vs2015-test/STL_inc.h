#pragma once
#ifndef __STL_INC_H
#define __STL_INC_H

// #undef min
// #undef max

#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <fstream>
#include <sstream>
#include <float.h>
#include <math.h>

#ifndef MAX
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef  MIN
#define MIN(a,b)         (((a) < (b)) ? (a) : (b))
#endif


#ifndef MATH_VALUES
#define MATH_VALUES
#define M_INV_255 0.00392156862745098039f
#define M_INV_PI 0.31830988618379067154f
#define M_INF_BIG FLT_MAX
#define M_INF_SMALL FLT_MIN
#define M_EPSILON FLT_EPSILON
#endif

#define GP_PI 3.141592653  

using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

struct IntersPoint
{
	double  t_Para;
	int     trangle_num;
	/*double  weight_;*/
};



#endif
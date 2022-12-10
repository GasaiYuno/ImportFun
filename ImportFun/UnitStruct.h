#pragma once
#ifndef _UNITSTRUCT_H_
#define _UNITSTRUCT_H_

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#pragma warning (disable: 4996)
#pragma warning (disable: 4018)

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif // !PI

#ifndef EPS
#define EPS 0.001
#endif // !EPS

//新框架的数据结构
struct Unit {
    double x, y, b;
    int iFlag;
};

#endif // !_UNITSTRUCT_H_
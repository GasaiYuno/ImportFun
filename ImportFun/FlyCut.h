#pragma once
#ifndef _FLYCUT_H_
#define _FLYCUT_H_

#include "UnitStruct.h"
#include "cavc\vector2.hpp"
#include "cavc\polyline.hpp"
#include "cavc/polylineintersects.hpp"
#include "cavc\mathutils.hpp"

using namespace std;
using namespace cavc;

//数据结构
struct Line
{
	cavc::Vector2<double> P1;
	cavc::Vector2<double> P2;
	//true为正走，false为反走
	bool mode = false;
	double ka, sangle, eangle, length, cornerlength;
};

struct FlyCutInfo
{
	//起点位置
	int startcutlocality;
	//容差，最大飞切间距
	double error, flycutspace;
};

struct  FlyCutInformation
{
	//轮廓个数
	int size1;
	//每个轮廓的图元数量
	int* size2;
	//每个轮廓的封闭属性
	bool* closelist;
	//每个轮廓的图元信息
	Unit* units;
	//返回的图元个数
	int* returnsize;
	//返回的图元信息
	Unit* returnunits;
	//剩下的图元个数
	int* remainsize;
	//剩下的图元信息
	Unit* remainunits;
};

//排序
class minxcomp {
public:
	bool operator()(Line l1, Line l2) {
		return l1.P1.x() < l2.P1.x();
	}
};
class maxxcomp {
public:
	bool operator()(Line l1, Line l2) {
		return l1.P1.x() > l2.P1.x();
	}
};
class minycomp {
public:
	bool operator()(Line l1, Line l2) {
		return l1.P1.y() < l2.P1.y();
	}
};
class maxycomp {
public:
	bool operator()(Line l1, Line l2) {
		return l1.P1.y() > l2.P1.y();
	}
};
#endif // !_FLYCUT_H_

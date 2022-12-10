#pragma once
#ifndef _BRIDGE_H_
#define _BRIDGE_H_

#include "UnitStruct.h"
#include "cavc/polylinecombine.hpp"
#include "cavc/polyline.hpp"
#include "cavc/polylineintersects.hpp"

using namespace std;
using namespace cavc;

struct BridgeInfoNew
{
	//相交的轮廓个数
	int size;
	//每个相交轮廓图元数量
	int* plinesize;
	//每个轮廓信息
	Unit* pline;
	//返回的轮廓个数
	int* returnnums;
	//返回的每个轮廓图元数量
	int* returnplinesize;
	//返回的轮廓信息
	Unit* returnpline;
	//起始点
	double sx, sy;
	//结束点
	double ex, ey;
	//宽度
	double wide;
	//桥接方式
	int state;
};

#endif //!_BRIDGE_H_



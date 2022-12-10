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

//���ݽṹ
struct Line
{
	cavc::Vector2<double> P1;
	cavc::Vector2<double> P2;
	//trueΪ���ߣ�falseΪ����
	bool mode = false;
	double ka, sangle, eangle, length, cornerlength;
};

struct FlyCutInfo
{
	//���λ��
	int startcutlocality;
	//�ݲ�����м��
	double error, flycutspace;
};

struct  FlyCutInformation
{
	//��������
	int size1;
	//ÿ��������ͼԪ����
	int* size2;
	//ÿ�������ķ������
	bool* closelist;
	//ÿ��������ͼԪ��Ϣ
	Unit* units;
	//���ص�ͼԪ����
	int* returnsize;
	//���ص�ͼԪ��Ϣ
	Unit* returnunits;
	//ʣ�µ�ͼԪ����
	int* remainsize;
	//ʣ�µ�ͼԪ��Ϣ
	Unit* remainunits;
};

//����
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

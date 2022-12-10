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
	//�ཻ����������
	int size;
	//ÿ���ཻ����ͼԪ����
	int* plinesize;
	//ÿ��������Ϣ
	Unit* pline;
	//���ص���������
	int* returnnums;
	//���ص�ÿ������ͼԪ����
	int* returnplinesize;
	//���ص�������Ϣ
	Unit* returnpline;
	//��ʼ��
	double sx, sy;
	//������
	double ex, ey;
	//���
	double wide;
	//�Žӷ�ʽ
	int state;
};

#endif //!_BRIDGE_H_



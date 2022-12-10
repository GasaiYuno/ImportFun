#pragma once
#ifndef _LEADLINE_H_
#define _LEADLINE_H_

#include "UnitStruct.h"
#include "cavc/polyline.hpp"
#include "cavc/polylineintersects.hpp"

using namespace std;
using namespace cavc;

//�¿���������ݽṹ
struct LeadInOutInfo
{
    int iTailInMode;     //0.��;1.ֱ��;2.Բ��;3.ֱ�߼�Բ��
    int iTailOutMode;    //0.��;1.ֱ��;2.Բ��;3.ֱ�߼�Բ��
    double dAngleIn;     //����Ƕ�
    double dAngleOut;    //�����Ƕ�
    double dLengthIn;    //����ֱ�߳���
    double dLengthOut;   //����ֱ�߳���
    double dRadiusIn;   //����Բ���뾶
    double dRadiusOut;   //����Բ���뾶 
    int iPostionMode;    //����λ��
    int iOption;         //ѡ��
};

/// leadline_error = 0; ��ʾ����
/// leadline_error = 1�� ��ʾ������ Բ�װ뾶�������߳���
/// leadline_error = 2�� ��ʾ�����߳���С��EPS��0.01��������
/// leadline_error = 3;  ��ʾ�����߽Ƕ�С��EPS��0.01��������
/// leadline_error = 4�� ��ʾ�����߳���С��EPS��0.01��������
/// leadline_error = 5�� ��ʾ�����߽Ƕ�С��EPS��0.01��������
/// leadline_error = 6;  ��ʾ�������в嵽���
/// leadline_error = 7;  ��ʾ�������в嵽���



#endif //!_LEADLINE_H_

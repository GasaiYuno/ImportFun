#pragma once
#ifndef _LEADLINE_H_
#define _LEADLINE_H_

#include "UnitStruct.h"
#include "cavc/polyline.hpp"
#include "cavc/polylineintersects.hpp"

using namespace std;
using namespace cavc;

//新框架引线数据结构
struct LeadInOutInfo
{
    int iTailInMode;     //0.无;1.直线;2.圆弧;3.直线加圆弧
    int iTailOutMode;    //0.无;1.直线;2.圆弧;3.直线加圆弧
    double dAngleIn;     //引入角度
    double dAngleOut;    //引出角度
    double dLengthIn;    //引入直线长度
    double dLengthOut;   //引出直线长度
    double dRadiusIn;   //引出圆弧半径
    double dRadiusOut;   //引出圆弧半径 
    int iPostionMode;    //引线位置
    int iOption;         //选项
};

/// leadline_error = 0; 表示正常
/// leadline_error = 1； 表示引入线 圆孔半径大于引线长度
/// leadline_error = 2； 表示引入线长度小于EPS（0.01），错误
/// leadline_error = 3;  表示引入线角度小于EPS（0.01），错误
/// leadline_error = 4； 表示引出线长度小于EPS（0.01），错误
/// leadline_error = 5； 表示引出线角度小于EPS（0.01），错误
/// leadline_error = 6;  表示引入线有插到零件
/// leadline_error = 7;  表示引出线有插到零件



#endif //!_LEADLINE_H_

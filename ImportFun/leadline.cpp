#include "leadline.h"

//计算角度
static Vector2<double> judgeAngle(double ang)
{
    double sinx = sin(ang);
    double cosx = cos(ang);
    double tanx = tan(ang);
    if (utils::fuzzyEqual(sinx, 0.0)) sinx = 0;
    if (utils::fuzzyEqual(cosx, 0.0)) cosx = 0;
    if ((sinx > 0 && cosx > 0) || (sinx < 0 && cosx > 0))
        return { 1, tanx };
    else if ((sinx > 0 && cosx < 0) || (sinx < 0 && cosx < 0))
        return { -1, -tanx };
    else if (utils::fuzzyEqual(sinx, 0.0) && cosx > 0)
        return { 1,0 };
    else if (utils::fuzzyEqual(sinx, 0.0) && cosx < 0)
        return { -1,0 };
    else if (utils::fuzzyEqual(cosx, 0.0) && sinx > 0)
        return { 0,1 };
    else
        return { 0,-1 };
}

//判断引入引出线是否与环相交
bool IsLineInLoop(const cavc::Polyline<double>& pline, const cavc::Vector2<double>& point) {

    bool flag = false;
    // windingNumber算法
    int wn = cavc::getWindingNumber(pline, point);

    if (wn != 0) {
        flag = true;
    }
    return flag;
}

//判断是否为正常引线
bool IsNormalLeadLine(cavc::Polyline<double> pline, cavc::Polyline<double> leadline, bool isinloop, bool inorout)
{
    cavc::Vector2<double> CutUnit, EndUnit;
    cavc::PlineIntersectsResult<double> intersects;
    cavc::findIntersects(pline, leadline, cavc::createApproxSpatialIndex(pline), intersects);
    int nums = 0;
    if (inorout) {
        CutUnit = leadline.lastVertex().pos();
        EndUnit = leadline[0].pos();
    }
    else
    {
        CutUnit = leadline[0].pos();
        EndUnit = leadline.lastVertex().pos();
    }
    for (auto& intersect : intersects.intersects) {
        if (!cavc::fuzzyEqual(intersect.pos, CutUnit, EPS))
            nums++;
    }
    for (auto& confident : intersects.coincidentIntersects) {
        if (!cavc::fuzzyEqual(confident.point1, CutUnit, EPS))
            nums++;
        if (!cavc::fuzzyEqual(confident.point2, CutUnit, EPS))
            nums++;
    }
    //如果引线与轮廓没有交点
    if (nums == 0) {
        //判断当内轮廓时引线是否在轮廓外，或者反之
        if ((IsLineInLoop(pline, EndUnit) && !isinloop && pline.isClosed()) ||
            (!IsLineInLoop(pline, EndUnit) && isinloop && pline.isClosed())) {
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

//引线自动伸缩
bool LeadLineShrink(cavc::Polyline<double> pline, cavc::Polyline<double>& leadline, bool isinloop, bool inorout)
{
    cavc::PlineIntersectsResult<double> result;
    cavc::Vector2<double> basepoint, translatepoint;

    double length = cavc::getPathLength(leadline);
    if (length <= 0.1)
        return false;

    if (inorout) {
        basepoint.x() = leadline.lastVertex().pos().x();
        basepoint.y() = leadline.lastVertex().pos().y();
    }
    else {
        basepoint.x() = leadline[0].pos().x();
        basepoint.y() = leadline[0].pos().y();
    }
    cavc::StaticSpatialIndex<double> Index = cavc::createApproxSpatialIndex(pline);
    cavc::findIntersects(pline, leadline, Index, result);
    int nums = result.intersects.size() + result.coincidentIntersects.size();
    if (nums == 0)
        cavc::invertDirection(pline);
    cavc::findIntersects(pline, leadline, Index, result);
    nums = result.intersects.size() + result.coincidentIntersects.size();
    for (auto& point : result.intersects) {
        if (cavc::fuzzyEqual(point.pos, basepoint, EPS))
            nums--;
    }
    for (auto& point : result.coincidentIntersects) {
        if ((cavc::fuzzyEqual(point.point1, basepoint, EPS)) &&
            (cavc::fuzzyEqual(point.point2, basepoint, EPS)))
            nums--;
    }
    bool isnormal = IsNormalLeadLine(pline, leadline, isinloop, inorout);
    if (nums > 0 || isnormal) {
        cavc::scalePolyline(leadline, 0.5);
        if (inorout) {
            translatepoint.x() = leadline.lastVertex().pos().x();
            translatepoint.y() = leadline.lastVertex().pos().y();
        }
        else {
            translatepoint.x() = leadline[0].pos().x();
            translatepoint.y() = leadline[0].pos().y();
        }
        cavc::translatePolyline(leadline, basepoint - translatepoint);
        if (LeadLineShrink(pline, leadline, isinloop, inorout))
            return true;
        else
            return false;
    }
    else
        return true;
}

//生成相对应的引线
void CreatInLineNew(cavc::Polyline<double> _pline, cavc::Vector2<double> startcut, LeadInOutInfo info, bool inloop,cavc::Polyline<double> &InLine,int & leadline_error) {

    cavc::ClosestPoint<double> result(_pline, startcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point();
    Vector2<double> prev, midv, nextv, startv;
    PlineVertex<double> currentUnit, preUnit, nextUnit;
    bool InLoop = inloop, ischangeline = false, iscww = false;
    int cutNo = result.index();
    double len, startAngel, ang, angle;

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    len = info.dLengthIn;
    if (info.dAngleIn <= 0)
        info.dAngleIn = 1.0;
    else if (info.dAngleIn >= 180)
        info.dAngleIn = 179.0;
    angle = info.dAngleIn / 180.0 * PI;

    if (std::fabs(angle - 0) < EPS)
    {
        leadline_error = 3;
        return ;
    }
    if (fabs(len - 0) < EPS)
    {
        leadline_error = 2;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    preUnit = _pline[preno];
    currentUnit = _pline[cutNo];
    nextUnit = _pline[nextno];

    //判断起点是否是顶点且轮廓为封闭
    if ((cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS)||
        cavc::fuzzyEqual(CutUnit,_pline[nextno].pos(),EPS)) && LoopLine.isClosed()) {
        //考虑cavc求解最近点有误差
        if (cavc::fuzzyEqual(CutUnit, _pline[nextno].pos(),EPS)) {
            cutNo = nextno;
            preno = cutNo - 1;
            nextno = cutNo + 1;
            if (preno < 0)
                preno = _pline.size() - 1;
            if (nextno == _pline.size())
                nextno = 0;
            preUnit = _pline[preno];
            currentUnit = _pline[cutNo];
            nextUnit = _pline[nextno];
        }
        prev = -segTangentVector(preUnit, currentUnit, currentUnit.pos());
        if (cavc::length(prev) == 0) { leadline_error = 2;  return; }
        cavc::normalize(prev);
        nextv = segTangentVector(currentUnit, nextUnit, currentUnit.pos());
        if (cavc::length(nextv) == 0) { leadline_error = 2;  return; }
        cavc::normalize(nextv);
        //判断两向量的方向
        ///180度后者0度
        if (cavc::fuzzyEqual(-prev, nextv,EPS)|| cavc::fuzzyEqual(prev, nextv,EPS)) {
            ischangeline = true;
        }
        ///如果两条向量带角度，则取两条向量的中点
        else {
            midv = cavc::midpoint(prev, nextv);
            ///如果是内轮廓引线在外部或者外轮廓引线在内部则顺着下一个点的方向的反向引入
            if ((InLoop && !IsLineInLoop(_pline, CutUnit + midv)) ||
                (!InLoop && IsLineInLoop(_pline, CutUnit + midv)))
                midv = -nextv;
        }
        startv = CutUnit + midv * len;
    }
    //如果起点不在顶点上，或者需要该引入方向，或者为不封闭图形
    if(ischangeline||(!cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS) && !cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) || !LoopLine.isClosed()) {
        nextv = segTangentVector(currentUnit, nextUnit, CutUnit);
        cavc::normalize(nextv);
        prev = -nextv;
        // 这里的angle是计算两个点形成的直线与向量（1，0）之间的角度
        startAngel = cavc::angle(CutUnit, CutUnit + prev);
        if (InLoop)
            if(iscww)
                ang = startAngel + angle;
            else
                ang = startAngel - angle;
        else
            if(iscww)
                ang = startAngel - angle;
            else
                ang = startAngel + angle;

        midv = judgeAngle(ang);
        cavc::normalize(midv);
        startv = CutUnit + midv * len;
    }
    pline.addVertex(startv.x(), startv.y(), 0);
    pline.addVertex(CutUnit.x(), CutUnit.y(), 0);
    pline.isClosed() = false;

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, true);
    IsNormalLeadLine(LoopLine, temp, InLoop, true);
        
    if (leadline_error != 0)
        InLine = pline;
    else
        InLine = temp;

    return ;
}

void CreatInArcNew(cavc::Polyline<double> _pline, cavc::Vector2<double> startcut, LeadInOutInfo info, bool inloop, cavc::Polyline<double>& InLine, int& leadline_error)
{
    cavc::ClosestPoint<double> result(_pline, startcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point(), prev, midv, nextv, startv, center;
    PlineVertex<double> currentUnit, preUnit, nextUnit;
    bool InLoop = inloop, ischangeline = false, iscww = false;
    int cutNo = result.index();
    double len, startAngel, ang, angle, b = 0, radius;

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    //计算基本变量
    len = info.dLengthIn;
    if (info.dAngleIn <= 0)
        info.dAngleIn = 1.0;
    else if (info.dAngleIn >= 180)
        info.dAngleIn = 179.0;
    angle = info.dAngleIn / 180.0 * PI;
    double __ang = 0;
    if (info.dAngleIn > 90) {
        __ang = 180 - info.dAngleIn;
    }
    else if (info.dAngleIn < 90) {
        __ang = info.dAngleIn;
    }
    else {
        __ang = 90;
    }
    radius = (len / 2.0) / sin(__ang / 180.0 * PI);
    //radius = len;

    if (fabs(angle - 0) < EPS)
    {
        leadline_error = 3;
        return;
    }

    if (fabs(radius - 0) < EPS)
    {
        leadline_error = 2;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    preUnit = _pline[preno];
    currentUnit = _pline[cutNo];
    nextUnit = _pline[nextno];

    //判断起点是否是顶点且轮廓为封闭
	if (cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS) || cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS) && LoopLine.isClosed())
	{
		if (cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) {
			cutNo = nextno;
			preno = cutNo - 1;
			nextno = cutNo + 1;
			if (preno < 0)
				preno = _pline.size() - 1;
			if (nextno == _pline.size())
				nextno = 0;
			preUnit = _pline[preno];
			currentUnit = _pline[cutNo];
			nextUnit = _pline[nextno];
		}
		prev = -segTangentVector(preUnit, currentUnit, currentUnit.pos());
		if (cavc::length(prev) == 0) { leadline_error = 2;  return; }
		cavc::normalize(prev);
		nextv = segTangentVector(currentUnit, nextUnit, currentUnit.pos());
		if (cavc::length(nextv) == 0) { leadline_error = 2;  return; }
		cavc::normalize(nextv);
		//判断两向量的方向
		///180度后者0度
		if (cavc::fuzzyEqual(-prev, nextv,EPS) || cavc::fuzzyEqual(prev, nextv,EPS)) {
			ischangeline = true;
		}
		///如果两条向量带角度，则取两条向量的中点
		else {
			midv = cavc::midpoint(prev, nextv);
			///如果是内轮廓引线在外部或者外轮廓引线在内部则顺着下一个点的方向的反向引入
			if ((InLoop && !IsLineInLoop(_pline, CutUnit + midv)) ||
				(!InLoop && IsLineInLoop(_pline, CutUnit + midv)))
				midv = -nextv;
		}
		startv = CutUnit + midv * len;
	}
    //如果起点不在顶点上，或者需要该引入方向，或者为不封闭图形
    if (ischangeline || (!cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS)&&!cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) || !LoopLine.isClosed())
    {
        nextv = segTangentVector(currentUnit, nextUnit, CutUnit);
        cavc::normalize(nextv);
        midv = cavc::perp(nextv);
        cavc::normalize(midv);
        midv = -midv;
        if (InLoop)
            midv = -midv;
        if (iscww)
            midv = -midv;
        center = CutUnit + midv * radius;
        startAngel = cavc::angle(CutUnit + midv, CutUnit);
        b = tan(angle / 2);
        if (InLoop)
            if (iscww)
                ang = startAngel + 2 * angle;
            else
                ang = startAngel - 2 * angle;
        else
        {
            if (iscww)
                ang = startAngel - 2 * angle;
            else
                ang = startAngel + 2 * angle;
            b = -b;
        }
        if (iscww)
            b = -b;
        nextv = judgeAngle(ang);
        cavc::normalize(nextv);
        startv = center + nextv * radius;
    }

    double _b = 0;
    Vector2<double> _midv;
    double _ang = 0;
    if (std::abs(b) > 1) {
        if (b < 0) {
            _b = -tan((info.dAngleIn - 90) / 180.0 * PI / 2);
        }
        else {
            _b = tan((info.dAngleIn - 90) / 180.0 * PI / 2);
        }
        if (InLoop)
            if (iscww)
                _ang = startAngel + PI;
            else
                _ang = startAngel - PI;
        else
        {
            if (iscww)
                _ang = startAngel - PI;
            else
                _ang = startAngel + PI;
        }

        _midv = judgeAngle(_ang);
        cavc::normalize(_midv);
        _midv = center + _midv * radius;

        pline.addVertex(startv.x(), startv.y(), _b);
        pline.addVertex(_midv.x(), _midv.y(), b < 0 ? -1 : 1);
        pline.addVertex(CutUnit.x(), CutUnit.y(), 0);
    }
    else {
        pline.addVertex(startv.x(), startv.y(), b);
        pline.addVertex(CutUnit.x(), CutUnit.y(), 0);
    }
    pline.isClosed() = false;

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, true);
    IsNormalLeadLine(LoopLine, temp, InLoop, true);

    if (leadline_error != 0)
        InLine = pline;
    else
        InLine = temp;
     
    return ;
}

void CreatInLineArcNew(cavc::Polyline<double> _pline, cavc::Vector2<double> startcut, LeadInOutInfo info,bool inloop, cavc::Polyline<double>& InLine, int& leadline_error)
{
    cavc::ClosestPoint<double> result(_pline, startcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point(), prev, midv, nextv, start1v, start2v, center;
    PlineVertex<double> currentUnit, preUnit, nextUnit;

    bool InLoop = inloop, ischangeline = false, iscww = false;
    double startAngel, ang, angle, radius, b = 0, len;
    int cutNo = result.index();

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    //计算基本变量
    len = info.dLengthIn;
    if (info.dAngleIn <= 0)
        info.dAngleIn = 1.0;
    else if (info.dAngleIn >= 180)
        info.dAngleIn = 179.0;
    angle = info.dAngleIn / 180.0 * PI;
    radius = info.dRadiusIn;

    if (fabs(angle - 0) < EPS)
    {
        leadline_error = 3;
        return;
    }

    if (fabs(radius - 0) < EPS)
    {
        leadline_error = 2;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    preUnit = _pline[preno];
    currentUnit = _pline[cutNo];
    nextUnit = _pline[nextno];

    //判断起点是否是顶点且轮廓为封闭
    if ((cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS)|| cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) && LoopLine.isClosed()) {
        if (cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) {
            cutNo = nextno;
            preno = cutNo - 1;
            nextno = cutNo + 1;
            if (preno < 0)
                preno = _pline.size() - 1;
            if (nextno == _pline.size())
                nextno = 0;
            preUnit = _pline[preno];
            currentUnit = _pline[cutNo];
            nextUnit = _pline[nextno];
        }
        prev = -segTangentVector(preUnit, currentUnit, currentUnit.pos());
        if (cavc::length(prev) == 0) { leadline_error = 2;  return; }
        cavc::normalize(prev);
        nextv = segTangentVector(currentUnit, nextUnit, currentUnit.pos());
        if (cavc::length(nextv) == 0) { leadline_error = 2;  return; }
        cavc::normalize(nextv);
        //判断两向量的方向
        ///180度后者0度
        if (cavc::fuzzyEqual(-prev, nextv,EPS) || cavc::fuzzyEqual(prev, nextv,EPS)) {
            ischangeline = true;
        }
        ///如果两条向量带角度，则取两条向量的中点
        else {
            midv = cavc::midpoint(prev, nextv);
            ///如果是内轮廓引线在外部或者外轮廓引线在内部则顺着下一个点的方向的反向引入
            if ((InLoop && !IsLineInLoop(_pline, CutUnit + midv)) ||
                (!InLoop && IsLineInLoop(_pline, CutUnit + midv))) {
                midv = -nextv;   
            }
            start1v = CutUnit + midv * len;

            pline.addVertex(start1v.x(), start1v.y(), b);
            pline.addVertex(CutUnit.x(), CutUnit.y(), 0);
            pline.isClosed() = false;
        }    
    }
    //如果起点不在顶点上，或者需要该引入方向，或者为不封闭图形
    if (ischangeline || (!cavc::fuzzyEqual(CutUnit, _pline[cutNo].pos(),EPS)&&!cavc::fuzzyEqual(CutUnit, nextUnit.pos(),EPS)) || !LoopLine.isClosed()) {   
        nextv = segTangentVector(currentUnit, nextUnit, CutUnit);
        cavc::normalize(nextv);
        midv = cavc::perp(nextv);
        cavc::normalize(midv);
        midv = -midv;
        if (InLoop)
            midv = -midv;
        if (iscww)
            midv = -midv;
        center = CutUnit + midv * radius;
        startAngel = cavc::angle(CutUnit + midv, CutUnit);
        b = tan(angle / 4);
        if (InLoop)
            if (iscww)
                ang = startAngel +  angle;
            else
                ang = startAngel -  angle;
        else
        {
            if (iscww)
                ang = startAngel -  angle;
            else
                ang = startAngel +  angle;
            b = -b;
        }
        if (iscww)
            b = -b;
        nextv = judgeAngle(ang);
        cavc::normalize(nextv);
        start2v = center + nextv * radius;
        prev = -segTangentVector({ start2v, b }, { CutUnit,0 }, start2v);
        cavc::normalize(prev);
        start1v = start2v + prev * len;

        pline.addVertex({ start1v,0 });
        pline.addVertex({ start2v, b });
        pline.addVertex({ CutUnit,0 });
        pline.isClosed() = false;
    }

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, true);
    IsNormalLeadLine(LoopLine, temp, InLoop, true);

    if (leadline_error != 0)
        InLine = pline;
    else
        InLine = temp;

    return ;
}

void CreatOutLineNew(cavc::Polyline<double> _pline, cavc::Vector2<double> endcut, LeadInOutInfo info,bool inloop, cavc::Polyline<double>& OutLine, int& leadline_error)
{
    cavc::ClosestPoint<double> result(_pline, endcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point();
    Vector2<double> prev, midv, nextv, startv;
    PlineVertex<double> currentUnit, preUnit, nextUnit;
    bool InLoop = inloop, ischangeline = false, iscww = false;
    int cutNo = result.index();
    double len, startAngel, ang, angle;

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    len = info.dLengthOut;
    if (info.dAngleOut <= 0)
        info.dAngleOut = 1.0;
    else if (info.dAngleOut >= 180)
        info.dAngleOut = 179.0;
    angle = info.dAngleOut / 180.0 * PI;

    if (std::fabs(angle - 0) < EPS)
    {
        leadline_error = 5;
        return;
    }
    if (fabs(len - 0) < EPS)
    {
        leadline_error = 4;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    //如果作用点与轮廓顶点相同
    //引出线和引入线的方式不同，引入先是顺着结束点的前一个点和结束点形成的方向引出
    //排除cavc求最近点的误差
    if (cavc::fuzzyEqual(_pline[cutNo].pos(), CutUnit, EPS)) {
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else if (cavc::fuzzyEqual(_pline[nextno].pos(), CutUnit, EPS)) {
        cutNo = nextno;
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else {
        preUnit = _pline[cutNo];
        currentUnit = _pline[nextno];
        nextUnit = _pline[nextno];
    }
    prev = cavc::segTangentVector(preUnit, currentUnit, CutUnit);
    cavc::normalize(prev);
    startAngel = cavc::angle({ 0,0 }, prev);
    if (InLoop)
        if(iscww)
            ang = startAngel - angle;
        else
            ang = startAngel + angle;
    else
        if(iscww)
            ang = startAngel + angle;
        else
            ang = startAngel - angle;

    midv = judgeAngle(ang);
    cavc::normalize(midv);
    startv = CutUnit + midv * len;
    pline.addVertex(CutUnit.x(), CutUnit.y(), 0);
    pline.addVertex(startv.x(), startv.y(), 0);

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, false);
    IsNormalLeadLine(LoopLine, temp, InLoop, false);


    if (leadline_error != 0)
        OutLine = pline;
    else
        OutLine = temp;

    if (cavc::getPathLength(OutLine) < EPS)
        leadline_error = 4;

    return;
}

void CreatOutArcNew(cavc::Polyline<double> _pline, cavc::Vector2<double> endcut, LeadInOutInfo info,bool inloop, cavc::Polyline<double>& OutLine, int& leadline_error)
{
    cavc::ClosestPoint<double> result(_pline, endcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point();
    Vector2<double> prev, midv, nextv, startv, center;
    PlineVertex<double> currentUnit, preUnit, nextUnit;
    bool InLoop = inloop, ischangeline = false, iscww = false;
    int cutNo = result.index();
    double len, startAngel, ang, angle, radius, b = 0;

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    len = info.dLengthOut;
    if (info.dAngleOut <= 0)
        info.dAngleOut = 1.0;
    else if (info.dAngleOut >= 180)
        info.dAngleOut = 179.0;
    angle = info.dAngleOut / 180.0 * PI;
    double __ang = 0;
    if (info.dAngleOut > 90) {
        __ang = 180 - info.dAngleOut;
    }
    else if (info.dAngleOut < 90) {
        __ang = info.dAngleOut;
    }
    else {
        __ang = 90;
    }
    radius = (len / 2.0) / sin(__ang / 180.0 * PI);
    //radius = len;

    if (std::fabs(angle - 0) < EPS)
    {
        leadline_error = 5;
        return;
    }
    if (fabs(len - 0) < EPS)
    {
        leadline_error = 4;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    //如果作用点与轮廓顶点相同
    //引出线和引入线的方式不同，引入先是顺着结束点的前一个点和结束点形成的方向引出
    //排除cavc求最近点的误差
    if (cavc::fuzzyEqual(_pline[cutNo].pos(), CutUnit, EPS)) {
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else if (cavc::fuzzyEqual(_pline[nextno].pos(), CutUnit, EPS)) {
        cutNo = nextno;
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else {
        preUnit = _pline[cutNo];
        currentUnit = _pline[nextno];
        nextUnit = _pline[nextno];
    }
    prev = cavc::segTangentVector(preUnit, currentUnit, CutUnit);
    cavc::normalize(prev);
    midv = cavc::perp(prev);
    cavc::normalize(midv);
    midv = -midv;
    if (iscww)
        midv = -midv;
    if (InLoop)
        midv = -midv;
    center = CutUnit + midv * radius;
    startAngel = cavc::angle(CutUnit + midv, CutUnit);
    b = tan(angle / 2);
    if (InLoop) {
        if (iscww)
            ang = startAngel - 2 * angle;
        else
            ang = startAngel + 2 * angle;
    }
    else {
        if (iscww)
            ang = startAngel + 2 * angle;
        else
            ang = startAngel - 2 * angle;
        b = -b;
    }
    if (iscww)
        b = -b;
    nextv = judgeAngle(ang);
    cavc::normalize(nextv);
    startv = center + nextv * radius;

    double _b = 0;
    Vector2<double> _midv;
    double _ang = 0;
    if (std::abs(b) > 1) {
        if (b < 0) {
            _b = -tan((info.dAngleOut - 90) / 180.0 * PI / 2);
        }
        else {
            _b = tan((info.dAngleOut - 90) / 180.0 * PI / 2);
        }
        if (InLoop) {
            if (iscww)
                _ang = startAngel -PI;
            else
                _ang = startAngel + PI;
        }
        else {
            if (iscww)
                _ang = startAngel + PI;
            else
                _ang = startAngel - PI;
        }

        _midv = judgeAngle(_ang);
        cavc::normalize(_midv);
        _midv = center + _midv * radius;

        pline.addVertex({ CutUnit,b < 0.0 ? -1.0 : 1.0 });
        pline.addVertex({ _midv,_b });
        pline.addVertex({ startv, 0 });
    }
    else {
        pline.addVertex({ CutUnit,b });
        pline.addVertex({ startv, 0 });
    }
    pline.isClosed() = false;

    /*pline.addVertex({ CutUnit,b });
    pline.addVertex({ startv, 0 });
    pline.isClosed() = false;*/

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, false);
    IsNormalLeadLine(LoopLine, temp, InLoop, false);

    if (leadline_error != 0)
        OutLine = pline;
    else
        OutLine = temp;

    if (cavc::getPathLength(OutLine) < EPS)
        leadline_error = 4;

    return ;
}

void CreatOutLineArcNew(cavc::Polyline<double> _pline, cavc::Vector2<double> endcut, LeadInOutInfo info,bool inloop, cavc::Polyline<double>& OutLine, int& leadline_error)
{
    cavc::ClosestPoint<double> result(_pline, endcut);
    cavc::Polyline<double> LoopLine = _pline, pline;
    Vector2<double> CutUnit = result.point();
    Vector2<double> prev, midv, nextv, start1v,start2v,center;
    PlineVertex<double> currentUnit, preUnit, nextUnit;
    bool InLoop = inloop, ischangeline = false, iscww = false;
    int cutNo = result.index();
    double len, startAngel, ang, angle, radius = info.dRadiusOut, b = 0;

    //如果图形封闭且顺时针排列
    if (LoopLine.isClosed() && cavc::getArea(LoopLine) < 0)
        iscww = true;
    len = info.dLengthOut;
    if (info.dAngleOut <= 0)
        info.dAngleOut = 1.0;
    else if (info.dAngleOut >= 180)
        info.dAngleOut = 179.0;
    angle = info.dAngleOut / 180.0 * PI;

    if (std::fabs(angle - 0) < EPS)
    {
        leadline_error = 5;
        return;
    }
    if (fabs(len - 0) < EPS)
    {
        leadline_error = 4;
        return;
    }

    //找到前后点
    int preno = cutNo - 1;
    int nextno = cutNo + 1;
    if (preno < 0)
        preno = _pline.size() - 1;
    if (nextno == _pline.size())
        nextno = 0;
    //如果作用点与轮廓顶点相同
    //引出线和引入线的方式不同，引入先是顺着结束点的前一个点和结束点形成的方向引出
    //排除cavc求最近点的误差
    if (cavc::fuzzyEqual(_pline[cutNo].pos(), CutUnit, EPS)) {
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else if (cavc::fuzzyEqual(_pline[nextno].pos(), CutUnit, EPS)) {
        cutNo = nextno;
        preno = cutNo - 1;
        nextno = cutNo + 1;
        if (preno < 0)
            preno = _pline.size() - 1;
        if (nextno == _pline.size())
            nextno = 0;
        preUnit = _pline[preno];
        currentUnit = _pline[cutNo];
        nextUnit = _pline[nextno];
    }
    else {
        preUnit = _pline[cutNo];
        currentUnit = _pline[nextno];
        nextUnit = _pline[nextno];
    }
    prev = cavc::segTangentVector(preUnit, currentUnit, CutUnit);
    cavc::normalize(prev);
    midv = cavc::perp(prev);
    cavc::normalize(midv);
    midv = -midv;
    if (iscww)
        midv = -midv;
    if (InLoop)
        midv = -midv;
    center = CutUnit + midv * radius;
    startAngel = cavc::angle(CutUnit + midv, CutUnit);
    b = tan(angle / 4);
    if (InLoop) {
        if (iscww)
            ang = startAngel -  angle;
        else
            ang = startAngel +  angle;
    }
    else {
        if (iscww)
            ang = startAngel + angle;
        else
            ang = startAngel - angle;
        b = -b;
    }
    if (iscww)
        b = -b;
    nextv = judgeAngle(ang);
    cavc::normalize(nextv);
    start1v = center + nextv * radius;
    prev = segTangentVector({ CutUnit,b }, { start1v, 0 }, start1v);
    cavc::normalize(prev);
    start2v = start1v + prev * len;

    pline.addVertex({ CutUnit,b });
    pline.addVertex({ start1v, 0 });
    pline.addVertex({ start2v, 0 });
    pline.isClosed() = false;

    cavc::Polyline<double> temp = pline;
    LeadLineShrink(LoopLine, temp, InLoop, false);
    IsNormalLeadLine(LoopLine, temp, InLoop, false);

    if (leadline_error != 0)
        OutLine = pline;
    else
        OutLine = temp;

    if (cavc::getPathLength(OutLine) < EPS)
        leadline_error = 4;

    return;
}

//生成引线，返回引线状态值
int* GetLeadLine(cavc::Polyline<double> pline, cavc::Vector2<double> startcut, cavc::Vector2<double> endcut, LeadInOutInfo info, bool inloop,
    cavc::Polyline<double>& InLine,cavc::Polyline<double>& OutLine,int& leadline_error)
{
    std::cout << info.iTailInMode << endl;
    int* retstate = new int[2];
    leadline_error = 0;
    switch (info.iTailInMode)
    {
    case 0:
        break;
    case 1:
        CreatInLineNew(pline, startcut, info,inloop,InLine, leadline_error);
        break;
    case 2:
        CreatInArcNew(pline, startcut, info,inloop, InLine, leadline_error);
        break;
    case 3:
        CreatInLineArcNew(pline, startcut, info,inloop, InLine, leadline_error);
        break;
    default:
        abort();
        break;
    }
    retstate[0] = leadline_error;
    leadline_error = 0;
    switch (info.iTailOutMode)
    {
    case 0:
        break;
    case 1:
        CreatOutLineNew(pline, endcut, info,inloop, OutLine, leadline_error);
        break;
    case 2:
        CreatOutArcNew(pline, endcut, info,inloop, OutLine, leadline_error);
        break;
    case 3:
        CreatOutLineArcNew(pline, endcut, info,inloop, OutLine, leadline_error);
        break;
    default:
        abort();
        break;
    }
    retstate[1] = leadline_error;
    return retstate;
}

bool Test_DiGui(cavc::Polyline<double> pline, cavc::Polyline<double>& leadline, bool inorout) {

    cavc::PlineIntersectsResult<double> result;
    cavc::PlineIntersectsResult<double> result2;
    cavc::Vector2<double> basepoint, translatepoint;

    if (cavc::getPathLength(leadline) <= 0.1)
        return false;

    if (inorout) {
        basepoint.x() = leadline.lastVertex().pos().x();
        basepoint.y() = leadline.lastVertex().pos().y();
    }
    else {
        basepoint.x() = leadline[0].pos().x();
        basepoint.y() = leadline[0].pos().y();
    }
    for (int i = 0; i < pline.size(); i++) {
        int s = i + 1;
        if (s == pline.size())
            s = 0;
        cavc::PlineVertex<double> p1(pline[i]), p2(pline[s]);
        for (int j = 0; j < leadline.size(); j++) {
            int w = j + 1;
            if (w == leadline.size())
                continue;
            cavc::PlineVertex<double> u1(leadline[j]), u2(leadline[w]);
            cavc::IntrPlineSegsResult<double> result = cavc::intrPlineSegsWD(p1, p2, u1, u2);
            if (result.intrType != cavc::PlineSegIntrType::NoIntersect)
            {
                cavc::scalePolyline(leadline, 0.5);
                if (inorout) {
                    translatepoint.x() = leadline.lastVertex().pos().x();
                    translatepoint.y() = leadline.lastVertex().pos().y();
                }
                else {
                    translatepoint.x() = leadline[0].pos().x();
                    translatepoint.y() = leadline[0].pos().y();
                }
                cavc::translatePolyline(leadline, basepoint - translatepoint);
                if (Test_DiGui(pline, leadline, inorout))
                    return true;
            }
            else
                continue;
        }
    }
    return false;


    /*cavc::StaticSpatialIndex<double> Index = cavc::createApproxSpatialIndex(pline);
    cavc::findIntersects(pline, leadline, Index, result);
    cavc::Polyline<double> templeadline = leadline;
    cavc::invertDirection(templeadline);
    cavc::findIntersects(pline, templeadline, Index, result2);
    if (result.hasIntersects()||result2.hasIntersects()) {
        cavc::scalePolyline(leadline, 0.5);
        if (inorout) {
            translatepoint.x() = leadline.lastVertex().pos().x();
            translatepoint.y() = leadline.lastVertex().pos().y();
        }
        else {
            translatepoint.x() = leadline[0].pos().x();
            translatepoint.y() = leadline[0].pos().y();
        }
        cavc::translatePolyline(leadline, basepoint - translatepoint);
        if (Test_DiGui(pline, leadline, inorout))
            return true;
    }*/
    
}

extern "C" __declspec(dllexport) int Intersect(int size, Unit* p, bool bClose, int* sizeline, Unit* leadline, bool leadbClose, bool inorout) {

    //变换误差
    //isleadline = true;

    cavc::Polyline<double> pline, leadpline;

    for (int i = 0; i < size; i++) {
        pline.addVertex(p[i].x, p[i].y, p[i].b);
    }
    pline.isClosed() = bClose;
    for (int i = 0; i < *sizeline; i++) {
        leadpline.addVertex(leadline[i].x, leadline[i].y, leadline[i].b);
    }
    leadpline.isClosed() = leadbClose;

    Test_DiGui(pline, leadpline, inorout);

    *sizeline = leadpline.size();
    for (int i = 0; i < *sizeline; i++) {
        leadline[i].x = leadpline[i].x();
        leadline[i].y = leadpline[i].y();
        leadline[i].b = leadpline[i].bulge();
    }

    //还原误差
    //isleadline = false;

    return 0;
}



extern "C" __declspec(dllexport) int GetLeadLineShrink(int size, Unit * units,bool bclose,bool inloop, Unit startcut, Unit endcut, LeadInOutInfo info, int* insize, Unit * Inline, int* outsize, Unit * Outline, int* retstate)
{
    //变换误差
    //isleadline = true;
    //变量
    cavc::Polyline<double> InLine, OutLine;
    int leadline_error = 0;

    clock_t starttime, endtime;
    starttime = clock();
    cavc::Polyline<double> pline;
    cavc::Polyline<double> inlead;
    cavc::Polyline<double> outlead;
    int* temp = new int[2];
    cavc::Vector2<double> startcutpoint(startcut.x, startcut.y), endcutpoint(endcut.x, endcut.y);
    for (int i = 0; i < size; i++) {
        cavc::PlineVertex<double> pv(units[i].x, units[i].y, units[i].b);
        pline.addVertex(pv);
    }
    pline.isClosed() = bclose;
    pline = cavc::pruneSingularities(pline, EPS);
    temp = GetLeadLine(pline, startcutpoint, endcutpoint, info, inloop, InLine, OutLine, leadline_error);
    std::cout << "生成引线" << endl;
    *insize = InLine.size();
    *outsize = OutLine.size();
    for (int i = 0; i < InLine.size(); i++) {
        Inline[i].x = InLine[i].x();
        Inline[i].y = InLine[i].y();
        Inline[i].b = InLine[i].bulge();
        Inline[i].iFlag = 8;
    }
    for (int i = 0; i < OutLine.size(); i++) {
        Outline[i].x = OutLine[i].x();
        Outline[i].y = OutLine[i].y();
        Outline[i].b = OutLine[i].bulge();
        Outline[i].iFlag = 9;
    }
     for (int i = 0; i < 2; i++) {
        cout << temp[i] << endl;
        retstate[i] = temp[i];
    }
    endtime = clock();
    std::cout << "总耗时为：  " << endtime - starttime << endl;
    delete[] temp;

    //还原误差
    //isleadline = false;

    return 0;
}


//测试
//int main() {
//
//    char flower[10] = "rose";
//    cout << &flower[0] << "s are red\n";
//
//
//	int size = 4;
//	Unit* units = new Unit[size];
//	bool isclose = true;
//	bool inloop = false;
//	Unit sp, ep;
//	int insize=2;
//	Unit* inLine=new Unit[insize];
//	int outsize = 3;
//	Unit* outLine = new Unit[outsize];
//	int* restate = new int[2];
//	LeadInOutInfo info;
//	units[0].x = 1569.46383202342; units[0].y = 1474.63367133062; units[0].b = 0;
//	units[1].x = 2115.42042126144; units[1].y = 1474.63367133062; units[1].b = 0;
//	units[2].x = 2115.42042126144; units[2].y = 1830.45329003385; units[2].b = 0;
//	units[3].x = 1569.46383202342; units[3].y = 1830.45329003385; units[3].b = 0;
//
//    inLine[0].x = 1842.44213867188; inLine[0].y = 1479.63367133062; inLine[0].b = -1;
//    inLine[1].x = 1842.44213867188; inLine[1].y = 1474.63367133062; inLine[1].b = 0;
//    //inLine[2].x = 0; inLine[2].y = 0; inLine[2].b =0;
//
//	sp.x = 500; sp.y = 362;
//	ep.x = 500; ep.y = 362;
//	info.dAngleIn = 60;
//	info.dAngleOut = 30;
//	info.dLengthIn = 5;
//	info.dLengthOut = 5;
//	info.iOption = 0;
//	info.iPostionMode = 0;
//	info.iTailInMode = 1;
//	info.iTailOutMode = 1;
//	//GetLeadLineShrink(size, units, isclose, inloop, sp, ep, info, &insize, inLine, &outsize, outLine, restate);
//    Intersect(size, units, true, &insize, inLine, false, true);
//    
//
//	return 0;
//}
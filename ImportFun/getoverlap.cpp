#include "getoverlap.h"

//线与线
int getOverlapEntsLineandLine(PlineVertex<double> v1, PlineVertex<double> v2, PlineVertex<double> u1, PlineVertex<double> u2, double scale, std::vector<cavc::Vector2<double>>& points) {

	//double l1, l2;
	bool isclock = false;
	cavc::Polyline<double> vp, up;
	vp.addVertex(v1);
	vp.addVertex(v2);
	up.addVertex(u1);
	up.addVertex(u2);

	//判断两条线的走向是否相同，如果不同反转两条线，v20220325
	cavc::Vector2<double> vv = cavc::unitPerp(v1.pos() - v2.pos());
	cavc::Vector2<double> uv = cavc::unitPerp(u1.pos() - u2.pos());
	if (!cavc::fuzzyEqual(vv, uv,(double)0.001)) {
		cavc::invertDirection(up);
		u1 = up.vertexes()[0];
		u2 = up.lastVertex();
		isclock = true;
	}
	//判断两向量的夹角
	double angle = vv.x() * uv.x() + vv.y() * uv.y() / std::sqrt((std::pow(vv.x(), 2) + std::pow(vv.y(), 2)) * (std::pow(uv.x(), 2) + std::pow(uv.y(), 2)));
	double ang = std::acos(angle);
	double _ang = 175.0 / 180.0 * PI;
	if (ang < _ang && !cavc::utils::fuzzyEqual(ang, 0.0, 0.001))
		return 0;

	cavc::ClosestPoint<double> result1(vp, u1.pos());
	cavc::ClosestPoint<double> result2(vp, u2.pos());
	double l3 = cavc::length(result1.point() - result2.point());
	//如果找到的两个最近点相邻跳出
	if (l3 < 0.01) return 0;
	//找到一边最近点且不是端点
	if (!cavc::fuzzyEqual(result1.point(), v1.pos())) {
		//找打另一边最近点且不是端点
		if (!cavc::fuzzyEqual(result2.point(), v2.pos())) {
			//判断距离是否小于重合线的参数
			if (result1.distance() <= scale && result2.distance() <= scale) {
				//overlapents.push_back(v1);
				//overlapents.push_back(v2);
				points.push_back(u1.pos());
				points.push_back(u2.pos());
				//爹边时长边
				return 1;
			}
		}
		else {
			cavc::ClosestPoint<double> result3(up, v2.pos());
			if (result1.distance() <= scale && result3.distance() <= scale) {
				//overlapents.push_back(v1);
				//overlapents.push_back(v2);
				//爹边右点是短边
				points.push_back(u1.pos());
				points.push_back(v2.pos());
				cavc::SplitResult<double> results = splitAtPoint(u1, u2, result3.point());
				//overlapents.push_back(results.splitVertex);
				//overlapents.push_back(u2);
				return 2;
			}
		}
	}
	//找到最近点且是端点
	if (cavc::fuzzyEqual(result1.point(), v1.pos())) {
		//找到另一个最近点且不是端点
		if (!cavc::fuzzyEqual(result2.point(), v2.pos())) {
			cavc::ClosestPoint<double> result3(up, v1.pos());
			if (result3.distance() <= scale && result2.distance() <= scale) {
				//overlapents.push_back(v1);
				//overlapents.push_back(v2);
				//爹边左边时段边
				points.push_back(v1.pos());
				points.push_back(u2.pos());
				cavc::SplitResult<double> results = splitAtPoint(u1, u2, result3.point());
				//overlapents.push_back(results.updatedStart);
				//overlapents.push_back(results.splitVertex);
				return 3;
			}
		}
		//找到另一个最近点且是端点
		if (cavc::fuzzyEqual(result2.point(), v2.pos())) {
			cavc::ClosestPoint<double> result3(up, v1.pos());
			cavc::ClosestPoint<double> result4(up, v2.pos());
			if (result3.distance() <= scale && result4.distance() <= scale) {
				points.push_back(v1.pos());
				points.push_back(v2.pos());
				//爹边是短边
				//overlapents.push_back(u1);
				//overlapents.push_back(u2);
				return 4;
			}
		}
	}
	return 0;
}

//线与弧,弧与弧
int getOverlapEntsArcandArc(PlineVertex<double> v1, PlineVertex<double> v2, PlineVertex<double> u1, PlineVertex<double> u2, double scale, std::vector<cavc::Vector2<double>>& points) {

	bool isclock = false;
	cavc::Polyline<double> vp, up;
	vp.addVertex(v1);
	vp.addVertex(v2);
	up.addVertex(u1);
	up.addVertex(u2);
	if (cavc::fuzzyZero(v1.pos() - v2.pos())) return 0;
	if (cavc::fuzzyZero(u1.pos() - u2.pos())) return 0;

	//判断两条线的走向是否相同，如果不同反转两条线，v20220908
	vp.isClosed() = true;
	up.isClosed() = true;
	double vc = cavc::getArea(vp);
	double uc = cavc::getArea(up);
	if (!((vc > 0 && uc > 0) || (vc < 0 && uc < 0))) {
		cavc::invertDirection(up);
		u1 = up.vertexes()[0];
		u2 = up.lastVertex();
		isclock = true;
	}
	vp.isClosed() = false;
	up.isClosed() = false;

	bool vx = v1.x() - v2.x() < 0 ? false : true;
	bool vy = v1.y() - v2.y() < 0 ? false : true;
	bool ux = u1.x() - u2.x() < 0 ? false : true;
	bool uy = u1.y() - u2.y() < 0 ? false : true;


	//弧的半径与圆心
	cavc::ArcRadiusAndCenter<double> centerandr = cavc::arcRadiusAndCenter(v1, v2);
	cavc::ArcRadiusAndCenter<double> centerandr_ = cavc::arcRadiusAndCenter(u1, u2);
	
	

	cavc::IntrPlineSegsResult<double> _result= cavc::intrPlineSegs(v1, v2, u1, u2);

	if (_result.intrType == cavc::PlineSegIntrType::ArcOverlap) {
		points.push_back(_result.point1);
		points.push_back(_result.point2);
		return 1;
	} 

	cavc::ClosestPoint<double> result1(vp, u1.pos());
	cavc::ClosestPoint<double> result2(vp, u2.pos());
	double l3 = cavc::length(result1.point() - result2.point());
	//如果找到的两个最近点相邻跳出
	if (l3 < 1E-3) return 0;
	//如果是圆退出
	if (cavc::fuzzyEqual(result1.point(), v2.pos(),0.1) &&
		cavc::fuzzyEqual(result2.point(), v1.pos(),0.1))
		return 0;
	cavc::ClosestPoint<double> lineclosepoint(up, centerandr.center);
	cavc::ClosestPoint<double> arcclosepoint(vp, lineclosepoint.point());
	cavc::ClosestPoint<double> lengthclosepoint(up, arcclosepoint.point());
	if (!cavc::fuzzyEqual(result1.point(), v1.pos())) {
		if (!cavc::fuzzyEqual(result2.point(), v2.pos())) {
			if (result1.distance() <= scale && result2.distance() <= scale && lengthclosepoint.distance() <= scale) {
				points.push_back(u1.pos());
				points.push_back(u2.pos());
				//爹边是长边
				return 1;
			}
		}
		else {
			cavc::ClosestPoint<double> result3(up, v2.pos());
			if (result1.distance() <= scale && result3.distance() <= scale && lengthclosepoint.distance() <= scale) {
				points.push_back(u1.pos());
				points.push_back(v2.pos());
				//爹边右边是短边
				return 2;
			}
		}
	}
	if (cavc::fuzzyEqual(result1.point(), v1.pos())) {
		if (!cavc::fuzzyEqual(result2.point(), v2.pos())) {
			cavc::ClosestPoint<double> result3(up, v1.pos());
			if (result3.distance() <= scale && result2.distance() <= scale && lengthclosepoint.distance() <= scale) {
				points.push_back(v1.pos());
				points.push_back(u2.pos());
				//爹边左边是短边
				return 3;
			}
		}
		else {
			cavc::ClosestPoint<double> result3(up, v1.pos());
			cavc::ClosestPoint<double> result4(up, v2.pos());
			if (result3.distance() <= scale && result4.distance() <= scale && lengthclosepoint.distance() <= scale) {
				points.push_back(v1.pos());
				points.push_back(v2.pos());
				//爹边是短边
				return 4;
			}
		}
	}
	return 0;
}

extern "C" __declspec(dllexport) int GetIsOverLap(Unit * line1,Unit * line2,double error,Unit *point) {

	cavc::PlineVertex<double> v1(line1[0].x, line1[0].y, line1[0].b),
		v2(line1[1].x, line1[1].y, 0), u1(line2[0].x, line2[0].y, line2[0].b),
		u2(line2[1].x, line2[1].y, 0);
	std::vector<cavc::Vector2<double>> overlappoint;

	int s = 0;
	if (v1.bulgeIsZero() && u1.bulgeIsZero())
		s = getOverlapEntsLineandLine(v1, v2, u1, u2, error, overlappoint);
	else if (!v1.bulgeIsZero() && !u1.bulgeIsZero())
		s = getOverlapEntsArcandArc(v1, v2, u1, u2, error, overlappoint);
	if (s != 0) {
		point[0].x = overlappoint[0].x(); point[0].y = overlappoint[0].y();
		point[1].x = overlappoint[1].x(); point[1].y = overlappoint[1].y();  
	}
	return s;
}

extern "C" __declspec(dllexport) bool IsTwoPointSame(double* unit1,double* unit2,double error) {

	cavc::Vector2<double> v1(unit1[0], unit1[1]);
	cavc::Vector2<double> v2(unit2[0], unit2[1]);

	return cavc::fuzzyEqual(v1, v2, error);
}

extern "C" __declspec(dllexport) bool IsIntersect(int size1,Unit* unit1,bool bclose1,int size2,Unit* unit2,bool bclose2) {

	cavc::Polyline<double> pline1, pline2;

	for (int i = 0; i < size1; i++) {
		pline1.addVertex(unit1[i].x, unit1[i].y, unit1[i].b);
	}
	for (int i = 0; i < size2; i++) {
		pline2.addVertex(unit2[i].x, unit2[i].y, unit2[i].b);
	}
	pline1.isClosed() = bclose1;
	pline2.isClosed() = bclose2;

	bool isintersect = false;
	std::vector<cavc::Vector2<double>> intersectpoints;
	for (int i = 0; i < pline1.size(); i++) {
		int s = i + 1;
		if (s == pline1.size())
			continue;
		cavc::PlineVertex<double> p1(pline1[i]), p2(pline1[s]);
		for (int j = 0; j < pline2.size(); j++) {
			int w = j + 1;
			if (w == pline2.size())
				continue;
			cavc::PlineVertex<double> u1(pline2[j]), u2(pline2[w]);
			if (cavc::intrPlineSegs(p1, p2, u1, u2).intrType != cavc::PlineSegIntrType::NoIntersect)
			{
				isintersect = true;
				intersectpoints.clear();
				return isintersect;
			}
			else if (p1.bulgeIsZero() && u1.bulgeIsZero()) {
				if (getOverlapEntsLineandLine(p1, p2, u1, u2, 1E-3, intersectpoints) != 0)
				{
					isintersect = true;
					intersectpoints.clear();
					return isintersect;
				}
			}
			else {
				if (getOverlapEntsArcandArc(p1, p2, u1, u2, 1E-3, intersectpoints) != 0)
				{
					isintersect = true;
					intersectpoints.clear();
					return isintersect;
				}
			}
		}
	}
	return isintersect;
}


////测试
//int main() {
//	cavc::findIntersects
//
//	cavc::PlineVertex<double> v1(2379.0294384393933,0,	1
//	),
//		v2(1429.0294384393933,	0 ,0
//		),
//		p1(1429.0294384393933, 0, 1
//		),
//		p2(2379.0294384393933, 0,0
//		);
//
//	std::vector<cavc::Vector2<double>> points;
//	getOverlapEntsArcandArc(p1, p2, v1, v2, 1E-3, points);
//
//	return 0;
//}
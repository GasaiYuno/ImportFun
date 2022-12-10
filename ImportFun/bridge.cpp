#include "bridge.h"

int getpoint(cavc::Vector2<double> spoint, cavc::Vector2<double> epoint, double wide, cavc::Vector2<double> point, cavc::Vector2<double>& point1, cavc::Vector2<double>& point2) {

	double k = (epoint.y() - spoint.y()) / (epoint.x() - spoint.x());

	double y = sqrt(wide * wide / (k * k + 1));
	double x = k * sqrt(wide * wide / (k * k + 1));

	point1.x() = point.x() - x / 2;
	point1.y() = point.y() + y / 2;
	point2.x() = point.x() + x / 2;
	point2.y() = point.y() - y / 2;

	return 0;
}

cavc::Polyline<double> return_matrix(cavc::Vector2<double> spoint, cavc::Vector2<double> epoint, double wide) {

	Polyline<double> matrix;
	Vector2<double> point1;
	Vector2<double> point2;
	Vector2<double> point3;
	Vector2<double> point4;


	getpoint(spoint, epoint, wide, spoint, point1, point2);
	getpoint(spoint, epoint, wide, epoint, point3, point4);


	matrix.addVertex(point1.x(), point1.y(), 0);
	matrix.addVertex(point3.x(), point3.y(), 0);
	matrix.addVertex(point4.x(), point4.y(), 0);
	matrix.addVertex(point2.x(), point2.y(), 0);
	matrix.isClosed() = true;

	if (cavc::getArea(matrix) < 0)
		cavc::invertDirection(matrix);

	return matrix;
}

//判断相交点的位置,true为在轮廓里面，false为在外面
bool GetIntersectPoints(std::vector<cavc::Vector2<double>> IntersectPoints, cavc::Polyline<double> combinepline, cavc::Vector2<double>& spoint, cavc::Vector2<double>& epoint) {
	int inpoint = 0;
	int outpoint = 0;
	double slength = 999;
	double elength = 999;
	cavc::Vector2<double> smid, emid;
	double l = 999;
	map<double, cavc::Vector2<double>> intersectlist;
	std::vector<cavc::Vector2<double>> newintersectpoint;

	for (int i = 0; i < IntersectPoints.size(); i++) {
		intersectlist[cavc::length(spoint - IntersectPoints[i])] = IntersectPoints[i];
	}
	for (map<double, cavc::Vector2<double>>::iterator iter = intersectlist.begin();
		iter != intersectlist.end();
		++iter) {
		newintersectpoint.push_back(iter->second);
	}
	for (int i = 0; i < newintersectpoint.size(); i++) {
		if (i + 1 > newintersectpoint.size()) return true;
		cavc::Vector2<double> midpoint(newintersectpoint[i].x() + (newintersectpoint[i + 1].x() - newintersectpoint[i].x()) / 2,
			newintersectpoint[i].y() + (newintersectpoint[i + 1].y() - newintersectpoint[i].y()) / 2);
		double lengths = cavc::length(spoint - midpoint);
		double lengthe = cavc::length(epoint - midpoint);
		if (slength > lengths) {
			slength = lengths;
			smid = midpoint;
		}
		if (elength > lengthe) {
			elength = lengthe;
			emid = midpoint;
		}
		if (cavc::getWindingNumber(combinepline, midpoint) != 0) { i = i + 2; }
		else { outpoint++; i = i + 2; }
	}
	int s = cavc::getWindingNumber(combinepline, spoint);
	int w = cavc::getWindingNumber(combinepline, epoint);
	//起点在轮廓内且终点在轮廓外
	if (s != 0 && w == 0) epoint = emid;

	//起点在轮廓外且终点在轮廓内
	if (s == 0 && w != 0) epoint = emid;

	if (cavc::getWindingNumber(combinepline, spoint) != 0 && cavc::getWindingNumber(combinepline, epoint) != 0 && outpoint > 0) return false;
	return true;
}

//桥接轮廓
void GetBridge(int state, cavc::Polyline<double>* plines, cavc::Vector2<double> spoint, cavc::Vector2<double> epoint, double wide, cavc::CombineResult<double>& Combineresult) {

	std::vector<cavc::Vector2<double>> IntersectPoints;
	cavc::Polyline<double> line;
	line.addVertex(spoint.x(), spoint.y(), 0);
	line.addVertex(epoint.x(), epoint.y(), 0);
	line.isClosed() = false;
	//只有一个轮廓进行相交
	if (state == 1) {
		//获取交点
		cavc::StaticSpatialIndex<double> index = cavc::createApproxSpatialIndex(plines[0]);
		cavc::PlineIntersectsResult<double> result;
		cavc::findIntersects(plines[0], line, index, result);
		if (result.hasIntersects()) {
			for (auto& point : result.intersects) {
				IntersectPoints.push_back(point.pos);
			}
		}

		//判断桥接方式
		if (GetIntersectPoints(IntersectPoints, plines[0], spoint, epoint)) {
			Polyline<double> matrix = return_matrix(spoint, epoint, wide);
			Combineresult = cavc::combinePolylines(plines[0], matrix, PlineCombineMode::Exclude);
		}
		else {
			Polyline<double> matrix = return_matrix(spoint, epoint, wide);
			Combineresult = cavc::combinePolylines(plines[0], matrix, PlineCombineMode::Union);
		}
	}
	//两个轮廓进行桥接外与外或者内与内
	if (state == 2) {
		Polyline<double> matrix = return_matrix(spoint, epoint, wide);
		CombineResult<double> result1 = combinePolylines(plines[0], matrix, PlineCombineMode::Union);
		CombineResult<double> result2 = combinePolylines(plines[1], matrix, PlineCombineMode::Union);
		Combineresult = combinePolylines(result1.remaining[0], result2.remaining[0], PlineCombineMode::Union);
	}
	//两个轮廓进行桥接外与内
	if (state == 3) {
		Polyline<double> matrix = return_matrix(spoint, epoint, wide);
		CombineResult<double> result2 = combinePolylines(plines[1], matrix, PlineCombineMode::Union);
		Combineresult = combinePolylines(plines[0], result2.remaining[0], PlineCombineMode::Exclude);
	}
	return;
}

int GetUnitNums(cavc::CombineResult<double> result) {
	int count = 0;
	for (int i = 0; i < result.remaining.size();i++) {
		count += result.remaining[i].size();
	}
	for (int i = 0; i < result.subtracted.size(); i++) {
		count += result.subtracted[i].size();
	}
	return count;
}

//接口函数
extern "C" __declspec(dllexport) int BridgeNew(BridgeInfoNew& info) {
	
	cavc::CombineResult<double> Combineresult;
	cavc::Polyline<double>* plines = new cavc::Polyline<double>[info.size];
	int count = 0;
	for (int i = 0; i < info.size; i++) {
		for (int j = 0; j < info.plinesize[i]; j++) {
			plines[i].addVertex(info.pline[count].x, info.pline[count].y, info.pline[count].b);
			count++;
		}
		plines[i].isClosed() = true;
		if (cavc::getArea(plines[i]) < 0)
			cavc::invertDirection(plines[i]);
	} 
	cavc::Vector2<double> spoint(info.sx, info.sy);
	cavc::Vector2<double> epoint(info.ex, info.ey);

	//桥接
	GetBridge(info.state, plines, spoint, epoint, info.wide, Combineresult);

	//检查桥接轮廓的合理性
	for (auto& pline : Combineresult.remaining) {
		bool isnan_ = false;
		for (auto& plv : pline.vertexes()) {
			if (isnan(plv.x()) || isnan(plv.y()) || isnan(plv.bulge()))
			{
				isnan_ = true;
				break;
			}
		}
		if (isnan_) {
			Combineresult.remaining.clear();
			break;
		}
	}
	for (auto& pline : Combineresult.subtracted) {
		bool isnan_ = false;
		for (auto& plv : pline.vertexes()) {
			if (isnan(plv.x()) || isnan(plv.y()) || isnan(plv.bulge()))
			{
				isnan_ = true;
				break;
			}
		}
		if (isnan_) {
			Combineresult.subtracted.clear();
			break;
		}
	}

	int countsize = 0;
	int countunit = 0;
	*info.returnnums = Combineresult.remaining.size() + Combineresult.subtracted.size();
	info.returnplinesize = new int[*info.returnnums];
	info.returnpline = new Unit[GetUnitNums(Combineresult)];
	for (int i = 0; i < Combineresult.remaining.size(); i++) {
		info.returnplinesize[countsize] = Combineresult.remaining[i].size();
		countsize++;
		for (int j = 0; j < Combineresult.remaining[i].size(); j++) {
			info.returnpline[countunit].x = Combineresult.remaining[i][j].x();
			info.returnpline[countunit].y = Combineresult.remaining[i][j].y();
			info.returnpline[countunit].b = Combineresult.remaining[i][j].bulge();
			countunit++;
		}
	}
	for (int i = 0; i < Combineresult.subtracted.size(); i++) {
		info.returnplinesize[countsize] = Combineresult.subtracted[i].size();
		countsize++;
		for (int j = 0; j < Combineresult.subtracted[i].size(); j++) {
			info.returnpline[countunit].x = Combineresult.subtracted[i][j].x();
			info.returnpline[countunit].y = Combineresult.subtracted[i][j].y();
			info.returnpline[countunit].b = Combineresult.subtracted[i][j].bulge();
			countunit++;
		}
	}

	delete[] plines;

	return 1;
}
//内存释放
extern "C" __declspec(dllexport) void freeBridge(BridgeInfoNew info) {
	info.pline = NULL;
	info.plinesize = NULL;
	info.returnnums = NULL;
	info.returnpline = NULL;
	info.returnplinesize = NULL;

	delete[] info.pline;
	delete[] info.plinesize;
	delete info.returnnums;
	delete[] info.returnpline;
	delete[] info.returnplinesize;
}
/// <summary>
/// 求出与轮廓与线相交返回最远的两个点
/// </summary>
/// <param name="size"></param>
/// <param name="units"></param>
/// <param name="isclose"></param>
/// <param name="startpoint"></param>
/// <param name="endpoint"></param>
/// <returns>true为返回的点为相同反之为不同</returns>
extern "C" __declspec(dllexport) bool IntersectPoint(int size, Unit * units, bool isclose, Unit * startpoint, Unit * endpoint) {
	//原图信息
	cavc::Polyline<double> polyline;
	for (int i = 0; i < size; i++)
		polyline.addVertex(units[i].x, units[i].y, units[i].b);
	polyline.isClosed() = isclose;
	//线信息
	cavc::Polyline<double> line;
	line.addVertex(startpoint->x, startpoint->y, startpoint->b);
	line.addVertex(endpoint->x, endpoint->y, endpoint->b);

	cavc::StaticSpatialIndex<double> index = cavc::createApproxSpatialIndex(polyline);
	cavc::PlineIntersectsResult<double> result;
	cavc::findIntersects(polyline, line, index, result);
	if (result.hasIntersects())
	{
		double length = -1;
		startpoint->x = result.intersects[0].pos.x();
		startpoint->y = result.intersects[0].pos.y();
		startpoint->b = 0;
		cavc::Vector2<double> sp(startpoint->x, startpoint->y);
		for (auto& point : result.intersects) {
			if (cavc::length(sp - point.pos) > length) {
				endpoint->x = point.pos.x();
				endpoint->y = point.pos.y();
				endpoint->b = 0;
			}
		}
		cavc::Vector2<double> ep(endpoint->x, endpoint->y);
		if (cavc::utils::fuzzyEqual(cavc::length(sp - ep), 0.0, 0.001))
			return true;
		else
			return false;
	}
	return false;
}

int main() {

	/*cavc::Polyline<double> polyline, line;
	polyline.addVertex(0, 0, 0);
	polyline.addVertex(1, 0, 0);
	polyline.addVertex(1, 1, 0);
	polyline.addVertex(0, 1, 0);
	polyline.isClosed() = true;

	line.addVertex(0.5, -0.1, 0);
	line.addVertex(0.501, -0.1, 0);
	line.addVertex(0.501, 1.1, 0);
	line.addVertex(0.5, 1.1, 0);
	line.isClosed() = true;

	cavc::CombineResult<double> result = cavc::combinePolylines(polyline, line, PlineCombineMode::Exclude);*/

	std::cout << std::acos(0.5) << std::endl;

	int s = 0;

	return 0;
}

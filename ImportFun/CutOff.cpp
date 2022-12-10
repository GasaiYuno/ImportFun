#include "CutOff.h"

std::vector<cavc::PlineIntersect<double>> GetIntersectpoints(cavc::Polyline<double> polyline, cavc::Polyline<double> line) {

	cavc::StaticSpatialIndex<double> tree = cavc::createApproxSpatialIndex(polyline);
	cavc::PlineIntersectsResult<double> intersectresult;
	cavc::findIntersects(polyline, line, tree, intersectresult);
	std::vector<cavc::PlineIntersect<double>> intersectpoints = intersectresult.intersects;

	return intersectpoints;
}
cavc::Polyline<double>* CutLoop(cavc::Polyline<double> polyline, std::vector<cavc::PlineIntersect<double>> intersectpoints, int& polylinessize) {

	cavc::Polyline<double> temppolyline = polyline;

	vector<int> intersectindex;
	vector<vector<int>> intersectindex_;

	if (intersectpoints.size() % 2 != 0)
		intersectpoints.pop_back();

	std::sort(intersectpoints.begin(), intersectpoints.end(), indexcomp());

	for (auto& intersectpoint : intersectpoints) {
		cavc::ClosestPoint<double> close(temppolyline, intersectpoint.pos);
		int nextindex = close.index() + 1 == temppolyline.size() ? 0 : close.index() + 1;
		cavc::SplitResult<double> result = cavc::splitAtPoint(temppolyline[close.index()], temppolyline[nextindex], intersectpoint.pos);
		temppolyline[close.index()] = result.updatedStart;
		if (nextindex == 0)
			temppolyline.vertexes().push_back(result.splitVertex);
		else
			temppolyline.vertexes().insert(temppolyline.vertexes().begin() + nextindex, result.splitVertex);
		int nextindex_ = close.index() + 1 == temppolyline.size() ? 0 : close.index() + 1;
		intersectindex.push_back(nextindex_);
	}
	int i = 0;
	vector<int> couple;
	couple.reserve(2);
	couple.push_back(0);
	couple.push_back(0);
	for (auto& index : intersectindex) {
		if (i == 0) {
			couple[i] = index;
		}
		if (i == 1) {
			couple[i] = index;
			intersectindex_.push_back(couple);
			i = 0;
			continue;
		}
		i++;
	}

	polylinessize = intersectpoints.size() / 2 + 1;
	cavc::Polyline<double>* newpolylines = new cavc::Polyline<double>[polylinessize];

	int s = 0;
	for (auto& index_ : intersectindex_) {
		cavc::Polyline<double> newpolyline;
		for (int i = index_[0]; i <= index_[1]; i++) {
			newpolyline.addVertex(temppolyline[i]);
		}
		newpolyline.lastVertex().bulge() = 0;
		newpolyline.isClosed() = true;
		newpolylines[s] = newpolyline;
		s++;
	}

	cavc::Polyline<double> newpolyline_;
	set<int> noneedindex;
	for (auto& index_ : intersectindex_) {
		for (int i = index_[0]+1; i < index_[1]; i++)
			noneedindex.insert(i);
	}
	for (int w = 0; w < temppolyline.size(); w++) {
		if (noneedindex.count(w)) {
			newpolyline_.lastVertex().bulge() = 0;
			continue;
		}
		else
			newpolyline_.addVertex(temppolyline[w]);
	}
	newpolyline_.isClosed() = true;
	newpolylines[polylinessize - 1] = newpolyline_;
	return newpolylines;
}


/// <summary>
/// 裁断封闭轮廓
/// </summary>
/// <param name="size">轮廓图元个数</param>
/// <param name="loop">轮廓图元</param>
/// <param name="line">裁断线的图元</param>
extern "C" __declspec(dllexport) void CutOff(int size, Unit* loop, Unit* line,int *rpolylinesize,int *runitsize,Unit* runits) {

	cavc::Polyline<double> polyline, cutline;
	for (int i = 0; i < size; i++)
	{
		polyline.addVertex(loop[i].x, loop[i].y, loop[i].b);
	}
	cutline.addVertex(line[0].x, line[0].y, line[0].b);
	cutline.addVertex(line[1].x, line[1].y, line[1].b);
	polyline.isClosed() = true;

	int newpolylinesize = 0;
	cavc::Polyline<double>* newpolylines;
	newpolylines = CutLoop(polyline, GetIntersectpoints(polyline, cutline), newpolylinesize);

	*rpolylinesize = newpolylinesize;
	int rcount = 0;
	for (int i = 0; i < newpolylinesize; i++)
	{
		for (int s = 0; s < newpolylines[i].size(); s++) {
			runits[rcount + s].x = newpolylines[i][s].x();
			runits[rcount + s].y = newpolylines[i][s].y();
			runits[rcount + s].b = newpolylines[i][s].bulge();
		}
		runitsize[i] = newpolylines[i].size();
		rcount += newpolylines[i].size();
	}
	std::cout << "返回" << *rpolylinesize << "个轮廓" << std::endl;
	int rcount_ = 0;
	for (int i = 0; i < *rpolylinesize; i++) {
		std::cout << "第" << i << "个轮廓信息：" << std::endl;
		for (int s = 0; s < runitsize[i]; s++) {
			std::cout << "x:" << runits[rcount_ + s].x << "y:" << runits[rcount_ + s].y << "b:" << runits[rcount_ + s].b << std::endl;
		}
		rcount_ += runitsize[i];
	}
	newpolylines = NULL;
	delete newpolylines;
}
/// <summary>
/// 求三点的夹角
/// </summary>
/// <param name="beforeunit"></param>  
/// <param name="persentunit"></param>
/// <param name="nextunit"></param>
/// <returns></returns>
extern "C" __declspec(dllexport) double GetCorner(Unit beforeunit,Unit persentunit,Unit nextunit) {

	cavc::PlineVertex<double> bv(beforeunit.x, beforeunit.y, beforeunit.b);
	cavc::PlineVertex<double> pv(persentunit.x, persentunit.y, persentunit.b);
	cavc::PlineVertex<double> nv(nextunit.x, nextunit.y, nextunit.b);

	cavc::Vector2<double> bv2 = cavc::segTangentVector(bv, pv, pv.pos());
	cavc::normalize(bv2);
	cavc::Vector2<double> nv2 = -cavc::segTangentVector(pv, nv, pv.pos());
	cavc::normalize(nv2);

	double angle = bv2.x() * nv2.x() + bv2.y() * nv2.y() / std::sqrt((std::pow(bv2.x(), 2) + std::pow(bv2.y(), 2)) * (std::pow(nv2.x(), 2) + std::pow(nv2.y(), 2)));

	return std::acos(angle);
}


extern "C" __declspec(dllexport) void GetLoopCentroid(int size,Unit* units, double* cx, double* cy)
//求质心
{
	Polyline<double>loop;
	for (int i = 0; i < size; i++)
		loop.addVertex(units[i].x, units[i].y, units[i].b);
	loop.isClosed() = true;
	
	int n;
	double x, y; //质心
	double totalArea = 0, eachArea = 0; //多边形总面积和各三角形面积

	x = y = 0.0;
	n = loop.size();
	for (int i = 0; i < n; i++) {
		eachArea = (loop[i].x() * loop[(i + 1) % n].y() - loop[(i + 1) % n].x() * loop[i].y()) / 2.0;//两个顶点与原点构成的三角形
		x += (loop[i].x() + loop[(i + 1) % n].x()) * eachArea; //每一步质心的除以3放在最后处理
		y += (loop[i].y() + loop[(i + 1) % n].y()) * eachArea;
		totalArea += eachArea; //总面积
	}

	*cx = x / (totalArea * 3.0);//求的多边形质心
	*cy = y / (totalArea * 3.0);

	return;
}


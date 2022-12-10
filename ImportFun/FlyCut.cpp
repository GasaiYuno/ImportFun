#include "FlyCut.h"

void LeftDown(cavc::Vector2<double> _basePoint, int i, double error, double flycutspace,bool isbox,
	std::vector<cavc::Vector2<double>> BasePoint, std::map<double, std::vector<Line>> LineList,
	std::vector<Line>& LastLine) {

	//变量
	LastLine.clear();

	//先去掉水平和垂直的线段
	auto LineList_ = LineList;
	if (LineList_.count(0.0) == 1)
		LineList_.erase(0.0);
	if (LineList_.count(90.0) == 1)
		LineList_.erase(90.0);

	std::map<double, std::vector<Line>> TrueList;
	std::map<double, std::vector<Line>> TrueListk;
	//遍历次数
	int count = 0;
	//线路容差
	double lineerror = error;
	//方向容差
	double directionerror = flycutspace;
	//基点
	cavc::Vector2<double> basePoint = _basePoint;

	//左下,右上
	int list = 0;
	if (i == 0 || i == 2) {
		//反向遍历
		for (std::map<double, std::vector<Line>>::reverse_iterator iter = LineList_.rbegin(); iter != LineList_.rend(); iter++) {
			TrueList[list] = iter->second;
			list++;
		}
	}
	//右下，左上
	if (i == 1 || i == 3) {
		TrueList = LineList_;
	}
	//水平和垂直最后切
	if (LineList.count(0.0) == 1) {
		TrueListk[0] = LineList.find(0.0)->second;
	}
	if (LineList.count(90.0) == 1) {
		TrueListk[1] = LineList.find(90.0)->second;
	}
	//先排序一般线段
	for (auto& l : TrueList) {
		double langle = l.second[0].ka;
		int state = i;
		//找到最后一个点的最近基点
		if (LastLine.size() > 0) {
			double plength = 999;
			for (int i = 0; i < BasePoint.size(); i++) {
				if (plength > cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P1)) {
					plength = cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P1);
					basePoint = BasePoint[i];
					state = i;
				}
			}
			//确定下一个基点
			int index1 = state + 1;
			int index2 = state - 1;
			if (index1 == 4)
				index1 = 0;
			if (index2 == -1)
				index2 = 3;
			state = cavc::length(LastLine[LastLine.size() - 1].P2 - BasePoint[index1]) <
				cavc::length(LastLine[LastLine.size() - 1].P2 - BasePoint[index2]) ? index1 : index2;
			basePoint = BasePoint[state];
		}
		//算出线段距离基点的相对距离
		std::map<double, std::vector<Line>> kLineList;
		for (auto& _l : l.second) {
			//基点到线条起点的向量
			cavc::Vector2<double> svector = basePoint - _l.P1;
			//基点到线条终点的向量
			cavc::Vector2<double> evector = basePoint - _l.P2;
			//基点到起点的角度
			_l.sangle = std::atan2(svector.y(), svector.x()) * 180 / cavc::utils::pi<double>();
			//基点到终点的角度
			_l.eangle = std::atan2(evector.y(), evector.x()) * 180 / cavc::utils::pi<double>();
			//基点到起点的距离
			double length = cavc::length(svector);
			if (_l.sangle - _l.eangle >= 0 || cavc::utils::fuzzyEqual(_l.sangle, 0.0, 0.01)) {
				_l.mode = true;
				if (cavc::utils::fuzzyEqual(_l.sangle, 0.0, 0.1)) {
					_l.sangle = _l.eangle;
				}
			}
			///水平直线
			if (langle == 0) {
				_l.length = _l.P1.y() - basePoint.y();
			}
			///垂直直线
			else if (langle == 90) {
				_l.length = _l.P1.x() - basePoint.x();
			}
			//斜线
			else {
				double l1 = basePoint.y() - std::tan(langle * cavc::utils::pi<double>() / 180.0) * basePoint.x();
				double l2 = _l.P1.y() - std::tan(langle * cavc::utils::pi<double>() / 180.0) * _l.P1.x();
				_l.length = (int)(l1 - l2) * std::cos(langle * cavc::utils::pi<double>() / 180.0);
			}
			if (kLineList.size() == 0) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			//找到相同的y值
			else if (kLineList.count(_l.length) == 1) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			else {
				int state = 0;
				for (auto& key : kLineList) {
					//找到存在误差的线
					if (std::abs(key.first - std::fabs(_l.length)) <= 0.1) {
						kLineList[key.first].push_back(_l);
						state = 1;
						break;
					}
				}
				//没有找到符合误差的线
				if (state == 0) {
					kLineList[std::fabs(_l.length)].push_back(_l);
				}
			}
		}
		int a = 0;
		cavc::Vector2<double> directionv;
		//通过相对距离排序
		for (auto& kl : kLineList) {
			//判断下一个距离上的线段的方向是否与前一个距离上的线段关系，相同就反序，不相同就默认
			if (a != 0) {
				cavc::Vector2<double> kdirectionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);
				//如果方向相同就反向
				if (cavc::fuzzyEqual(kdirectionv, directionv, 0.01)) {
					for (auto& al : kl.second) {
						auto tl = al.P1;
						al.P1 = al.P2;
						al.P2 = tl;
					}
				}
			}
			//开始通过方向排序
			switch (state)
			{
			case 0:
				if (kl.second[0].P2.x() - kl.second[0].P1.x() > 0) {
					std::sort(kl.second.begin(), kl.second.end(), minxcomp());
				}
				else {
					std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
				}
				break;
			case 1:
				if (kl.second[0].P2.y() - kl.second[0].P1.y() > 0) {
					std::sort(kl.second.begin(), kl.second.end(), minycomp());
				}
				else {
					std::sort(kl.second.begin(), kl.second.end(), maxycomp());
				}
				break;
			case 2:
				if (kl.second[0].P2.x() - kl.second[0].P1.x() < 0) {
					std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
				}
				else {
					std::sort(kl.second.begin(), kl.second.end(), minxcomp());
				}
				break;
			case 3:
				if (kl.second[0].P2.y() - kl.second[0].P1.y() < 0) {
					std::sort(kl.second.begin(), kl.second.end(), maxycomp());
				}
				else {
					std::sort(kl.second.begin(), kl.second.end(), minycomp());
				}
				break;
			default:
				break;
			}
			//将排序的内容写入最终表中
			for (auto& al : kl.second) {
				LastLine.push_back(al);
			}

			directionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);

			a++;
		}
		count++;
	}
	//水平和垂直线段排序基本和上边一样
	int countk = count;
	basePoint = BasePoint[0];
	for (auto& l : TrueListk) {
		double langle = l.second[0].ka;
		int state = 0;
		if ((LastLine.size() > 0 && countk != 0)) {
			double plength = 999;
			for (int i = 0; i < BasePoint.size(); i++) {
				if (plength > cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P2)) {
					plength = cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P2);
					basePoint = BasePoint[i];
					state = i;
				}
			}
			basePoint = BasePoint[state];
		}
		else if (isbox) {
			basePoint = BasePoint[i];
		}
		std::map<double, std::vector<Line>> kLineList;
		for (auto& _l : l.second) {
			//基点到线条起点的向量
			cavc::Vector2<double> svector = basePoint - _l.P1;
			//基点到线条终点的向量
			cavc::Vector2<double> evector = basePoint - _l.P2;
			//基点到起点的角度
			_l.sangle = std::atan2(svector.y(), svector.x()) * 180 / cavc::utils::pi<double>();
			//基点到终点的角度
			_l.eangle = std::atan2(evector.y(), evector.x()) * 180 / cavc::utils::pi<double>();
			//基点到起点的距离
			///水平直线
			if (langle == 0) {
				_l.length = _l.P1.y() - basePoint.y();
			}
			///垂直直线
			else if (langle == 90) {
				_l.length = _l.P1.x() - basePoint.x();
			}
			if (kLineList.size() == 0) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			//找到相同的y值
			else if (kLineList.count(_l.length) == 1) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			else {
				int state = 0;
				for (auto& key : kLineList) {
					//找到存在误差的线
					if (std::abs(key.first - std::fabs(_l.length)) <= 0.1) {
						kLineList[key.first].push_back(_l);
						state = 1;
						break;
					}
				}
				//没有找到符合误差的线
				if (state == 0) {
					kLineList[std::fabs(_l.length)].push_back(_l);
				}
			}
		}
		int a = 0;
		Line* line = &kLineList.begin()->second[0];
		cavc::Vector2<double> directionv = cavc::unitPerp(basePoint - line->P1);
		for (auto& kl : kLineList)
		{
			if (a != 0) {
				cavc::Vector2<double> kdirectionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);
				//如果方向相同就反向
				if (cavc::fuzzyEqual(kdirectionv, directionv, 0.01)) {
					for (auto& al : kl.second) {
						auto tl = al.P1;
						al.P1 = al.P2;
						al.P2 = tl;
					}
				}
			}
			if (a == 0) {
				cavc::Vector2<double> kdirectionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);
				//如果方向不同就反向
				if (!cavc::fuzzyEqual(kdirectionv, directionv, 0.01)) {
					for (auto& al : kl.second) {
						auto tl = al.P1;
						al.P1 = al.P2;
						al.P2 = tl;
					}
				}
			}

			switch (state)
			{
			case 0:
				if (kl.second[0].P2.x() - kl.second[0].P1.x() > 0 && !cavc::utils::fuzzyEqual(kl.second[0].P2.x(), kl.second[0].P1.x(), 0.0001)) {
					std::sort(kl.second.begin(), kl.second.end(), minxcomp());
				}
				else {
					//如果x相等比较y
					if (cavc::utils::fuzzyEqual(kl.second[0].P2.x(), kl.second[0].P1.x(), 0.0001)) {
						if (kl.second[0].P2.y() - kl.second[0].P1.y() > 0)
							std::sort(kl.second.begin(), kl.second.end(), minycomp());
						else
							std::sort(kl.second.begin(), kl.second.end(), maxycomp());
					}
					else
						std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
				}
				break;
			case 1:
				if (kl.second[0].P2.y() - kl.second[0].P1.y() > 0 && !cavc::utils::fuzzyEqual(kl.second[0].P2.y(), kl.second[0].P1.y(), 0.0001)) {
					std::sort(kl.second.begin(), kl.second.end(), minycomp());
				}
				else {
					//如果y相等比较x
					if (cavc::utils::fuzzyEqual(kl.second[0].P2.y(), kl.second[0].P1.y(), 0.0001)) {
						if (kl.second[0].P2.x() - kl.second[0].P1.x() > 0)
							std::sort(kl.second.begin(), kl.second.end(), minxcomp());
						else
							std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
					}
					else
						std::sort(kl.second.begin(), kl.second.end(), maxycomp());
				}
				break;
			case 2:
				if (kl.second[0].P2.x() - kl.second[0].P1.x() < 0 && !cavc::utils::fuzzyEqual(kl.second[0].P2.x(), kl.second[0].P1.x(), 0.0001)) {
					std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
				}
				else {
					//如果x相等比较y
					if (cavc::utils::fuzzyEqual(kl.second[0].P2.x(), kl.second[0].P1.x(), 0.0001)) {
						if (kl.second[0].P2.y() - kl.second[0].P1.y() < 0)
							std::sort(kl.second.begin(), kl.second.end(), maxycomp());
						else
							std::sort(kl.second.begin(), kl.second.end(), minycomp());
					}
					else
						std::sort(kl.second.begin(), kl.second.end(), minxcomp());
				}
				break;
			case 3:
				if (kl.second[0].P2.y() - kl.second[0].P1.y() < 0 && !cavc::utils::fuzzyEqual(kl.second[0].P2.y(), kl.second[0].P1.y(), 0.0001)) {
					std::sort(kl.second.begin(), kl.second.end(), maxycomp());
				}
				else {
					//如果y相等比较x
					if (cavc::utils::fuzzyEqual(kl.second[0].P2.y(), kl.second[0].P1.y(), 0.0001)) {
						if (kl.second[0].P2.x() - kl.second[0].P1.x() < 0)
							std::sort(kl.second.begin(), kl.second.end(), maxxcomp());
						else
							std::sort(kl.second.begin(), kl.second.end(), minxcomp());
					}
					else
						std::sort(kl.second.begin(), kl.second.end(), minycomp());
				}
				break;
			default:
				break;
			}

			for (auto& al : kl.second) {
				LastLine.push_back(al);
			}

			directionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);

			a++;
		}
		countk++;
	}
}

void LoopToLine(FlyCutInfo info, std::vector<cavc::Vector2<double>> BasePoint, 
	std::vector<cavc::Polyline<double>> PList, std::vector<Line>& LastLine) {
	
	//变量
	std::map<double, std::vector<Line>> LineList;
	std::map<double, std::vector<Line>> BoxList;
	std::vector<Line> boxlastline;
	cavc::Polyline<double> boxpolyline;
	LineList.clear();

	boxpolyline.addVertex(BasePoint[0].x(), BasePoint[0].y(), 0);
	boxpolyline.addVertex(BasePoint[1].x(), BasePoint[1].y(), 0);
	boxpolyline.addVertex(BasePoint[2].x(), BasePoint[2].y(), 0);
	boxpolyline.addVertex(BasePoint[3].x(), BasePoint[3].y(), 0);
	boxpolyline.isClosed() = true;

	cavc::StaticSpatialIndex<double> tree = cavc::createApproxSpatialIndex(boxpolyline);

	for (auto& pline : PList) {
		//拆分为线条
		for (int i = 0; i < pline.size(); i++) {
			int s = i + 1;
			if (s == pline.size() && pline.isClosed())
				s = 0;
			else if (s == pline.size() && !pline.isClosed())
				continue;
			Line l;
			l.P1 = pline.vertexes()[i].pos();
			l.P2 = pline.vertexes()[s].pos();
			cavc::Vector2<double> KP = l.P2 - l.P1;
			double Ka = std::atan2(KP.y(), KP.x()) * 180 / cavc::utils::pi<double>();
			if (Ka < 0)
				Ka = Ka + 180;
			else if (Ka == 180)
				Ka = 0;
			int iKa = round(Ka);
			l.ka = iKa;
			//是否在包路和上的线段
			if (iKa == 0.00 || iKa == 90.00) {
				cavc::PlineIntersectsResult<double> result;
				cavc::Polyline<double> linep;
				linep.addVertex(pline.vertexes()[i]);
				linep.addVertex(pline.vertexes()[s]);
				linep.isClosed() = false;
				cavc::findIntersects(boxpolyline, linep, tree, result);
				if (result.coincidentIntersects.size() != 0)
					BoxList[iKa].push_back(l);
				else
					LineList[iKa].push_back(l);
			}	
			else
				LineList[iKa].push_back(l);
		}
	}

	//如果不是边框
	if (BoxList.count(0.00) == 0) {
		if (BoxList[90.00].size() != 0) {
			LineList[90.00].insert(LineList[90.00].end(), BoxList[90.00].begin(), BoxList[90.00].end());
			BoxList[90.00].clear();
		}
	}
	if (BoxList.count(90.00) == 0) {
		if (BoxList[0.00].size() != 0) {
			LineList[0.00].insert(LineList[0.00].end(), BoxList[0.00].begin(), BoxList[0.00].end());
			BoxList[0.00].clear();
		}
	}

	if (LineList.size() == 2 && LineList.count(0.00) == 1 && LineList.count(90.00) == 1) {
		if (BoxList[90.00].size() != 0) {
			LineList[90.00].insert(LineList[90.00].end(), BoxList[90.00].begin(), BoxList[90.00].end());
			BoxList[90.00].clear();
		}
		if (BoxList[0.00].size() != 0) {
			LineList[0.00].insert(LineList[0.00].end(), BoxList[0.00].begin(), BoxList[0.00].end());
			BoxList[0.00].clear();
		} 
	}

	switch (info.startcutlocality)
	{
		//默认
	case 0:
		LeftDown(BasePoint[0], 0, info.error, info.flycutspace,false, BasePoint, LineList, LastLine);
		break;
		//左下
	case 1:
		LeftDown(BasePoint[0], 0, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//右下
	case 2:
		LeftDown(BasePoint[1], 1, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//左上
	case 3:
		LeftDown(BasePoint[3], 3, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//右上
	case 4:
		LeftDown(BasePoint[2], 2, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
	default:
		break;
	}

	//把边框添加到最后
	if (BoxList.size() != 0) {
		if (BoxList[0.00].size() != 0 || BoxList[90.00].size() != 0) {
			cavc::ClosestPoint<double> close(boxpolyline, LastLine[LastLine.size() - 1].P2);
			double length1 = cavc::length(boxpolyline[close.index()].pos() - LastLine[LastLine.size() - 1].P2);
			int nextindex = close.index() + 1 == boxpolyline.size() ? 0 : close.index() + 1;
			double length2 = cavc::length(boxpolyline[nextindex].pos() - LastLine[LastLine.size() - 1].P2);
			if (length1 <= length2)
				LeftDown(BasePoint[close.index()], close.index(), info.error, info.flycutspace, true, BasePoint, BoxList, boxlastline);
			else
				LeftDown(BasePoint[nextindex], nextindex, info.error, info.flycutspace, true, BasePoint, BoxList, boxlastline);
			int count1 = BoxList[0.00].size();
			int count2 = BoxList[90.00].size();
			std::vector<Line> newlastline;
			for (int i = 0; i < count1 / 2; i++) 
				newlastline.push_back(boxlastline[i]);
			for (int i = count1 + count2 / 2; i < count1 + count2; i++)
				newlastline.push_back(boxlastline[i]);
			for (int i = count1 / 2; i < count1; i++)
				newlastline.push_back(boxlastline[i]);
			for (int i = count1; i < count1 + count2 / 2; i++)
				newlastline.push_back(boxlastline[i]);
			LastLine.insert(LastLine.end(), newlastline.begin(), newlastline.end());
		}
	}	
}

extern "C" __declspec(dllexport) int FlyCutLine(FlyCutInfo info, FlyCutInformation & fcinfo)
{
	//变量
	std::vector<cavc::Polyline<double>> PList;
	std::vector<cavc::Vector2<double>> BasePoint;
	std::vector<Line> LastLine;
	std::vector<std::vector<cavc::PlineVertex<double>>> RemainList;

	int count = 0;
	cavc::AABB<double> maxmin;
	maxmin.xMin = 999;
	maxmin.yMin = 999;
	maxmin.xMax = -999;
	maxmin.yMax = -999;
	
	//把轮廓信息加入轮廓列表中
	int allcount = 0;
	for (int i = 0; i < fcinfo.size1; i++) {
		allcount = fcinfo.size2[i] + allcount;
		bool isb = false;
		cavc::Polyline<double> pline;
		for (int j = 0; j < fcinfo.size2[i]; j++) {
			pline.addVertex(fcinfo.units[count].x, fcinfo.units[count].y, fcinfo.units[count].b);
			count++;
		}
		if (!isb) {
			pline.isClosed() = fcinfo.closelist[i];
			PList.push_back(pline);
		}

		//获取轮廓的包路和
		cavc::AABB<double> pmaxmin = cavc::getExtents(pline);
		if (pmaxmin.xMin < maxmin.xMin)
			maxmin.xMin = pmaxmin.xMin;
		if (pmaxmin.xMax > maxmin.xMax)
			maxmin.xMax = pmaxmin.xMax;
		if (pmaxmin.yMin < maxmin.yMin)
			maxmin.yMin = pmaxmin.yMin;
		if (pmaxmin.yMax > maxmin.yMax)
			maxmin.yMax = pmaxmin.yMax;
	}

	//把包路和填入基点中
	cavc::Vector2<double> point1(maxmin.xMin, maxmin.yMin);
	cavc::Vector2<double> point2(maxmin.xMax, maxmin.yMin);
	cavc::Vector2<double> point3(maxmin.xMax, maxmin.yMax);
	cavc::Vector2<double> point4(maxmin.xMin, maxmin.yMax);
	BasePoint.push_back(point1);
	BasePoint.push_back(point2);
	BasePoint.push_back(point3);
	BasePoint.push_back(point4);

	LoopToLine(info, BasePoint, PList, LastLine);

	*fcinfo.returnsize = LastLine.size() * 2;
	fcinfo.returnunits = new Unit[LastLine.size() * 2];
	int count2 = 0;
	
	for (auto& line : LastLine) {
		fcinfo.returnunits[count2].x = line.P1.x();
		fcinfo.returnunits[count2].y = line.P1.y();
		fcinfo.returnunits[count2].b = 0;
		fcinfo.returnunits[count2].iFlag = 0;
		fcinfo.returnunits[count2 + 1].x = line.P2.x();
		fcinfo.returnunits[count2 + 1].y = line.P2.y();
		fcinfo.returnunits[count2 + 1].b = 0;
		fcinfo.returnunits[count2 + 1].iFlag = 3;
		count2 = count2 + 2;
	}

	*fcinfo.remainsize = RemainList.size() * 2;
	fcinfo.remainunits = new Unit[RemainList.size() * 2];
	int count3 = 0;

	for (auto& remain : RemainList) {
		fcinfo.remainunits[count3].x = remain[0].x();
		fcinfo.remainunits[count3].y = remain[0].y();
		fcinfo.remainunits[count3].b = remain[0].bulge();
		fcinfo.remainunits[count3].iFlag = 0;
		fcinfo.remainunits[count3 + 1].x = remain[1].x();
		fcinfo.remainunits[count3 + 1].y = remain[1].y();
		fcinfo.remainunits[count3 + 1].b = remain[1].bulge();
		fcinfo.remainunits[count3 + 1].iFlag = 3;
		count3 = count3 + 2;
	}


	PList.clear();
	LastLine.clear();
	BasePoint.clear();
	RemainList.clear();

	return 0;
}
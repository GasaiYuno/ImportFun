#include "FlyCut.h"

void LeftDown(cavc::Vector2<double> _basePoint, int i, double error, double flycutspace,bool isbox,
	std::vector<cavc::Vector2<double>> BasePoint, std::map<double, std::vector<Line>> LineList,
	std::vector<Line>& LastLine) {

	//����
	LastLine.clear();

	//��ȥ��ˮƽ�ʹ�ֱ���߶�
	auto LineList_ = LineList;
	if (LineList_.count(0.0) == 1)
		LineList_.erase(0.0);
	if (LineList_.count(90.0) == 1)
		LineList_.erase(90.0);

	std::map<double, std::vector<Line>> TrueList;
	std::map<double, std::vector<Line>> TrueListk;
	//��������
	int count = 0;
	//��·�ݲ�
	double lineerror = error;
	//�����ݲ�
	double directionerror = flycutspace;
	//����
	cavc::Vector2<double> basePoint = _basePoint;

	//����,����
	int list = 0;
	if (i == 0 || i == 2) {
		//�������
		for (std::map<double, std::vector<Line>>::reverse_iterator iter = LineList_.rbegin(); iter != LineList_.rend(); iter++) {
			TrueList[list] = iter->second;
			list++;
		}
	}
	//���£�����
	if (i == 1 || i == 3) {
		TrueList = LineList_;
	}
	//ˮƽ�ʹ�ֱ�����
	if (LineList.count(0.0) == 1) {
		TrueListk[0] = LineList.find(0.0)->second;
	}
	if (LineList.count(90.0) == 1) {
		TrueListk[1] = LineList.find(90.0)->second;
	}
	//������һ���߶�
	for (auto& l : TrueList) {
		double langle = l.second[0].ka;
		int state = i;
		//�ҵ����һ������������
		if (LastLine.size() > 0) {
			double plength = 999;
			for (int i = 0; i < BasePoint.size(); i++) {
				if (plength > cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P1)) {
					plength = cavc::length(BasePoint[i] - LastLine[LastLine.size() - 1].P1);
					basePoint = BasePoint[i];
					state = i;
				}
			}
			//ȷ����һ������
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
		//����߶ξ���������Ծ���
		std::map<double, std::vector<Line>> kLineList;
		for (auto& _l : l.second) {
			//���㵽������������
			cavc::Vector2<double> svector = basePoint - _l.P1;
			//���㵽�����յ������
			cavc::Vector2<double> evector = basePoint - _l.P2;
			//���㵽���ĽǶ�
			_l.sangle = std::atan2(svector.y(), svector.x()) * 180 / cavc::utils::pi<double>();
			//���㵽�յ�ĽǶ�
			_l.eangle = std::atan2(evector.y(), evector.x()) * 180 / cavc::utils::pi<double>();
			//���㵽���ľ���
			double length = cavc::length(svector);
			if (_l.sangle - _l.eangle >= 0 || cavc::utils::fuzzyEqual(_l.sangle, 0.0, 0.01)) {
				_l.mode = true;
				if (cavc::utils::fuzzyEqual(_l.sangle, 0.0, 0.1)) {
					_l.sangle = _l.eangle;
				}
			}
			///ˮƽֱ��
			if (langle == 0) {
				_l.length = _l.P1.y() - basePoint.y();
			}
			///��ֱֱ��
			else if (langle == 90) {
				_l.length = _l.P1.x() - basePoint.x();
			}
			//б��
			else {
				double l1 = basePoint.y() - std::tan(langle * cavc::utils::pi<double>() / 180.0) * basePoint.x();
				double l2 = _l.P1.y() - std::tan(langle * cavc::utils::pi<double>() / 180.0) * _l.P1.x();
				_l.length = (int)(l1 - l2) * std::cos(langle * cavc::utils::pi<double>() / 180.0);
			}
			if (kLineList.size() == 0) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			//�ҵ���ͬ��yֵ
			else if (kLineList.count(_l.length) == 1) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			else {
				int state = 0;
				for (auto& key : kLineList) {
					//�ҵ�����������
					if (std::abs(key.first - std::fabs(_l.length)) <= 0.1) {
						kLineList[key.first].push_back(_l);
						state = 1;
						break;
					}
				}
				//û���ҵ�����������
				if (state == 0) {
					kLineList[std::fabs(_l.length)].push_back(_l);
				}
			}
		}
		int a = 0;
		cavc::Vector2<double> directionv;
		//ͨ����Ծ�������
		for (auto& kl : kLineList) {
			//�ж���һ�������ϵ��߶εķ����Ƿ���ǰһ�������ϵ��߶ι�ϵ����ͬ�ͷ��򣬲���ͬ��Ĭ��
			if (a != 0) {
				cavc::Vector2<double> kdirectionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);
				//���������ͬ�ͷ���
				if (cavc::fuzzyEqual(kdirectionv, directionv, 0.01)) {
					for (auto& al : kl.second) {
						auto tl = al.P1;
						al.P1 = al.P2;
						al.P2 = tl;
					}
				}
			}
			//��ʼͨ����������
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
			//�����������д�����ձ���
			for (auto& al : kl.second) {
				LastLine.push_back(al);
			}

			directionv = cavc::unitPerp(kl.second[0].P1 - kl.second[0].P2);

			a++;
		}
		count++;
	}
	//ˮƽ�ʹ�ֱ�߶�����������ϱ�һ��
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
			//���㵽������������
			cavc::Vector2<double> svector = basePoint - _l.P1;
			//���㵽�����յ������
			cavc::Vector2<double> evector = basePoint - _l.P2;
			//���㵽���ĽǶ�
			_l.sangle = std::atan2(svector.y(), svector.x()) * 180 / cavc::utils::pi<double>();
			//���㵽�յ�ĽǶ�
			_l.eangle = std::atan2(evector.y(), evector.x()) * 180 / cavc::utils::pi<double>();
			//���㵽���ľ���
			///ˮƽֱ��
			if (langle == 0) {
				_l.length = _l.P1.y() - basePoint.y();
			}
			///��ֱֱ��
			else if (langle == 90) {
				_l.length = _l.P1.x() - basePoint.x();
			}
			if (kLineList.size() == 0) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			//�ҵ���ͬ��yֵ
			else if (kLineList.count(_l.length) == 1) {
				kLineList[std::fabs(_l.length)].push_back(_l);
				continue;
			}
			else {
				int state = 0;
				for (auto& key : kLineList) {
					//�ҵ�����������
					if (std::abs(key.first - std::fabs(_l.length)) <= 0.1) {
						kLineList[key.first].push_back(_l);
						state = 1;
						break;
					}
				}
				//û���ҵ�����������
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
				//���������ͬ�ͷ���
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
				//�������ͬ�ͷ���
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
					//���x��ȱȽ�y
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
					//���y��ȱȽ�x
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
					//���x��ȱȽ�y
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
					//���y��ȱȽ�x
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
	
	//����
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
		//���Ϊ����
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
			//�Ƿ��ڰ�·���ϵ��߶�
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

	//������Ǳ߿�
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
		//Ĭ��
	case 0:
		LeftDown(BasePoint[0], 0, info.error, info.flycutspace,false, BasePoint, LineList, LastLine);
		break;
		//����
	case 1:
		LeftDown(BasePoint[0], 0, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//����
	case 2:
		LeftDown(BasePoint[1], 1, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//����
	case 3:
		LeftDown(BasePoint[3], 3, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
		//����
	case 4:
		LeftDown(BasePoint[2], 2, info.error, info.flycutspace, false, BasePoint, LineList, LastLine);
		break;
	default:
		break;
	}

	//�ѱ߿���ӵ����
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
	//����
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
	
	//��������Ϣ���������б���
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

		//��ȡ�����İ�·��
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

	//�Ѱ�·�����������
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
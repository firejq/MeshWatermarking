#pragma once
#ifndef _WATER_MARK_
#define _WATER_MARK_
#include"EigenDeform.h"

class WaterMark
{
public:
	WaterMark();
	~WaterMark();

	//void Init(PolygonMesh::Mesh * _mesh);

	//PolygonMesh::Mesh * my_embed_wm(PolygonMesh::Mesh * _mesh);
	//void my_extract_wm(PolygonMesh::Mesh * _mesh);

	void createA();//�������ԭʼˮӡ���У�����a /*ΪVecA��ֵ*/
	void createWB();//��������Ƕ��ˮӡ���У�����b' /*ΪVecB��ֵ������д�뵽�ļ�Wb.txt��*/
	void createP();//�����������P /*ΪP��ֵ����д�뵽�ļ�P.txt��*/

	void setM(int);
	//int getM();
	void setC(int);
	//int getC();
	void setAlpha(double);
	void setKey(int);

	//private:
	vector<int> vecA;//��
	vector<int> vecB;//��
	vector<int> P;//��

	int m;//ԭʼˮӡλ��
	int c;//��Ƭ����
	int key;//��
	double alpha;//��
	//int chip_rate;//��

};

#endif
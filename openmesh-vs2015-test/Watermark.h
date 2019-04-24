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

	void createA();//生成随机原始水印序列，数组a /*为VecA赋值*/
	void createWB();//生成最终嵌入水印序列，数组b' /*为VecB赋值，并将写入到文件Wb.txt中*/
	void createP();//生成随机序列P /*为P赋值，并写入到文件P.txt中*/

	void setM(int);
	//int getM();
	void setC(int);
	//int getC();
	void setAlpha(double);
	void setKey(int);

	//private:
	vector<int> vecA;//？
	vector<int> vecB;//？
	vector<int> P;//？

	int m;//原始水印位数
	int c;//码片速率
	int key;//？
	double alpha;//？
	//int chip_rate;//？

};

#endif
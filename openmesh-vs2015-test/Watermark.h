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

	//vector<int> all_labels;//存储每个顶点的分割编号
	//vector<bool> is_boundray;//标记每个顶点是否为边界，是为1
	//int num_of_cluster;//网格分割块数

	//map<int, vector<int>>Area_map;//用于记录分割情况
	//vector< vector<int> > readArea;//用于读取分割情况

	//							   //存储所有区域的非边界点个数，由于每个区域的水印位数应该相同
	//							   //因此要找出所有区域中，非边界点个数最少的区域
	//double *all_non_bdy;//存储每个区域的非边界顶点个数

	//					//记录分割结果


	//					/*计算曲率时要用到的参数*/
	//PolygonMesh::Mesh * myMesh;//传入的网格
	//vector<double> m_curvatureKmin;
	//vector<double> m_curvatureKmax;

	//int num_vtx;//网格顶点数
	//int cnt_non_bdy;//非边界顶点个数

	//vector<RowVectorX>     m_eigenVector;//特征向量--行向量表示
	//MatrixX3              R_matrix;
	//MatrixXX              E_matrix;

	//vector<int> map_vec;//区域中顶点的新旧编号对应关系
	//					//vector<int> newmap;//map_vec中，交换键值和值
};

#endif
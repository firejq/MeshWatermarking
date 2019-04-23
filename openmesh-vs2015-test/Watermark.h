#pragma once
#ifndef _WATER_MARK_
#define _WATER_MARK_
#include"EigenDeform.h"
#include<map>

using namespace GeometryProcess;
using namespace std;

class WaterMark
{
public:
	WaterMark();
	~WaterMark();

	void Init(PolygonMesh::Mesh * _mesh);

	PolygonMesh::Mesh * my_embed_wm(PolygonMesh::Mesh * _mesh);
	void my_extract_wm(PolygonMesh::Mesh * _mesh);

	//void segmentationMesh(PolygonMesh::Mesh * _mesh);
	//void visualSegmentations(PolygonMesh::Mesh * _mesh);

	//void find_Boundary(PolygonMesh::Mesh * _mesh, int area);
	//void cal_area_Lap(PolygonMesh::Mesh * _mesh, int area);
	//void cal_area_Rmatrix(PolygonMesh::Mesh * _mesh, int area);
	//void embed_by_area(PolygonMesh::Mesh * _mesh);
	//void extract_by_area(PolygonMesh::Mesh * _mesh);

	//void saveSegmentation();

	///*计算主曲率*/
	//void calCurvature();//使用该函数前要调用Init函数
	//					/*计算边确定的3*3矩阵*/
	//void calEdgeVec(OpenMesh::HalfedgeHandle he, Matrix33 & estiMat);
	///*计算半边所确定的两三角形面法向夹角及面积*/
	//void calAngleTris(OpenMesh::HalfedgeHandle he, double & _angle, double & _area);
	///*计算三角形面积*/
	//double calTriArea(OpenMesh::FaceHandle faceH);
	///*计算估计后切平面3×3矩阵的特征值*/
	//void calculateCurvature(const Matrix33 & tanMat, double & Kmin, double & Kmax);
	///*可视化点*/
	//void visualCoordiate(vector<double> & diffCoor);
	///*对曲率做一个归一化*/
	//void guassNormaize(vector<double> & _cur, double & mean_);

	void createA();//生成随机原始水印序列，数组a
	void createWB();//生成最终嵌入水印虚拟，数组b' 
	void createP();//生成随机序列P

	void setM(int);
	int getM();

	void setC(int);
	int getC();

	void setAlpha(double);
	void setKey(int);

	//private:
	vector<int> vecA;//？
	vector<int> vecB;//？
	vector<int> P;
	int m;//原始水印位数
	int c;//码片速率
	int key;
	double alpha;
	int chip_rate;

	vector<int> all_labels;//存储每个顶点的分割编号
	vector<bool> is_boundray;//标记每个顶点是否为边界，是为1
	int num_of_cluster;//网格分割块数

	map<int, vector<int>>Area_map;//用于记录分割情况
	vector< vector<int> > readArea;//用于读取分割情况

								   //存储所有区域的非边界点个数，由于每个区域的水印位数应该相同
								   //因此要找出所有区域中，非边界点个数最少的区域
	double *all_non_bdy;//存储每个区域的非边界顶点个数

						//记录分割结果


						/*计算曲率时要用到的参数*/
	PolygonMesh::Mesh * myMesh;//传入的网格
	vector<double> m_curvatureKmin;
	vector<double> m_curvatureKmax;

	int num_vtx;//网格顶点数
	int cnt_non_bdy;//非边界顶点个数

	vector<RowVectorX>     m_eigenVector;//特征向量--行向量表示
	MatrixX3              R_matrix;
	MatrixXX              E_matrix;

	vector<int> map_vec;//区域中顶点的新旧编号对应关系
						//vector<int> newmap;//map_vec中，交换键值和值
};

#endif
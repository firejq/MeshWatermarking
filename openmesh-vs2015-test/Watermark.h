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

	///*����������*/
	//void calCurvature();//ʹ�øú���ǰҪ����Init����
	//					/*�����ȷ����3*3����*/
	//void calEdgeVec(OpenMesh::HalfedgeHandle he, Matrix33 & estiMat);
	///*��������ȷ�������������淨��нǼ����*/
	//void calAngleTris(OpenMesh::HalfedgeHandle he, double & _angle, double & _area);
	///*�������������*/
	//double calTriArea(OpenMesh::FaceHandle faceH);
	///*������ƺ���ƽ��3��3���������ֵ*/
	//void calculateCurvature(const Matrix33 & tanMat, double & Kmin, double & Kmax);
	///*���ӻ���*/
	//void visualCoordiate(vector<double> & diffCoor);
	///*��������һ����һ��*/
	//void guassNormaize(vector<double> & _cur, double & mean_);

	void createA();//�������ԭʼˮӡ���У�����a
	void createWB();//��������Ƕ��ˮӡ���⣬����b' 
	void createP();//�����������P

	void setM(int);
	int getM();

	void setC(int);
	int getC();

	void setAlpha(double);
	void setKey(int);

	//private:
	vector<int> vecA;//��
	vector<int> vecB;//��
	vector<int> P;
	int m;//ԭʼˮӡλ��
	int c;//��Ƭ����
	int key;
	double alpha;
	int chip_rate;

	vector<int> all_labels;//�洢ÿ������ķָ���
	vector<bool> is_boundray;//���ÿ�������Ƿ�Ϊ�߽磬��Ϊ1
	int num_of_cluster;//����ָ����

	map<int, vector<int>>Area_map;//���ڼ�¼�ָ����
	vector< vector<int> > readArea;//���ڶ�ȡ�ָ����

								   //�洢��������ķǱ߽�����������ÿ�������ˮӡλ��Ӧ����ͬ
								   //���Ҫ�ҳ����������У��Ǳ߽��������ٵ�����
	double *all_non_bdy;//�洢ÿ������ķǱ߽綥�����

						//��¼�ָ���


						/*��������ʱҪ�õ��Ĳ���*/
	PolygonMesh::Mesh * myMesh;//���������
	vector<double> m_curvatureKmin;
	vector<double> m_curvatureKmax;

	int num_vtx;//���񶥵���
	int cnt_non_bdy;//�Ǳ߽綥�����

	vector<RowVectorX>     m_eigenVector;//��������--��������ʾ
	MatrixX3              R_matrix;
	MatrixXX              E_matrix;

	vector<int> map_vec;//�����ж�����¾ɱ�Ŷ�Ӧ��ϵ
						//vector<int> newmap;//map_vec�У�������ֵ��ֵ
};

#endif
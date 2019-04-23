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

	//vector<int> all_labels;//�洢ÿ������ķָ���
	//vector<bool> is_boundray;//���ÿ�������Ƿ�Ϊ�߽磬��Ϊ1
	//int num_of_cluster;//����ָ����

	//map<int, vector<int>>Area_map;//���ڼ�¼�ָ����
	//vector< vector<int> > readArea;//���ڶ�ȡ�ָ����

	//							   //�洢��������ķǱ߽�����������ÿ�������ˮӡλ��Ӧ����ͬ
	//							   //���Ҫ�ҳ����������У��Ǳ߽��������ٵ�����
	//double *all_non_bdy;//�洢ÿ������ķǱ߽綥�����

	//					//��¼�ָ���


	//					/*��������ʱҪ�õ��Ĳ���*/
	//PolygonMesh::Mesh * myMesh;//���������
	//vector<double> m_curvatureKmin;
	//vector<double> m_curvatureKmax;

	//int num_vtx;//���񶥵���
	//int cnt_non_bdy;//�Ǳ߽綥�����

	//vector<RowVectorX>     m_eigenVector;//��������--��������ʾ
	//MatrixX3              R_matrix;
	//MatrixXX              E_matrix;

	//vector<int> map_vec;//�����ж�����¾ɱ�Ŷ�Ӧ��ϵ
	//					//vector<int> newmap;//map_vec�У�������ֵ��ֵ
};

#endif
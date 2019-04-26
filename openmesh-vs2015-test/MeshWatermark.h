#pragma once
#ifndef MESHWATERMARK_H_
#define MESHWATERMARK_H_

#include "Init.h"

class MeshWatermark
{
public:
	MeshWatermark();
	~MeshWatermark() { }
	
	void Init(PolygonMesh::Mesh * _mesh, string watermark_path);
	void Init(PolygonMesh::Mesh * _mesh, string watermark_path, int c, double a);
	/*嵌入水印*/
	bool embedWatermark();
	/*提取水印*/
	bool extractWatermark(string extr_watermark);

private:
	/*计算拉普拉斯矩阵L = D - A，并进行特征值分解*/
	void calLap_Matrix();
	/*计算频谱系数矩阵E . R = V,
	E_matrix为单位化的特征向量矩阵，V_matrix是顶点坐标矩阵*/
	bool calR_Matrix();
	/*将得出的特征向量单位化*/
	void normVec();

	PolygonMesh::Mesh * m_oriMesh;//输入的原始网格
								  //PolygonMesh::Mesh * m_wmMesh;//输入的水印网格

	MatrixXX            m_verPostion;//顶点坐标
	int                 m_vertexNum;//网格顶点个数-用于resize矩阵
	int chip_rate;//@firejq：码片速率
	double alpha;

	vector<RowVectorX>     m_eigenVector;//特征向量--行向量表示
	vector<RowVectorX>     m_normVector;//单位化的特征向量
	vector<double> m_eigenValue;//特征值
	MatrixX3              R_matrix;//频谱系数矩阵，nx3
	MatrixX3              V_matrix;//顶点坐标矩阵，nx3
	MatrixXX			    E_matrix;//单位化特征向量矩阵，nxn
	MatrixXX				Lap_matrix;//拉普拉斯矩阵，L = D - A
									   //MatrixXX              m_eigenValue;//nxn 对角矩阵

	MatrixX3              wR_matrix;//频谱系数矩阵，nx3
	MatrixX3              wV_matrix;//顶点坐标矩阵，nx3
	MatrixXX			    wE_matrix;//单位化特征向量矩阵，nxn
	MatrixXX				wLap_matrix;//拉普拉斯矩阵，L = D - A
	Engine            * m_engine;//matlab引擎

	string				wm_path;//水印图像的路径


	/*内部类WatermarkSeq：生成水印序列的相关方法和属性*/
	class WatermarkSeq {
	public:
		void initWB(string watermark_path);//初始化水印图像的信息
		void createWB(int m_vertexNum);//生成最终嵌入水印序列，数组b' /*为VecB赋值，并将写入到文件Wb.txt中*/
		void createP();//生成随机序列P /*为P赋值，并写入到文件P.txt中*/

		void setM(int);
		//int getM();
		void setC(int);
		//int getC();
		void setAlpha(double);
		void setKey(int);

		//private:
		vector<int> vecB;//水印序列，长度为传入的二值图像的长度
		vector<int> P;//伪随机序列(PRNS)，长度为m*c（即网格顶点个数)

		int m;//原始水印位数
		int c;//码片速率
		int key;//生成PRNS的密钥，暂时用不到
		double alpha;//调制幅度(modulation amplitude)
	
		string wmImgPath;//输入的图片水印的位置
		int wmSeqLength;//输入的图片水印序列的总长度
		int wmImgWidth;//输入的图片水印的宽
		int wmImgLength;//输入的图片水印的长

	
	};
};

#endif
#include "EigenDeform.h"
#include "WaterMark.h"
#include <fstream>
#include <math.h>
#include <string>

using namespace GeometryProcess::PolygonMesh;
//using namespace GeometryProcess::PolyhedralMesh;
using namespace std;

//This codes were writted for "Eigen Deformation of 3D Models"
namespace GeometryProcess
{

	EigenDeformation::EigenDeformation()
	{
		m_oriMesh = NULL;
		m_engine = NULL;//初始化matlab引擎
		m_bdyCenter.setZero();
		m_meanDis = 0.0;//copy from yuan for segementation
	}
	//add by shiqun
	void EigenDeformation::Init(PolygonMesh::Mesh * _mesh)
	{
		chip_rate = 7;
		m_oriMesh = _mesh;
		//m_wmMesh =_mesh;
		m_vertexNum = m_oriMesh->n_vertices();
	}
	//

	//void EigenDeformation::Init(BaseEntity * _meshEntity)
	//{
	//	//	m_meshEntity = _meshEntity;
	//	m_oriMesh = ((PolygonMesh::PGMeshEntity *)_meshEntity)->get_mesh();
	//	m_vertexNum = m_oriMesh->n_vertices();

	//	int vex_index = 0;

	//	for (auto v_it = m_oriMesh->vertices_begin(); v_it != m_oriMesh->vertices_end(); ++v_it)//初始化顶点坐标
	//	{
	//		auto  point = m_oriMesh->point(v_it.handle());//原始坐标①
	//		m_verPostion(vex_index, 0) = point[0];
	//		m_verPostion(vex_index, 1) = point[1];
	//		m_verPostion(vex_index, 2) = point[2];
	//		vex_index++;
	//	}
	//	//m_verPos_iter = m_verPostion;//赋初始顶点值  
	//}

	//copy from yuan 使用面积计算拉普拉斯矩阵
	bool EigenDeformation::calLapMatrix()
	{
		const int nEdges = m_oriMesh->n_edges();
		const int nFaces = m_oriMesh->n_faces();
		const int totalElements = m_vertexNum + nEdges * 2;//laplacian矩阵的对角线上元素为顶点的度，每行除了对角线元素外的

		double * xCoordinate = new double[m_vertexNum];
		double * yCoordinate = new double[m_vertexNum];
		double * zCoordinate = new double[m_vertexNum];

		double * px = xCoordinate;
		double * py = yCoordinate;
		double * pz = zCoordinate;
		//获取网格上的点的坐标
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, ++px, ++py, ++pz)
		{
			auto t_p = m_oriMesh->point(*vit);
			px[0] = t_p[0];
			py[0] = t_p[1];
			pz[0] = t_p[2];
		}
		double* tTri = new double[nFaces * 3];
		int fcounter = 0;
		for (auto fit = m_oriMesh->faces_begin(); fit != m_oriMesh->faces_end(); ++fit)
		{
			for (auto fvit = m_oriMesh->fv_iter(*fit); fvit != m_oriMesh->fv_end(*fit); ++fvit)
			{
				tTri[fcounter] = (*fvit).idx() + 1;
				++fcounter;
			}
		}
		double* Tri = new double[nFaces * 3];
		for (int i = 0; i < nFaces * 3; ++i)
		{
			Tri[i] = tTri[(i%nFaces) * 3 + (i / nFaces)];//变为按列的数组--没有看明白
		}
		if ((!m_engine && !(m_engine = engOpen(NULL))))// 使用matlab算laplacian矩阵及特征值
		{
			return false;
		}
		engSetVisible(m_engine, 1);
		mxArray *XCoordinate = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);//把网格数据传入matlab
		mxArray *YCoordinate = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);//创建m_vertexNum行，1列的实数组
		mxArray *ZCoordinate = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);
		mxArray *TRIV = mxCreateDoubleMatrix(nFaces, 3, mxREAL);
		memcpy((void *)mxGetPr(XCoordinate), (void *)xCoordinate, m_vertexNum * sizeof(double));
		memcpy((void *)mxGetPr(YCoordinate), (void *)yCoordinate, m_vertexNum * sizeof(double));
		memcpy((void *)mxGetPr(ZCoordinate), (void *)zCoordinate, m_vertexNum * sizeof(double));
		memcpy((void *)mxGetPr(TRIV), (void *)Tri, (nFaces * 3) * sizeof(double));
		engPutVariable(m_engine, "X", XCoordinate);
		engPutVariable(m_engine, "Y", YCoordinate);
		engPutVariable(m_engine, "Z", ZCoordinate);
		engPutVariable(m_engine, "TRIV", TRIV);

		engEvalString(m_engine, "cd('D:\\code\\GeometryProcessing-1\\mat_code')");//shenzhen
		engEvalString(m_engine, "shape = struct('X',{X},'Y',{Y},'Z',{Z},'TRIV',{TRIV});");//add by shiqun,struct 是matlab自带的函数


																						  //buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		//modified by shiqun
		/*
		原函数是给matlab一个形状，然后计算前300个特征值和特征向量，computeLacian函数的第二个
		输入参数只能是数字，不能用m_vertexNum。下面的函数调用的是matlab的代码，其输入的参数
		是一个数值，而不是一个变量
		*/
		engEvalString(m_engine, "[eigVector,eigValue,W,D] = computeLaplacian(shape,233);");

		printf("%s", buffer);//当matlab代码出错时，用来输出调试信息


		mxArray * eigenVectorObj = NULL;
		eigenVectorObj = engGetVariable(m_engine, "eigVector");
		mxArray * eigenValueObj = NULL;
		eigenValueObj = engGetVariable(m_engine, "eigValue");

		/*add by shiqun
		eigenValue矩阵是一个一维矩阵，矩阵元素个数为1180，每个元素是对应的特征值
		*/
		double  * eigenVector = NULL;
		eigenVector = (double*)mxGetData(eigenVectorObj);
		double  * eigenValue = NULL;
		eigenValue = (double*)mxGetData(eigenValueObj);
		//		
		RowVectorX eVector_iter;
		eVector_iter.resize(m_vertexNum);

		//add by shiqun for test
		//std::cout<<sizeof(eigenVector)<<std::endl;

		for (int eig_index = 0; eig_index < m_vertexNum; eig_index++)//把eigenvector赋给m_eigenVector
		{
			eVector_iter.setZero();
			for (int ver_pos = 0; ver_pos < m_vertexNum; ver_pos++)
			{
				eVector_iter[ver_pos] = eigenVector[eig_index * m_vertexNum + ver_pos];
			}
			m_eigenVector.push_back(eVector_iter);
		}

		mxDestroyArray(XCoordinate);
		mxDestroyArray(YCoordinate);
		mxDestroyArray(ZCoordinate);
		mxDestroyArray(TRIV);
		mxDestroyArray(eigenValueObj);
		mxDestroyArray(eigenVectorObj);

		if (xCoordinate)
		{
			delete[] xCoordinate;
			xCoordinate = NULL;
		}
		if (yCoordinate)
		{
			delete[] yCoordinate;
			yCoordinate = NULL;
		}
		if (zCoordinate)
		{
			delete[] zCoordinate;
			zCoordinate = NULL;
		}
		if (tTri)
		{
			delete[] tTri;
			tTri = NULL;
		}
		if (Tri)
		{
			delete[] Tri;
			Tri = NULL;
		}
		/*if (m_engine)
		{
		engClose(m_engine);
		m_engine = NULL;
		}*/
		return true;
	}

	//L = D - A
	void EigenDeformation::calLap_Matrix()
	{
		MatrixXX D_matrix;//nxn度矩阵
		MatrixXX A_matrix;//nxn邻接矩阵
		D_matrix.setZero(m_vertexNum, m_vertexNum);
		A_matrix.setZero(m_vertexNum, m_vertexNum);
		Lap_matrix.setZero(m_vertexNum, m_vertexNum);//nxn拉普拉斯矩阵

													 //计算度矩阵
		int v_index = 0;
		for (auto v_it = m_oriMesh->vertices_begin(); v_it != m_oriMesh->vertices_end(); ++v_it, v_index++)
		{
			int vnum = m_oriMesh->valence(*v_it);//获得顶点i的度
			D_matrix(v_index, v_index) = vnum;
		}

		//计算邻接矩阵
		for (auto v_it = m_oriMesh->vertices_begin(); v_it != m_oriMesh->vertices_end(); ++v_it)
		{
			OpenMesh::VertexHandle vertex1 = *v_it;
			int v1_idx = vertex1.idx();
			for (auto vv_it = m_oriMesh->vv_begin(*v_it); vv_it != m_oriMesh->vv_end(*v_it); ++vv_it)
			{
				OpenMesh::VertexHandle vertex2 = *vv_it;
				int v2_idx = vertex2.idx();
				A_matrix(v1_idx, v2_idx) = 1;
			}
		}

		//计算拉普拉斯矩阵
		Lap_matrix = D_matrix - A_matrix;

		// 使用matlab算laplacian矩阵及特征值
		int l_index = 0;
		double *PL = new double[m_vertexNum*m_vertexNum];//@firejq:???
		for (int l_col = 0; l_col < m_vertexNum; l_col++)
		{
			for (int l_row = 0; l_row < m_vertexNum; l_row++)
			{
				PL[l_index++] = Lap_matrix(l_row, l_col);
			}
		}

		if ((!m_engine && !(m_engine = engOpen(NULL))))// 使用matlab算laplacian矩阵及特征值
		{
			std::cout << "Matlab engine inits failed" << std::endl;
			return;
		}
		engSetVisible(m_engine, 1);

		mxArray *pL = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);
		memcpy((void *)mxGetPr(pL), (void *)PL, m_vertexNum*m_vertexNum * sizeof(double));
		engPutVariable(m_engine, "L", pL);

		engEvalString(m_engine, "cd('D:\\firejq\\毕设\\Watermarking\\mat_code')");//调用matalab代码

																				//buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "[eigVector,eigValue] = calcLaplacian(L);");
		//engEvalString(m_engine,"[eigVector,eigValue] = eigenDeform(L);");//只计算前6个最大的特征值

		printf("%s", buffer);//当matlab代码出错时，用来输出调试信息

		mxArray * eigenVectorObj = NULL;
		eigenVectorObj = engGetVariable(m_engine, "eigVector");
		mxArray * eigenValueObj = NULL;
		eigenValueObj = engGetVariable(m_engine, "eigValue");

		/*add by shiqun
		eigenValue矩阵是一个一维矩阵，矩阵元素个数为1180，每个元素是对应的特征值
		*/
		double  * eigenVector = NULL;
		eigenVector = (double*)mxGetData(eigenVectorObj);
		double  * eigenValue = NULL;
		eigenValue = (double*)mxGetData(eigenValueObj);
		//		
		RowVectorX eVector_iter;
		eVector_iter.resize(m_vertexNum);

		//获取特征向量
		for (int eig_index = 0; eig_index < m_vertexNum; eig_index++)//把eigenvector赋给m_eigenVector
		{
			eVector_iter.setZero();
			for (int ver_pos = 0; ver_pos < m_vertexNum; ver_pos++)
			{
				eVector_iter[ver_pos] = eigenVector[eig_index * m_vertexNum + ver_pos];
			}
			m_eigenVector.push_back(eVector_iter);
		}

		//获取特征值
		for (int val_index = 0; val_index < m_vertexNum; val_index++)
		{
			m_eigenValue.push_back(eigenValue[val_index]);
		}
	}
	//add by shiqun
	//将特征向量单位化

	void::EigenDeformation::normVec()
	{
		m_normVector = m_eigenVector;//将特征向量复制给将要单位化的向量
		VectorXd  normV(m_vertexNum);//声明一个Eigen中的vector类型的临时变量，用来存储一个单独的特征向量
									 //每次将一个特征向量复制给中间变量normV，然后normV调用单位化函数，再将normV赋给m_normVector
		for (int index = 0; index < m_vertexNum; index++)//这里把条件写成了index<1，这种粗心造成的错误真是够够的了。。
		{
			normV = m_normVector[index];
			normV.normalize();
			m_normVector.push_back(normV);
		}
	}

	//add by shiqun
	//将顶点坐标投影到特征向量上,修改频谱系数，并将频谱系数反变换到原始坐标中
	//计算频谱系数
	bool EigenDeformation::calR_Matrix()
	{
		//顶点坐标
		VectorXd xCoordinate(m_vertexNum);
		VectorXd yCoordinate(m_vertexNum);
		VectorXd zCoordinate(m_vertexNum);

		//将顶点坐标赋值给VectorXd
		int j = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, j++)
		{
			auto t_p = m_oriMesh->point(*vit);

			xCoordinate[j] = t_p[0];
			yCoordinate[j] = t_p[1];
			zCoordinate[j] = t_p[2];
		}

		//计算顶点坐标矩阵V_matrix
		V_matrix.setZero(m_vertexNum, 3);
		V_matrix.col(0) = xCoordinate;//V_matrix是顶点坐标矩阵，nx3
		V_matrix.col(1) = yCoordinate;
		V_matrix.col(2) = zCoordinate;

		int vetex_index = 0;
		double *vetex_coor = new double[m_vertexNum * 3];
		for (int col_index = 0; col_index < 3; col_index++)
		{
			for (int row_index = 0; row_index < m_vertexNum; row_index++)
			{
				vetex_coor[vetex_index++] = V_matrix(row_index, col_index);
			}
		}

		//计算特征向量矩阵E_matrix
		E_matrix.setZero(m_vertexNum, m_vertexNum);
		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix是特征向量矩阵，nxn
		}

		int eig_index = 0;
		double *eigVec = new double[m_vertexNum*m_vertexNum];
		for (int axi_col = 0; axi_col < m_vertexNum; axi_col++)
		{
			for (int axi_row = 0; axi_row < m_vertexNum; axi_row++)
			{
				eigVec[eig_index++] = E_matrix(axi_row, axi_col);
			}
		}

		//将单位化的特征向量写入文件中
		ofstream NormEfile;
		NormEfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Ne.txt", ios_base::out);
		if (NormEfile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				for (int j = 0; j < m_vertexNum; j++)
				{
					NormEfile << E_matrix(i, j) << " ";
				}
			}
		}
		NormEfile.close();

		/*调用matlab计算频谱系数矩阵R_matrix*/
		mxArray *Vetex_Coor = mxCreateDoubleMatrix(m_vertexNum, 3, mxREAL);//把网格顶点坐标矩阵传入matlab
		mxArray *Eigen_vetor = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);

		memcpy((void *)mxGetPr(Eigen_vetor), (void *)eigVec, (m_vertexNum * m_vertexNum) * sizeof(double));
		memcpy((void *)mxGetPr(Vetex_Coor), (void *)vetex_coor, (m_vertexNum * 3) * sizeof(double));

		engPutVariable(m_engine, "A", Eigen_vetor);//矩阵方程为AX = B,其中A = E_matrix， X = R_matrix， B = V_matrix
		engPutVariable(m_engine, "B", Vetex_Coor);

		//buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "cd('D:\\firejq\\毕设\\Watermarking\\mat_code')");//shenzhen
		engEvalString(m_engine, "R = calcR(A, B);");

		printf("%s", buffer);//当matlab代码出错时，用来输出调试信息

		mxArray * r_matrix = NULL;
		r_matrix = engGetVariable(m_engine, "R");
		double  * Rmatrix = NULL;
		Rmatrix = (double*)mxGetData(r_matrix);

		R_matrix.setZero(m_vertexNum, 3);
		int r_index = 0;
		for (int r_col = 0; r_col < 3; r_col++)
		{
			for (int r_row = 0; r_row < m_vertexNum; r_row++)
			{
				R_matrix(r_row, r_col) = Rmatrix[r_index++];
			}
		}

		mxDestroyArray(Vetex_Coor);
		mxDestroyArray(Eigen_vetor);
		mxDestroyArray(r_matrix);
		// 测试时不用关闭
		/*if (m_engine)
		{
		engClose(m_engine);
		m_engine = NULL;
		}*/

		return true;
	}

	bool EigenDeformation::embedWatermark()
	{
		//频谱系数坐标
		VectorXd Rs(m_vertexNum);
		VectorXd Rt(m_vertexNum);
		VectorXd Ru(m_vertexNum);
		Rs = R_matrix.col(0);
		Rt = R_matrix.col(1);
		Ru = R_matrix.col(2);

		//修改后的频谱系数坐标
		VectorXd Rs_hat(m_vertexNum);
		VectorXd Rt_hat(m_vertexNum);
		VectorXd Ru_hat(m_vertexNum);

		//反变换后的顶点坐标
		VectorXd newxCoordinate(m_vertexNum);
		VectorXd newyCoordinate(m_vertexNum);
		VectorXd newzCoordinate(m_vertexNum);
		newxCoordinate.setZero();
		newyCoordinate.setZero();
		newzCoordinate.setZero();

		WaterMark mytest;
		mytest.setM((int)ceil((double)m_vertexNum / chip_rate));//设置原始水印位数
		mytest.setC(chip_rate);//设置码片速率
		mytest.setAlpha(0.005);//？
		mytest.setKey(7);//？

		mytest.createA();
		mytest.createWB();
		mytest.createP();


		for (int i_for_get = 0; i_for_get < m_vertexNum; i_for_get++)
		{
			Rs_hat[i_for_get] = Rs[i_for_get] + mytest.vecB[i_for_get] * mytest.P[i_for_get] * mytest.alpha;//Rs[i]' = Rs[i] + bi' . pi . a 
			Rt_hat[i_for_get] = Rt[i_for_get] + mytest.vecB[i_for_get] * mytest.P[i_for_get] * mytest.alpha;
			Ru_hat[i_for_get] = Ru[i_for_get] + mytest.vecB[i_for_get] * mytest.P[i_for_get] * mytest.alpha;
		}

		/**
		将原始网格的频谱系数R矩阵写入文件中，
		这样提取水印的时候就不必再读入原始网格，只需计算出读入的水印网格中的wR矩阵，
		然后从文件中读取出R矩阵，两者做差值，再进行判断步骤即可。
		**/
		//Rs
		ofstream Rsfile;
		Rsfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Rs.txt", ios_base::out);
		if (Rsfile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				Rsfile << Rs[i] << " ";
			}
		}
		Rsfile.close();
		//Rt
		ofstream Rtfile;
		Rtfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Rt.txt", ios_base::out);
		if (Rtfile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				Rtfile << Rt[i] << " ";
			}
		}
		Rtfile.close();
		//Ru
		ofstream Rufile;
		Rufile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Ru.txt", ios_base::out);
		if (Rufile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				Rufile << Ru[i] << " ";
			}
		}
		Rufile.close();

		VectorXd  normV(m_vertexNum);//Eigen类型的数组，用作临时变量
		for (int i_for_set = 0; i_for_set< m_vertexNum; i_for_set++)
		{
			normV = m_normVector[i_for_set];
			newxCoordinate = newxCoordinate + Rs_hat[i_for_set] * normV;//X = ∑Rs[i] . ei
			newyCoordinate = newyCoordinate + Rt_hat[i_for_set] * normV;
			newzCoordinate = newzCoordinate + Ru_hat[i_for_set] * normV;
		}


		//将修改后的顶点坐标赋回给网格上的顶点
		int vi = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, vi++)
		{
			auto& t_p = m_oriMesh->point(*vit);

			t_p[0] = newxCoordinate[vi];
			t_p[1] = newyCoordinate[vi];
			t_p[2] = newzCoordinate[vi];
		}

		m_oriMesh->update_face_normals();//Update normal vectors for all faces

		return true;
	}

	//提取水印
	//在嵌入之前将水印信息写入文件中，
	//提取水印时，将提取出的数据也写入文件中，然后比较两个文件中的信息是否一致
	bool EigenDeformation::extractWatermark()
	{
		//原始网格与水印网格都做特征值分解，然后特征值相减

		//顶点坐标
		VectorXd wxCoordinate(m_vertexNum);
		VectorXd wyCoordinate(m_vertexNum);
		VectorXd wzCoordinate(m_vertexNum);

		//将顶点坐标赋值给VectorXd
		int k = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, k++)
		{
			auto t_p = m_oriMesh->point(*vit);

			wxCoordinate[k] = t_p[0];
			wyCoordinate[k] = t_p[1];
			wzCoordinate[k] = t_p[2];
		}

		V_matrix.setZero(m_vertexNum, 3);

		V_matrix.col(0) = wxCoordinate;//V_matrix是顶点坐标矩阵，nx3
		V_matrix.col(1) = wyCoordinate;
		V_matrix.col(2) = wzCoordinate;

		int vetex_index = 0;
		double *vetex_coor = new double[m_vertexNum * 3];
		for (int col_index = 0; col_index < 3; col_index++)
		{
			for (int row_index = 0; row_index < m_vertexNum; row_index++)
			{
				vetex_coor[vetex_index++] = V_matrix(row_index, col_index);
			}
		}

		E_matrix.setZero(m_vertexNum, m_vertexNum);
		//调用matlab计算特征向量
		//for ( int e_index = 0; e_index< m_vertexNum; e_index++)
		//{ 
		//	E_matrix.col(e_index) = m_normVector[e_index];//E_matrix是特征向量矩阵，nxn
		//}

		//int eig_index = 0;
		//double *eigVec = new double[m_vertexNum*m_vertexNum];
		//for (int axi_col = 0; axi_col < m_vertexNum; axi_col++ )
		//{
		//	for ( int axi_row =0; axi_row < m_vertexNum; axi_row++)
		//	{	
		//		eigVec[eig_index++] = E_matrix(axi_row,axi_col);
		//	}
		//}

		////从文件中将原始网格的E矩阵读取出来
		double *pE = new double[m_vertexNum * m_vertexNum];
		ifstream NERfile;
		NERfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Ne.txt", ios_base::in);
		if (NERfile)
		{
			// read into memory
			NERfile.seekg(0, NERfile.end);
			int length = NERfile.tellg();
			NERfile.seekg(0, NERfile.beg);

			char *buffer = new char[length];
			NERfile.read(buffer, length);
			NERfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (iss >> pE[i++]);
			delete[] buffer;
			// print or use it.
		}
		NERfile.close();

		int ecount = 0;
		for (int i = 0; i < m_vertexNum; i++)
		{
			for (int j = 0; j < m_vertexNum; j++)
			{
				E_matrix(i, j) = pE[ecount++];
			}
		}

		int eig_index = 0;
		double *eigVec = new double[m_vertexNum*m_vertexNum];
		for (int axi_col = 0; axi_col < m_vertexNum; axi_col++)
		{
			for (int axi_row = 0; axi_row < m_vertexNum; axi_row++)
			{
				eigVec[eig_index++] = E_matrix(axi_row, axi_col);
			}
		}

		//调用matlab
		if ((!m_engine && !(m_engine = engOpen(NULL))))
		{
			return false;
		}
		engSetVisible(m_engine, 1);

		mxArray *Vetex_Coor = mxCreateDoubleMatrix(m_vertexNum, 3, mxREAL);//把网格顶点坐标矩阵传入matlab
		mxArray *Eigen_vetor = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);

		memcpy((void *)mxGetPr(Eigen_vetor), (void *)eigVec, (m_vertexNum * m_vertexNum) * sizeof(double));
		memcpy((void *)mxGetPr(Vetex_Coor), (void *)vetex_coor, (m_vertexNum * 3) * sizeof(double));

		engPutVariable(m_engine, "A", Eigen_vetor);//矩阵方程为AX = B,其中A = E_matrix， X = R_matrix， B = V_matrix
		engPutVariable(m_engine, "B", Vetex_Coor);

		//buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "cd('D:\\firejq\\毕设\\Watermarking\\mat_code')");//shenzhen
		engEvalString(m_engine, "R = calcR(A, B);");

		printf("%s", buffer);//当matlab代码出错时，用来输出调试信息

		mxArray * r_matrix = NULL;
		r_matrix = engGetVariable(m_engine, "R");
		double  * Rmatrix = NULL;
		Rmatrix = (double*)mxGetData(r_matrix);

		R_matrix.setZero(m_vertexNum, 3);
		int r_index = 0;
		for (int r_col = 0; r_col < 3; r_col++)
		{
			for (int r_row = 0; r_row < m_vertexNum; r_row++)
			{
				R_matrix(r_row, r_col) = Rmatrix[r_index++];
			}
		}

		mxDestroyArray(Vetex_Coor);
		mxDestroyArray(Eigen_vetor);
		mxDestroyArray(r_matrix);
		// 测试时不用关闭
		if (m_engine)
		{
			engClose(m_engine);
			m_engine = NULL;
		}

		//频谱系数坐标
		VectorXd wRs(m_vertexNum);
		VectorXd wRt(m_vertexNum);
		VectorXd wRu(m_vertexNum);

		wRs = R_matrix.col(0);
		wRt = R_matrix.col(1);
		wRu = R_matrix.col(2);

		//将原始网格的频谱系数从文件中读取出来
		double *pRs = new double[m_vertexNum];
		double *pRt = new double[m_vertexNum];
		double *pRu = new double[m_vertexNum];

		//Rs
		ifstream Rsfile;
		Rsfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Rs.txt", ios_base::in);
		if (Rsfile)
		{
			// read into memory
			Rsfile.seekg(0, Rsfile.end);
			int length = Rsfile.tellg();
			Rsfile.seekg(0, Rsfile.beg);

			char *buffer = new char[length];
			Rsfile.read(buffer, length);
			Rsfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (iss >> pRs[i++]);
			delete[] buffer;
			// print or use it.
		}
		Rsfile.close();

		//Rt
		ifstream Rtfile;
		Rtfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Rt.txt", ios_base::in);
		if (Rtfile)
		{
			// read into memory
			Rtfile.seekg(0, Rtfile.end);
			int length = Rtfile.tellg();
			Rtfile.seekg(0, Rtfile.beg);

			char *buffer = new char[length];
			Rtfile.read(buffer, length);
			Rtfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (iss >> pRt[i++]);
			delete[] buffer;
			// print or use it.
		}
		Rtfile.close();

		//Ru
		ifstream Rufile;
		Rufile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Ru.txt", ios_base::in);
		if (Rufile)
		{
			// read into memory
			Rufile.seekg(0, Rufile.end);
			int length = Rufile.tellg();
			Rufile.seekg(0, Rufile.beg);

			char *buffer = new char[length];
			Rufile.read(buffer, length);
			Rufile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (iss >> pRu[i++]);
			delete[] buffer;
			// print or use it.
		}
		Rufile.close();

		//将P读取出来
		double *P = new double[m_vertexNum];
		ifstream Pfile;
		Pfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\P.txt", ios_base::in);
		if (Pfile)
		{
			// read into memory
			Pfile.seekg(0, Pfile.end);
			int length = Pfile.tellg();
			Pfile.seekg(0, Pfile.beg);

			char *buffer = new char[length];
			Pfile.read(buffer, length);
			Pfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (i < m_vertexNum)
			{
				iss >> P[i++];
			}
			delete[] buffer;
			// print or use it.
		}
		Pfile.close();

		double *Qx = new double[m_vertexNum];
		double *Qy = new double[m_vertexNum];
		double *Qz = new double[m_vertexNum];
		double *Q = new double[m_vertexNum];

		for (int qi = 0; qi < m_vertexNum; qi++)
		{
			Qx[qi] = chip_rate * (wRs[qi] - pRs[qi]) * P[qi];
			Qy[qi] = chip_rate * (wRt[qi] - pRt[qi]) * P[qi];
			Qz[qi] = chip_rate * (wRu[qi] - pRu[qi]) * P[qi];
			Q[qi] = 1 / 3.0 * (Qx[qi] + Qy[qi] + Qz[qi]);
		}


		double * wB = new double[m_vertexNum];
		for (int i = 0; i < m_vertexNum; i++)
		{
			wB[i] = _copysign(1.0, Q[i]);
		}

		//@firejq: 将提取出的水印序列写到文件，方便调试时查看
		ofstream Wb_extract;
		Wb_extract.open("Wb_extract.txt", ios_base::out);
		if (Wb_extract)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				Wb_extract << wB[i] << " ";
			}
		}
		Wb_extract.close();



		//将提原始的水印序列读出
		double *RB = new double[m_vertexNum];
		ifstream Wbfile;
		Wbfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Wb.txt", ios_base::in);
		if (Wbfile)
		{
			// read into memory
			Wbfile.seekg(0, Wbfile.end);
			int length = Wbfile.tellg();
			Wbfile.seekg(0, Wbfile.beg);

			char *buffer = new char[length];
			Wbfile.read(buffer, length);
			Wbfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (i<m_vertexNum)
			{
				iss >> RB[i++];
			}
			delete[] buffer;
			// print or use it.
		}
		Wbfile.close();

		//比较读取出的Wb和提取出的wB
		double Bcounter = 0.0;
		for (int bi = 0; bi < m_vertexNum; bi++)
		{
			if (RB[bi] == wB[bi])
			{
				Bcounter++;
			}
		}

		double corr = Bcounter / (double)m_vertexNum;
		cout << "The result is: " << corr << endl;


		//将提取出的wB存入文件
		//先转换为二值序列
		for (int k = 0; k < 1186; k++)
		{
			if (wB[k] == -1)
			{
				wB[k] = 0;
			}
		}
		ofstream EWbfile;
		EWbfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\AMesh\\Extract\\FlowerWb.txt", ios_base::out);
		int i = 0;
		if (EWbfile)
		{
			for (int r = 0; r < 33; r++)
			{
				for (int c = 0; c < 36; c++)
				{
					if (i < 1186)
					{
						EWbfile << wB[i++] << " ";
					}
					else
					{
						EWbfile << "1" << " ";
					}
				}
				EWbfile << endl;
			}
		}
		EWbfile.close();


		delete[] pRs;
		delete[]pRt;
		delete[]pRu;
		delete[]P;
		delete[]Qx;
		delete[]Qy;
		delete[]Qz;
		delete[]Q;
		delete[]wB;
		delete[]RB;

		return true;
	}

	//改变特征值或特征向量来嵌入水印
	bool EigenDeformation::embedByL()
	{
		WaterMark wm;
		wm.setM((int)ceil((double)m_vertexNum / chip_rate));//设置原始水印位数
		wm.setC(chip_rate);//设置码片速率
		wm.setAlpha(0.005);
		wm.setKey(7);

		wm.createA();
		wm.createWB();
		wm.createP();

		E_matrix.setZero(m_vertexNum, m_vertexNum);
		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix是特征向量矩阵，nxn
		}

		for (int i = 0; i < m_vertexNum*0.75; i++)
		{
			E_matrix(i, i) = E_matrix(i, i) + wm.vecB[i] * wm.P[i] * wm.alpha;
			E_matrix(i, i + 1) = E_matrix(i, i + 1) + wm.vecB[i] * wm.P[i] * wm.alpha;
			E_matrix(i + 1, i) = E_matrix(i + 1, i) + wm.vecB[i] * wm.P[i] * wm.alpha;
			//E_matrix(i-1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
			//E_matrix(i+1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
		}

		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			m_normVector[e_index] = E_matrix.col(e_index);//E_matrix是特征向量矩阵，nxn
		}

		//将修改后的单位化的特征向量写入文件中
		ofstream WEfile;
		WEfile.open("D:\\code\\GeometryProcessing-1\\Txt\\WNe.txt", ios_base::out);
		if (WEfile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				for (int j = 0; j < m_vertexNum; j++)
				{
					WEfile << E_matrix(i, j) << " ";
				}
			}
		}
		WEfile.close();

		//频谱系数坐标
		VectorXd Rs(m_vertexNum);
		VectorXd Rt(m_vertexNum);
		VectorXd Ru(m_vertexNum);

		Rs = R_matrix.col(0);
		Rt = R_matrix.col(1);
		Ru = R_matrix.col(2);

		//反变换后的顶点坐标
		VectorXd newxCoordinate(m_vertexNum);
		VectorXd newyCoordinate(m_vertexNum);
		VectorXd newzCoordinate(m_vertexNum);

		newxCoordinate.setZero();
		newyCoordinate.setZero();
		newzCoordinate.setZero();

		VectorXd  normV(m_vertexNum);//Eigen类型的数组，用作临时变量
		for (int i_for_set = 0; i_for_set< m_vertexNum; i_for_set++)
		{
			normV = m_normVector[i_for_set];
			newxCoordinate = newxCoordinate + Rs[i_for_set] * normV;//X = ∑Rs[i] . ei
			newyCoordinate = newyCoordinate + Rt[i_for_set] * normV;
			newzCoordinate = newzCoordinate + Ru[i_for_set] * normV;
		}


		//将修改后的顶点坐标赋回给网格上的顶点
		int vi = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, vi++)
		{
			auto& t_p = m_oriMesh->point(*vit);

			t_p[0] = newxCoordinate[vi];
			t_p[1] = newyCoordinate[vi];
			t_p[2] = newzCoordinate[vi];
		}

		m_oriMesh->update_face_normals();
		//int eig_index = 0;
		//double *eigVec = new double[m_vertexNum*m_vertexNum];
		//for (int axi_col = 0; axi_col < m_vertexNum; axi_col++ )
		//{
		//	for ( int axi_row =0; axi_row < m_vertexNum; axi_row++)
		//	{	
		//		eigVec[eig_index++] = E_matrix(axi_row,axi_col);
		//	}
		//}


		//double *PN = new double[m_vertexNum];
		//for (int val_index = 0; val_index < m_vertexNum; val_index++ )
		//{
		//	PN[val_index] = m_eigenValue[val_index];
		//} 


		//if((!m_engine && !(m_engine = engOpen(NULL))))// 使用matlab算新的laplacian矩阵及特征值
		//{
		//	return false;
		//}
		//engSetVisible(m_engine,1);	

		//mxArray *pN = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);
		//mxArray *Eigen_vetor= mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);
		//
		//memcpy((void *) mxGetPr(pN), (void *) PN, m_vertexNum * sizeof(double));
		//memcpy((void *) mxGetPr(Eigen_vetor), (void *) eigVec,  (m_vertexNum * m_vertexNum) * sizeof(double));

		//engEvalString(m_engine, "cd('D:\\code\\GeometryProcessing-1\\mat_code')");//shenzhen

		//engPutVariable(m_engine, "C", pN);//矩阵方程为LX = λX, L =λX/X, 其中λ是特征值 = C ，X是特征向量 = D

		//engPutVariable(m_engine, "D",Eigen_vetor);

		//   engEvalString(m_engine,"[La,resC] = calcNewL(C, D);");

		////buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
		//char buffer[255];
		//buffer[254] = '\0';
		//engOutputBuffer(m_engine, buffer, 255);

		//printf("%s", buffer);//当matlab代码出错时，用来输出调试信息
		//
		//mxArray * l_matrix = NULL;
		//l_matrix = engGetVariable(m_engine, "La");
		//double  * Lmatrix = NULL;
		//Lmatrix = (double*)mxGetData(l_matrix);

		//MatrixXX  newL;
		//newL.setZero(m_vertexNum,m_vertexNum);

		//int l_index = 0;
		//for (  int l_col = 0; l_col < m_vertexNum; l_col ++)
		//{
		//	for ( int l_row = 0; l_row < m_vertexNum; l_row++)
		//	{
		//		newL(l_row, l_col) = Lmatrix[l_index++];
		//	}
		//}
		return true;
	}

	bool EigenDeformation::extractByL()
	{
		calLap_Matrix();
		normVec();

		E_matrix.setZero(m_vertexNum, m_vertexNum);
		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix是特征向量矩阵，nxn
		}


		////从文件中将原始网格的E矩阵读取出来
		double *pE = new double[m_vertexNum * m_vertexNum];
		ifstream NERfile;
		NERfile.open("D:\\code\\GeometryProcessing-1\\Txt\\WNe.txt", ios_base::in);
		if (NERfile)
		{
			// read into memory
			NERfile.seekg(0, NERfile.end);
			int length = NERfile.tellg();
			NERfile.seekg(0, NERfile.beg);

			char *buffer = new char[length];
			NERfile.read(buffer, length);
			NERfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (iss >> pE[i++]);
			delete[] buffer;
			// print or use it.
		}
		NERfile.close();

		MatrixXX	RE_matrix;
		RE_matrix.setZero(m_vertexNum, m_vertexNum);
		int ecount = 0;
		for (int i = 0; i < m_vertexNum; i++)
		{
			for (int j = 0; j < m_vertexNum; j++)
			{
				RE_matrix(i, j) = pE[ecount++];
			}
		}

		//将P读取出来
		double *P = new double[m_vertexNum];
		ifstream Pfile;
		Pfile.open("D:\\code\\GeometryProcessing-1\\Txt\\P.txt", ios_base::in);
		if (Pfile)
		{
			// read into memory
			Pfile.seekg(0, Pfile.end);
			int length = Pfile.tellg();
			Pfile.seekg(0, Pfile.beg);

			char *buffer = new char[length];
			Pfile.read(buffer, length);
			Pfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (i < m_vertexNum)
			{
				iss >> P[i++];
			}
			delete[] buffer;
			// print or use it.
		}
		Pfile.close();

		double * result = new double[m_vertexNum / 2];
		for (int ei = 0; ei<m_vertexNum / 2; ei++)
		{
			result[ei] = chip_rate*(RE_matrix(ei, ei) - E_matrix(ei, ei))*P[ei];
		}

		double * wB = new double[m_vertexNum];
		for (int i = 0; i < m_vertexNum / 2; i++)
		{
			wB[i] = _copysign(1.0, result[i]);
		}

		//将提原始的水印序列读出
		double *RB = new double[m_vertexNum];
		ifstream Wbfile;
		Wbfile.open("D:\\code\\GeometryProcessing-1\\Txt\\Wb.txt", ios_base::in);
		if (Wbfile)
		{
			// read into memory
			Wbfile.seekg(0, Wbfile.end);
			int length = Wbfile.tellg();
			Wbfile.seekg(0, Wbfile.beg);

			char *buffer = new char[length];
			Wbfile.read(buffer, length);
			Wbfile.close();

			// parse into array
			std::istringstream iss(buffer);
			int i = 0;
			while (i<m_vertexNum)
			{
				iss >> RB[i++];
			}
			delete[] buffer;
			// print or use it.
		}
		Wbfile.close();


		//比较读取出的Wb和提取出的wB
		double Bcounter = 0.0;
		for (int bi = 0; bi < m_vertexNum / 2; bi++)
		{
			if (RB[bi] == wB[bi])
			{
				Bcounter++;
			}
		}

		double corr = Bcounter / (m_vertexNum / 2);
		cout << "The result is: " << corr << endl;

		return true;
	}

}
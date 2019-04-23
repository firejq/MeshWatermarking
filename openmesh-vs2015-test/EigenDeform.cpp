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
		m_engine = NULL;//��ʼ��matlab����
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

	//	for (auto v_it = m_oriMesh->vertices_begin(); v_it != m_oriMesh->vertices_end(); ++v_it)//��ʼ����������
	//	{
	//		auto  point = m_oriMesh->point(v_it.handle());//ԭʼ�����
	//		m_verPostion(vex_index, 0) = point[0];
	//		m_verPostion(vex_index, 1) = point[1];
	//		m_verPostion(vex_index, 2) = point[2];
	//		vex_index++;
	//	}
	//	//m_verPos_iter = m_verPostion;//����ʼ����ֵ  
	//}

	//copy from yuan ʹ���������������˹����
	bool EigenDeformation::calLapMatrix()
	{
		const int nEdges = m_oriMesh->n_edges();
		const int nFaces = m_oriMesh->n_faces();
		const int totalElements = m_vertexNum + nEdges * 2;//laplacian����ĶԽ�����Ԫ��Ϊ����Ķȣ�ÿ�г��˶Խ���Ԫ�����

		double * xCoordinate = new double[m_vertexNum];
		double * yCoordinate = new double[m_vertexNum];
		double * zCoordinate = new double[m_vertexNum];

		double * px = xCoordinate;
		double * py = yCoordinate;
		double * pz = zCoordinate;
		//��ȡ�����ϵĵ������
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
			Tri[i] = tTri[(i%nFaces) * 3 + (i / nFaces)];//��Ϊ���е�����--û�п�����
		}
		if ((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab��laplacian��������ֵ
		{
			return false;
		}
		engSetVisible(m_engine, 1);
		mxArray *XCoordinate = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);//���������ݴ���matlab
		mxArray *YCoordinate = mxCreateDoubleMatrix(m_vertexNum, 1, mxREAL);//����m_vertexNum�У�1�е�ʵ����
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
		engEvalString(m_engine, "shape = struct('X',{X},'Y',{Y},'Z',{Z},'TRIV',{TRIV});");//add by shiqun,struct ��matlab�Դ��ĺ���


																						  //buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		//modified by shiqun
		/*
		ԭ�����Ǹ�matlabһ����״��Ȼ�����ǰ300������ֵ������������computeLacian�����ĵڶ���
		�������ֻ�������֣�������m_vertexNum������ĺ������õ���matlab�Ĵ��룬������Ĳ���
		��һ����ֵ��������һ������
		*/
		engEvalString(m_engine, "[eigVector,eigValue,W,D] = computeLaplacian(shape,233);");

		printf("%s", buffer);//��matlab�������ʱ���������������Ϣ


		mxArray * eigenVectorObj = NULL;
		eigenVectorObj = engGetVariable(m_engine, "eigVector");
		mxArray * eigenValueObj = NULL;
		eigenValueObj = engGetVariable(m_engine, "eigValue");

		/*add by shiqun
		eigenValue������һ��һά���󣬾���Ԫ�ظ���Ϊ1180��ÿ��Ԫ���Ƕ�Ӧ������ֵ
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

		for (int eig_index = 0; eig_index < m_vertexNum; eig_index++)//��eigenvector����m_eigenVector
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
		MatrixXX D_matrix;//nxn�Ⱦ���
		MatrixXX A_matrix;//nxn�ڽӾ���
		D_matrix.setZero(m_vertexNum, m_vertexNum);
		A_matrix.setZero(m_vertexNum, m_vertexNum);
		Lap_matrix.setZero(m_vertexNum, m_vertexNum);//nxn������˹����

													 //����Ⱦ���
		int v_index = 0;
		for (auto v_it = m_oriMesh->vertices_begin(); v_it != m_oriMesh->vertices_end(); ++v_it, v_index++)
		{
			int vnum = m_oriMesh->valence(*v_it);//��ö���i�Ķ�
			D_matrix(v_index, v_index) = vnum;
		}

		//�����ڽӾ���
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

		//����������˹����
		Lap_matrix = D_matrix - A_matrix;

		// ʹ��matlab��laplacian��������ֵ
		int l_index = 0;
		double *PL = new double[m_vertexNum*m_vertexNum];//@firejq:???
		for (int l_col = 0; l_col < m_vertexNum; l_col++)
		{
			for (int l_row = 0; l_row < m_vertexNum; l_row++)
			{
				PL[l_index++] = Lap_matrix(l_row, l_col);
			}
		}

		if ((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab��laplacian��������ֵ
		{
			std::cout << "Matlab engine inits failed" << std::endl;
			return;
		}
		engSetVisible(m_engine, 1);

		mxArray *pL = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);
		memcpy((void *)mxGetPr(pL), (void *)PL, m_vertexNum*m_vertexNum * sizeof(double));
		engPutVariable(m_engine, "L", pL);

		engEvalString(m_engine, "cd('D:\\firejq\\����\\Watermarking\\mat_code')");//����matalab����

																				//buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "[eigVector,eigValue] = calcLaplacian(L);");
		//engEvalString(m_engine,"[eigVector,eigValue] = eigenDeform(L);");//ֻ����ǰ6����������ֵ

		printf("%s", buffer);//��matlab�������ʱ���������������Ϣ

		mxArray * eigenVectorObj = NULL;
		eigenVectorObj = engGetVariable(m_engine, "eigVector");
		mxArray * eigenValueObj = NULL;
		eigenValueObj = engGetVariable(m_engine, "eigValue");

		/*add by shiqun
		eigenValue������һ��һά���󣬾���Ԫ�ظ���Ϊ1180��ÿ��Ԫ���Ƕ�Ӧ������ֵ
		*/
		double  * eigenVector = NULL;
		eigenVector = (double*)mxGetData(eigenVectorObj);
		double  * eigenValue = NULL;
		eigenValue = (double*)mxGetData(eigenValueObj);
		//		
		RowVectorX eVector_iter;
		eVector_iter.resize(m_vertexNum);

		//��ȡ��������
		for (int eig_index = 0; eig_index < m_vertexNum; eig_index++)//��eigenvector����m_eigenVector
		{
			eVector_iter.setZero();
			for (int ver_pos = 0; ver_pos < m_vertexNum; ver_pos++)
			{
				eVector_iter[ver_pos] = eigenVector[eig_index * m_vertexNum + ver_pos];
			}
			m_eigenVector.push_back(eVector_iter);
		}

		//��ȡ����ֵ
		for (int val_index = 0; val_index < m_vertexNum; val_index++)
		{
			m_eigenValue.push_back(eigenValue[val_index]);
		}
	}
	//add by shiqun
	//������������λ��

	void::EigenDeformation::normVec()
	{
		m_normVector = m_eigenVector;//�������������Ƹ���Ҫ��λ��������
		VectorXd  normV(m_vertexNum);//����һ��Eigen�е�vector���͵���ʱ�����������洢һ����������������
									 //ÿ�ν�һ�������������Ƹ��м����normV��Ȼ��normV���õ�λ���������ٽ�normV����m_normVector
		for (int index = 0; index < m_vertexNum; index++)//���������д����index<1�����ִ�����ɵĴ������ǹ������ˡ���
		{
			normV = m_normVector[index];
			normV.normalize();
			m_normVector.push_back(normV);
		}
	}

	//add by shiqun
	//����������ͶӰ������������,�޸�Ƶ��ϵ��������Ƶ��ϵ�����任��ԭʼ������
	//����Ƶ��ϵ��
	bool EigenDeformation::calR_Matrix()
	{
		//��������
		VectorXd xCoordinate(m_vertexNum);
		VectorXd yCoordinate(m_vertexNum);
		VectorXd zCoordinate(m_vertexNum);

		//���������긳ֵ��VectorXd
		int j = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, j++)
		{
			auto t_p = m_oriMesh->point(*vit);

			xCoordinate[j] = t_p[0];
			yCoordinate[j] = t_p[1];
			zCoordinate[j] = t_p[2];
		}

		//���㶥���������V_matrix
		V_matrix.setZero(m_vertexNum, 3);
		V_matrix.col(0) = xCoordinate;//V_matrix�Ƕ����������nx3
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

		//����������������E_matrix
		E_matrix.setZero(m_vertexNum, m_vertexNum);
		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix��������������nxn
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

		//����λ������������д���ļ���
		ofstream NormEfile;
		NormEfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Ne.txt", ios_base::out);
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

		/*����matlab����Ƶ��ϵ������R_matrix*/
		mxArray *Vetex_Coor = mxCreateDoubleMatrix(m_vertexNum, 3, mxREAL);//�����񶥵����������matlab
		mxArray *Eigen_vetor = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);

		memcpy((void *)mxGetPr(Eigen_vetor), (void *)eigVec, (m_vertexNum * m_vertexNum) * sizeof(double));
		memcpy((void *)mxGetPr(Vetex_Coor), (void *)vetex_coor, (m_vertexNum * 3) * sizeof(double));

		engPutVariable(m_engine, "A", Eigen_vetor);//���󷽳�ΪAX = B,����A = E_matrix�� X = R_matrix�� B = V_matrix
		engPutVariable(m_engine, "B", Vetex_Coor);

		//buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "cd('D:\\firejq\\����\\Watermarking\\mat_code')");//shenzhen
		engEvalString(m_engine, "R = calcR(A, B);");

		printf("%s", buffer);//��matlab�������ʱ���������������Ϣ

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
		// ����ʱ���ùر�
		/*if (m_engine)
		{
		engClose(m_engine);
		m_engine = NULL;
		}*/

		return true;
	}

	bool EigenDeformation::embedWatermark()
	{
		//Ƶ��ϵ������
		VectorXd Rs(m_vertexNum);
		VectorXd Rt(m_vertexNum);
		VectorXd Ru(m_vertexNum);
		Rs = R_matrix.col(0);
		Rt = R_matrix.col(1);
		Ru = R_matrix.col(2);

		//�޸ĺ��Ƶ��ϵ������
		VectorXd Rs_hat(m_vertexNum);
		VectorXd Rt_hat(m_vertexNum);
		VectorXd Ru_hat(m_vertexNum);

		//���任��Ķ�������
		VectorXd newxCoordinate(m_vertexNum);
		VectorXd newyCoordinate(m_vertexNum);
		VectorXd newzCoordinate(m_vertexNum);
		newxCoordinate.setZero();
		newyCoordinate.setZero();
		newzCoordinate.setZero();

		WaterMark mytest;
		mytest.setM((int)ceil((double)m_vertexNum / chip_rate));//����ԭʼˮӡλ��
		mytest.setC(chip_rate);//������Ƭ����
		mytest.setAlpha(0.005);//��
		mytest.setKey(7);//��

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
		��ԭʼ�����Ƶ��ϵ��R����д���ļ��У�
		������ȡˮӡ��ʱ��Ͳ����ٶ���ԭʼ����ֻ�����������ˮӡ�����е�wR����
		Ȼ����ļ��ж�ȡ��R������������ֵ���ٽ����жϲ��輴�ɡ�
		**/
		//Rs
		ofstream Rsfile;
		Rsfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Rs.txt", ios_base::out);
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
		Rtfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Rt.txt", ios_base::out);
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
		Rufile.open("D:\\firejq\\����\\Watermarking\\Txt\\Ru.txt", ios_base::out);
		if (Rufile)
		{
			for (int i = 0; i < m_vertexNum; i++)
			{
				Rufile << Ru[i] << " ";
			}
		}
		Rufile.close();

		VectorXd  normV(m_vertexNum);//Eigen���͵����飬������ʱ����
		for (int i_for_set = 0; i_for_set< m_vertexNum; i_for_set++)
		{
			normV = m_normVector[i_for_set];
			newxCoordinate = newxCoordinate + Rs_hat[i_for_set] * normV;//X = ��Rs[i] . ei
			newyCoordinate = newyCoordinate + Rt_hat[i_for_set] * normV;
			newzCoordinate = newzCoordinate + Ru_hat[i_for_set] * normV;
		}


		//���޸ĺ�Ķ������긳�ظ������ϵĶ���
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

	//��ȡˮӡ
	//��Ƕ��֮ǰ��ˮӡ��Ϣд���ļ��У�
	//��ȡˮӡʱ������ȡ��������Ҳд���ļ��У�Ȼ��Ƚ������ļ��е���Ϣ�Ƿ�һ��
	bool EigenDeformation::extractWatermark()
	{
		//ԭʼ������ˮӡ����������ֵ�ֽ⣬Ȼ������ֵ���

		//��������
		VectorXd wxCoordinate(m_vertexNum);
		VectorXd wyCoordinate(m_vertexNum);
		VectorXd wzCoordinate(m_vertexNum);

		//���������긳ֵ��VectorXd
		int k = 0;
		for (auto vit = m_oriMesh->vertices_begin(); vit != m_oriMesh->vertices_end(); ++vit, k++)
		{
			auto t_p = m_oriMesh->point(*vit);

			wxCoordinate[k] = t_p[0];
			wyCoordinate[k] = t_p[1];
			wzCoordinate[k] = t_p[2];
		}

		V_matrix.setZero(m_vertexNum, 3);

		V_matrix.col(0) = wxCoordinate;//V_matrix�Ƕ����������nx3
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
		//����matlab������������
		//for ( int e_index = 0; e_index< m_vertexNum; e_index++)
		//{ 
		//	E_matrix.col(e_index) = m_normVector[e_index];//E_matrix��������������nxn
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

		////���ļ��н�ԭʼ�����E�����ȡ����
		double *pE = new double[m_vertexNum * m_vertexNum];
		ifstream NERfile;
		NERfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Ne.txt", ios_base::in);
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

		//����matlab
		if ((!m_engine && !(m_engine = engOpen(NULL))))
		{
			return false;
		}
		engSetVisible(m_engine, 1);

		mxArray *Vetex_Coor = mxCreateDoubleMatrix(m_vertexNum, 3, mxREAL);//�����񶥵����������matlab
		mxArray *Eigen_vetor = mxCreateDoubleMatrix(m_vertexNum, m_vertexNum, mxREAL);

		memcpy((void *)mxGetPr(Eigen_vetor), (void *)eigVec, (m_vertexNum * m_vertexNum) * sizeof(double));
		memcpy((void *)mxGetPr(Vetex_Coor), (void *)vetex_coor, (m_vertexNum * 3) * sizeof(double));

		engPutVariable(m_engine, "A", Eigen_vetor);//���󷽳�ΪAX = B,����A = E_matrix�� X = R_matrix�� B = V_matrix
		engPutVariable(m_engine, "B", Vetex_Coor);

		//buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
		char buffer[255];
		buffer[254] = '\0';
		engOutputBuffer(m_engine, buffer, 255);

		engEvalString(m_engine, "cd('D:\\firejq\\����\\Watermarking\\mat_code')");//shenzhen
		engEvalString(m_engine, "R = calcR(A, B);");

		printf("%s", buffer);//��matlab�������ʱ���������������Ϣ

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
		// ����ʱ���ùر�
		if (m_engine)
		{
			engClose(m_engine);
			m_engine = NULL;
		}

		//Ƶ��ϵ������
		VectorXd wRs(m_vertexNum);
		VectorXd wRt(m_vertexNum);
		VectorXd wRu(m_vertexNum);

		wRs = R_matrix.col(0);
		wRt = R_matrix.col(1);
		wRu = R_matrix.col(2);

		//��ԭʼ�����Ƶ��ϵ�����ļ��ж�ȡ����
		double *pRs = new double[m_vertexNum];
		double *pRt = new double[m_vertexNum];
		double *pRu = new double[m_vertexNum];

		//Rs
		ifstream Rsfile;
		Rsfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Rs.txt", ios_base::in);
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
		Rtfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Rt.txt", ios_base::in);
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
		Rufile.open("D:\\firejq\\����\\Watermarking\\Txt\\Ru.txt", ios_base::in);
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

		//��P��ȡ����
		double *P = new double[m_vertexNum];
		ifstream Pfile;
		Pfile.open("D:\\firejq\\����\\Watermarking\\Txt\\P.txt", ios_base::in);
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

		//@firejq: ����ȡ����ˮӡ����д���ļ����������ʱ�鿴
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



		//����ԭʼ��ˮӡ���ж���
		double *RB = new double[m_vertexNum];
		ifstream Wbfile;
		Wbfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Wb.txt", ios_base::in);
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

		//�Ƚ϶�ȡ����Wb����ȡ����wB
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


		//����ȡ����wB�����ļ�
		//��ת��Ϊ��ֵ����
		for (int k = 0; k < 1186; k++)
		{
			if (wB[k] == -1)
			{
				wB[k] = 0;
			}
		}
		ofstream EWbfile;
		EWbfile.open("D:\\firejq\\����\\Watermarking\\Txt\\AMesh\\Extract\\FlowerWb.txt", ios_base::out);
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

	//�ı�����ֵ������������Ƕ��ˮӡ
	bool EigenDeformation::embedByL()
	{
		WaterMark wm;
		wm.setM((int)ceil((double)m_vertexNum / chip_rate));//����ԭʼˮӡλ��
		wm.setC(chip_rate);//������Ƭ����
		wm.setAlpha(0.005);
		wm.setKey(7);

		wm.createA();
		wm.createWB();
		wm.createP();

		E_matrix.setZero(m_vertexNum, m_vertexNum);
		for (int e_index = 0; e_index< m_vertexNum; e_index++)
		{
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix��������������nxn
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
			m_normVector[e_index] = E_matrix.col(e_index);//E_matrix��������������nxn
		}

		//���޸ĺ�ĵ�λ������������д���ļ���
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

		//Ƶ��ϵ������
		VectorXd Rs(m_vertexNum);
		VectorXd Rt(m_vertexNum);
		VectorXd Ru(m_vertexNum);

		Rs = R_matrix.col(0);
		Rt = R_matrix.col(1);
		Ru = R_matrix.col(2);

		//���任��Ķ�������
		VectorXd newxCoordinate(m_vertexNum);
		VectorXd newyCoordinate(m_vertexNum);
		VectorXd newzCoordinate(m_vertexNum);

		newxCoordinate.setZero();
		newyCoordinate.setZero();
		newzCoordinate.setZero();

		VectorXd  normV(m_vertexNum);//Eigen���͵����飬������ʱ����
		for (int i_for_set = 0; i_for_set< m_vertexNum; i_for_set++)
		{
			normV = m_normVector[i_for_set];
			newxCoordinate = newxCoordinate + Rs[i_for_set] * normV;//X = ��Rs[i] . ei
			newyCoordinate = newyCoordinate + Rt[i_for_set] * normV;
			newzCoordinate = newzCoordinate + Ru[i_for_set] * normV;
		}


		//���޸ĺ�Ķ������긳�ظ������ϵĶ���
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


		//if((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab���µ�laplacian��������ֵ
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

		//engPutVariable(m_engine, "C", pN);//���󷽳�ΪLX = ��X, L =��X/X, ���Ц�������ֵ = C ��X���������� = D

		//engPutVariable(m_engine, "D",Eigen_vetor);

		//   engEvalString(m_engine,"[La,resC] = calcNewL(C, D);");

		////buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
		//char buffer[255];
		//buffer[254] = '\0';
		//engOutputBuffer(m_engine, buffer, 255);

		//printf("%s", buffer);//��matlab�������ʱ���������������Ϣ
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
			E_matrix.col(e_index) = m_normVector[e_index];//E_matrix��������������nxn
		}


		////���ļ��н�ԭʼ�����E�����ȡ����
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

		//��P��ȡ����
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

		//����ԭʼ��ˮӡ���ж���
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


		//�Ƚ϶�ȡ����Wb����ȡ����wB
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
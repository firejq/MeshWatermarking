#include "MeshWatermark.h"

MeshWatermark::MeshWatermark()
{
	m_oriMesh = NULL;
	m_engine = NULL;//��ʼ��matlab����
	//m_bdyCenter.setZero();
	//m_meanDis = 0.0;//copy from yuan for segementation
}

void MeshWatermark::Init(PolygonMesh::Mesh * _mesh, string watermark_path, int c, double a)
{
	chip_rate = c;
	alpha = a;
	m_oriMesh = _mesh;
	//m_wmMesh =_mesh;
	m_vertexNum = m_oriMesh->n_vertices();
	wm_path = watermark_path;
}

void MeshWatermark::Init(PolygonMesh::Mesh * _mesh, string watermark_path)
{
	chip_rate = 7;
	alpha = 0.005;
	m_oriMesh = _mesh;
	//m_wmMesh =_mesh;
	m_vertexNum = m_oriMesh->n_vertices();
	wm_path = watermark_path;
}

//L = D - A
void MeshWatermark::calLap_Matrix()
{
	std::cout << "Calculating Laplacian Matrix..." << std::endl;

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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//����matalab����

																			//buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
	char buffer[255];
	buffer[254] = '\0';
	engOutputBuffer(m_engine, buffer, 255);

	engEvalString(m_engine, "[eigVector,eigValue] = calcLaplacian(L);");

	printf("%s", buffer);//��matlab�������ʱ���������������Ϣ

	mxArray * eigenVectorObj = NULL;
	eigenVectorObj = engGetVariable(m_engine, "eigVector");
	mxArray * eigenValueObj = NULL;
	eigenValueObj = engGetVariable(m_engine, "eigValue");

	/*
	eigenValue������һ��һά���󣬾���Ԫ�ظ���Ϊ1180��ÿ��Ԫ���Ƕ�Ӧ������ֵ
	*/
	double  * eigenVector = NULL;
	eigenVector = (double*)mxGetData(eigenVectorObj);
	double  * eigenValue = NULL;
	eigenValue = (double*)mxGetData(eigenValueObj);

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

//������������λ��
void::MeshWatermark::normVec()
{
	std::cout << "Normalizing Eigen Vectors..." << std::endl;

	m_normVector = m_eigenVector;//�������������Ƹ���Ҫ��λ��������
	VectorXd  normV(m_vertexNum);//����һ��Eigen�е�vector���͵���ʱ�����������洢һ����������������

	//ÿ�ν�һ�������������Ƹ��м����normV��Ȼ��normV���õ�λ���������ٽ�normV����m_normVector
	for (int index = 0; index < m_vertexNum; index++)
	{
		normV = m_normVector[index];
		normV.normalize();
		m_normVector.push_back(normV);
	}
}

//����������ͶӰ�����������ϣ��޸�Ƶ��ϵ����Ȼ��Ƶ��ϵ�����任��ԭʼ������
//����Ƶ��ϵ��
bool MeshWatermark::calR_Matrix()
{
	std::cout << "Calculating spectral coefficient matrix..." << std::endl;

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
	for (int e_index = 0; e_index < m_vertexNum; e_index++)
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
	NormEfile.open(ROOT_PATH + "file\\tmp\\ori_NormEigenVector.txt", ios_base::out);
	if (NormEfile)
	{
		for (int i = 0; i < m_vertexNum; i++)
		{
			for (int j = 0; j < m_vertexNum; j++)
			{
				NormEfile << E_matrix(i, j) << " ";
			}
			//NormEfile << "\n";
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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//shenzhen
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

bool MeshWatermark::embedWatermark()
{
	calLap_Matrix(); // @firejq: ���������������˹���󲢵���matlab��������ֵ�ֽ�
	normVec();//������������λ��
	calR_Matrix();//@firejq: ����������ͶӰ�����������ϵõ�Ƶ��ϵ������Ȼ�����ˮӡ�޸�Ƶ��ϵ�����ٽ�Ƶ��ϵ�����任��ԭʼ������

	std::cout << "Embing watermark..." << std::endl;

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

	// ����ˮӡ����
	WatermarkSeq wmSeq;
	wmSeq.setM((int)ceil((double)m_vertexNum / chip_rate));//����ԭʼˮӡλ��
	wmSeq.setC(chip_rate);//������Ƭ����
	wmSeq.setAlpha(alpha);//
	wmSeq.setKey(7);//
	wmSeq.initWB(wm_path);//��ʼ��ˮӡͼ�����Ϣ
	wmSeq.createWB(m_vertexNum);//����ˮӡ����
	wmSeq.createP();//����α������� (PRNS)

	//��ˮӡ���и�ֵ��Ƶ��ϵ������Rs/Rt/Ru�У��õ��޸ĺ��Ƶ��ϵ������Rs_hat/Rt_hat/Ru_hat
	//vecB�ĳ�����1188����m_vertexNum��һ����������len_VecB>len_m_vertexNum���ᶪʧһ����ˮӡ��Ϣ�����Ҫ��len_VecB<=len_m_vertexNum
	for (int i_for_get = 0; i_for_get < m_vertexNum; i_for_get++)
	{
		Rs_hat[i_for_get] = Rs[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;//Rs[i]' = Rs[i] + bi' . pi . a 
		Rt_hat[i_for_get] = Rt[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;
		Ru_hat[i_for_get] = Ru[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;
	}

	/**
	��ԭʼ�����Ƶ��ϵ��R����д���ļ��У�
	������ȡˮӡ��ʱ��Ͳ����ٶ���ԭʼ����ֻ�����������ˮӡ�����е�wR����
	Ȼ����ļ��ж�ȡ��R������������ֵ���ٽ����жϲ��輴�ɡ�
	**/
	//Rs
	ofstream Rsfile;
	Rsfile.open(ROOT_PATH + "file\\tmp\\ori_Rs.txt", ios_base::out);
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
	Rtfile.open(ROOT_PATH + "file\\tmp\\ori_Rt.txt", ios_base::out);
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
	Rufile.open(ROOT_PATH + "file\\tmp\\ori_Ru.txt", ios_base::out);
	if (Rufile)
	{
		for (int i = 0; i < m_vertexNum; i++)
		{
			Rufile << Ru[i] << " ";
		}
	}
	Rufile.close();
	/*------------------------------------------------*/


	/*���޸ĺ�Ķ������긳�ظ������ϵĶ���*/
	VectorXd  normV(m_vertexNum);//Eigen���͵����飬������ʱ����
	for (int i_for_set = 0; i_for_set < m_vertexNum; i_for_set++)
	{
		normV = m_normVector[i_for_set];
		newxCoordinate = newxCoordinate + Rs_hat[i_for_set] * normV;//X = ��Rs[i] . ei
		newyCoordinate = newyCoordinate + Rt_hat[i_for_set] * normV;
		newzCoordinate = newzCoordinate + Ru_hat[i_for_set] * normV;
	}

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
bool MeshWatermark::extractWatermark(string extr_watermark)
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

	//���ļ��н�ԭʼ�����E�����ȡ����
	std::cout << "Reading NormEigenVector of original mesh..." << std::endl;
	double *pE = new double[m_vertexNum * m_vertexNum];
	ifstream NERfile;
	NERfile.open(ROOT_PATH + "file\\tmp\\ori_NormEigenVector.txt", ios_base::in);
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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//shenzhen
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

	std::cout << "Reading spectral coefficient matrix of original mesh..." << std::endl;
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
	Rsfile.open(ROOT_PATH + "file\\tmp\\ori_Rs.txt", ios_base::in);
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
	Rtfile.open(ROOT_PATH + "file\\tmp\\ori_Rt.txt", ios_base::in);
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
	Rufile.open(ROOT_PATH + "file\\tmp\\ori_Ru.txt", ios_base::in);
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

	std::cout << "Reading PRNS used for embing..." << std::endl;
	//��P��ȡ����
	double *P = new double[m_vertexNum];
	ifstream Pfile;
	Pfile.open(ROOT_PATH + "file\\tmp\\emb_PRNS.txt", ios_base::in);
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


	std::cout << "Calculating the extracted watermark sequence..." << std::endl;
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
	Wb_extract.open(ROOT_PATH + "file\\tmp\\extr_WmSeq.txt", ios_base::out);
	if (Wb_extract)
	{
		for (int i = 0; i < m_vertexNum; i++)
		{
			Wb_extract << wB[i] << " ";
		}
	}
	Wb_extract.close();
	std::cout << "The watermark sequence is extracted successfully and saved in " +
		ROOT_PATH + "file\\extr_WmSeq.txt" << std::endl;

	//��ԭʼ��ˮӡ���ж���
	double *RB = new double[m_vertexNum];
	ifstream Wbfile;
	Wbfile.open(ROOT_PATH + "file\\tmp\\emb_WmSeq.txt", ios_base::in);
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
		while (i < m_vertexNum)
		{
			iss >> RB[i++];
		}
		delete[] buffer;
		// print or use it.
	}
	Wbfile.close();

	//�Ƚ϶�ȡ����Wb����ȡ����wB
	std::cout << "Comparing embWatermarkSeq with extrWatermarkSeq..." << std::endl;
	double Bcounter = 0.0;
	for (int bi = 0; bi < m_vertexNum; bi++)
	{
		//std::cout << RB[bi] << ',' << wB[bi] << std::endl;
		if (RB[bi] == wB[bi])
		{
			Bcounter++;
		}
	}
	double corr = Bcounter / (double)m_vertexNum;
	//for debug
	//std::cout << Bcounter << std::endl;
	//std::cout << m_vertexNum << std::endl;
	//std::cout << corr << std::endl;

	std::cout << "The comparing result is: ��" << (corr == 1 ? "true" : "false") << "��" << std::endl;

	WatermarkSeq wmSeq;
	wmSeq.initWB(wm_path);//��ʼ��ˮӡͼ�����Ϣ

	//����ȡ����wB�����ļ�
	//��ת��Ϊ��ֵ����
	for (int k = 0; k < wmSeq.wmSeqLength; k++)
	{
		if (wB[k] == -1)
		{
			wB[k] = 0;
		}
	}
	ofstream EWbfile;
	EWbfile.open(extr_watermark, ios_base::out);
	int i = 0;
	if (EWbfile)
	{
		for (int r = 0; r < wmSeq.wmImgWidth; r++)
		{
			for (int c = 0; c < wmSeq.wmImgLength; c++)
			{
				if (i < m_vertexNum)
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
	std::cout << "The watermark image is saved in " + extr_watermark << std::endl;

	delete[] pRs;
	delete[] pRt;
	delete[] pRu;
	delete[] P;
	delete[] Qx;
	delete[] Qy;
	delete[] Qz;
	delete[] Q;
	delete[] wB;
	delete[] RB;

	return true;
}

/*�ڲ���WatermarkSeq����*/
//ԭʼˮӡλ��
void MeshWatermark::WatermarkSeq::setM(int m_m)
{
	//cout << m_m << endl;
	m = m_m;
}

//��Ƭ����
void MeshWatermark::WatermarkSeq::setC(int m_c)
{
	c = m_c;
}


void MeshWatermark::WatermarkSeq::setAlpha(double a)
{
	alpha = a;
}

void MeshWatermark::WatermarkSeq::setKey(int k)
{
	key = k;
}

void MeshWatermark::WatermarkSeq::initWB(string watermark_path)
{
	//�Ѹĳ��Զ���ȡˮӡ��ֵͼƬ�ĳ��Ϳ�����ˮӡ���г���
	wmImgPath = watermark_path;

	ifstream Wfile;
	Wfile.open(wmImgPath, ios_base::in);
	if (!Wfile)
	{
		std::cout << "Load watermark image failed!" << std::endl;
		exit(1);
	}
	// read into memory
	Wfile.seekg(0, Wfile.end);
	int length = Wfile.tellg();
	Wfile.seekg(0, Wfile.beg);

	char *buffer = new char[length];
	Wfile.read(buffer, length);
	Wfile.close();

	// parse into string
	std::istringstream iss(buffer);
	const string content = iss.str();
	//std::cout << "---" << std::endl;
	//std::cout << count(content.begin(), content.end(), ' ') << std::endl;//1188
	//std::cout << count(content.begin(), content.end(), '\n') << std::endl;//33
	//std::cout << "---" << std::endl;
	delete[] buffer;
	Wfile.close();

	wmSeqLength = count(content.begin(), content.end(), ' ');//1188
	wmImgWidth = count(content.begin(), content.end(), '\n');//36
	wmImgLength = wmSeqLength / wmImgWidth;

}

/*Ϊˮӡ����VecB��ֵ������д�뵽�ļ���*/
void MeshWatermark::WatermarkSeq::createWB(int m_vertexNum)
{
	if (wmSeqLength > m_vertexNum) {
		std::cout << "��warning��The length of watermark sequance is bigger than the number of mesh vertexes. "
			<< "The extracted result will be incomplete!" << std::endl;
	}

	double *B = new double[wmSeqLength];
	ifstream Wfile;
	//cout << watermark_path << endl;
	Wfile.open(wmImgPath, ios_base::in);
	if (Wfile)
	{
		// read into memory
		Wfile.seekg(0, Wfile.end);
		int length = Wfile.tellg();
		Wfile.seekg(0, Wfile.beg);

		char *buffer = new char[length];
		Wfile.read(buffer, length);
		Wfile.close();

		// parse into array
		std::istringstream iss(buffer);
		int i = 0;
		while (i < wmSeqLength)
		{
			iss >> B[i++];
		}
		delete[] buffer;
		// print or use it.
	}
	Wfile.close();

	for (int i = 0; i < wmSeqLength; i++)
	{
		//cout << B[i] << ",";
		vecB.push_back(B[i] == 0 ? -1 : 1);
	}

	//@firejq
	//��ˮӡ����С�ڶ�������������Ƚ�ˮӡ������1��ȫ��ʹ��ˮӡ���ȵ��ڶ������
	//����ȡˮӡ��ʱ����ȥ���ⲿ�ֲ�ȫ�ĳ��ȼ���
	if (wmSeqLength < m_vertexNum) {
		for (int i = 0; i < m_vertexNum - wmSeqLength; i++)
		{
			vecB.push_back(1);
		}
	}

	//��ˮӡ����д���ļ���
	ofstream Wbfile;
	Wbfile.open(ROOT_PATH + "file\\tmp\\emb_WmSeq.txt", ios_base::out);
	if (Wbfile)
	{
		int length = std::max(wmSeqLength, m_vertexNum);
		for (int i = 0; i < length; i++)
		{
			Wbfile << vecB[i] << " ";
		}
	}
	Wbfile.close();
}

/*Ϊα�������P��ֵ����д�뵽�ļ���*/
void MeshWatermark::WatermarkSeq::createP()
{
	srand(key);//���������ֵ
	for (int i = 0; i < m*c; i++)
	{
		P.push_back(rand() % 2 == 0 ? -1 : 1);
	}

	//��Pд���ļ���
	ofstream Pfile;
	Pfile.open(ROOT_PATH + "file\\tmp\\emb_PRNS.txt", ios_base::out);
	if (Pfile)
	{
		for (int i = 0; i < m*c; i++)
		{
			Pfile << P[i] << " ";
		}
	}
	Pfile.close();
}
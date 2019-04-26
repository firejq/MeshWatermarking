#include "MeshWatermark.h"

MeshWatermark::MeshWatermark()
{
	m_oriMesh = NULL;
	m_engine = NULL;//初始化matlab引擎
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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//调用matalab代码

																			//buffer用来接收调试信息，当matlab代码有错时，可以输出buffer查看错误信息
	char buffer[255];
	buffer[254] = '\0';
	engOutputBuffer(m_engine, buffer, 255);

	engEvalString(m_engine, "[eigVector,eigValue] = calcLaplacian(L);");

	printf("%s", buffer);//当matlab代码出错时，用来输出调试信息

	mxArray * eigenVectorObj = NULL;
	eigenVectorObj = engGetVariable(m_engine, "eigVector");
	mxArray * eigenValueObj = NULL;
	eigenValueObj = engGetVariable(m_engine, "eigValue");

	/*
	eigenValue矩阵是一个一维矩阵，矩阵元素个数为1180，每个元素是对应的特征值
	*/
	double  * eigenVector = NULL;
	eigenVector = (double*)mxGetData(eigenVectorObj);
	double  * eigenValue = NULL;
	eigenValue = (double*)mxGetData(eigenValueObj);

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

//将特征向量单位化
void::MeshWatermark::normVec()
{
	std::cout << "Normalizing Eigen Vectors..." << std::endl;

	m_normVector = m_eigenVector;//将特征向量复制给将要单位化的向量
	VectorXd  normV(m_vertexNum);//声明一个Eigen中的vector类型的临时变量，用来存储一个单独的特征向量

	//每次将一个特征向量复制给中间变量normV，然后normV调用单位化函数，再将normV赋给m_normVector
	for (int index = 0; index < m_vertexNum; index++)
	{
		normV = m_normVector[index];
		normV.normalize();
		m_normVector.push_back(normV);
	}
}

//将顶点坐标投影到特征向量上，修改频谱系数，然后将频谱系数反变换到原始坐标中
//计算频谱系数
bool MeshWatermark::calR_Matrix()
{
	std::cout << "Calculating spectral coefficient matrix..." << std::endl;

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
	for (int e_index = 0; e_index < m_vertexNum; e_index++)
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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//shenzhen
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

bool MeshWatermark::embedWatermark()
{
	calLap_Matrix(); // @firejq: 计算网格的拉普拉斯矩阵并调用matlab进行特征值分解
	normVec();//将特征向量单位化
	calR_Matrix();//@firejq: 将顶点坐标投影到特征向量上得到频谱系数矩阵，然后根据水印修改频谱系数，再将频谱系数反变换到原始坐标中

	std::cout << "Embing watermark..." << std::endl;

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

	// 创建水印序列
	WatermarkSeq wmSeq;
	wmSeq.setM((int)ceil((double)m_vertexNum / chip_rate));//设置原始水印位数
	wmSeq.setC(chip_rate);//设置码片速率
	wmSeq.setAlpha(alpha);//
	wmSeq.setKey(7);//
	wmSeq.initWB(wm_path);//初始化水印图像的信息
	wmSeq.createWB(m_vertexNum);//创建水印序列
	wmSeq.createP();//创建伪随机序列 (PRNS)

	//将水印序列赋值到频谱系数坐标Rs/Rt/Ru中，得到修改后的频谱系数坐标Rs_hat/Rt_hat/Ru_hat
	//vecB的长度是1188，和m_vertexNum不一样，因此如果len_VecB>len_m_vertexNum，会丢失一部分水印信息，因此要求len_VecB<=len_m_vertexNum
	for (int i_for_get = 0; i_for_get < m_vertexNum; i_for_get++)
	{
		Rs_hat[i_for_get] = Rs[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;//Rs[i]' = Rs[i] + bi' . pi . a 
		Rt_hat[i_for_get] = Rt[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;
		Ru_hat[i_for_get] = Ru[i_for_get] + wmSeq.vecB[i_for_get] * wmSeq.P[i_for_get] * wmSeq.alpha;
	}

	/**
	将原始网格的频谱系数R矩阵写入文件中，
	这样提取水印的时候就不必再读入原始网格，只需计算出读入的水印网格中的wR矩阵，
	然后从文件中读取出R矩阵，两者做差值，再进行判断步骤即可。
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


	/*将修改后的顶点坐标赋回给网格上的顶点*/
	VectorXd  normV(m_vertexNum);//Eigen类型的数组，用作临时变量
	for (int i_for_set = 0; i_for_set < m_vertexNum; i_for_set++)
	{
		normV = m_normVector[i_for_set];
		newxCoordinate = newxCoordinate + Rs_hat[i_for_set] * normV;//X = ∑Rs[i] . ei
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

//提取水印
//在嵌入之前将水印信息写入文件中，
//提取水印时，将提取出的数据也写入文件中，然后比较两个文件中的信息是否一致
bool MeshWatermark::extractWatermark(string extr_watermark)
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

	//从文件中将原始网格的E矩阵读取出来
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

	engEvalString(m_engine, ("cd('" + ROOT_PATH + "mat_code')").c_str());//shenzhen
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

	std::cout << "Reading spectral coefficient matrix of original mesh..." << std::endl;
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
	//将P读取出来
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

	//@firejq: 将提取出的水印序列写到文件，方便调试时查看
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

	//将原始的水印序列读出
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

	//比较读取出的Wb和提取出的wB
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

	std::cout << "The comparing result is: 【" << (corr == 1 ? "true" : "false") << "】" << std::endl;

	WatermarkSeq wmSeq;
	wmSeq.initWB(wm_path);//初始化水印图像的信息

	//将提取出的wB存入文件
	//先转换为二值序列
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

/*内部类WatermarkSeq定义*/
//原始水印位数
void MeshWatermark::WatermarkSeq::setM(int m_m)
{
	//cout << m_m << endl;
	m = m_m;
}

//码片速率
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
	//把改成自动读取水印二值图片的长和宽，计算水印序列长度
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

/*为水印序列VecB赋值，并将写入到文件中*/
void MeshWatermark::WatermarkSeq::createWB(int m_vertexNum)
{
	if (wmSeqLength > m_vertexNum) {
		std::cout << "【warning】The length of watermark sequance is bigger than the number of mesh vertexes. "
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
	//若水印长度小于顶点个数，这里先将水印长度用1补全，使得水印长度等于顶点个数
	//在提取水印的时候再去除这部分补全的长度即可
	if (wmSeqLength < m_vertexNum) {
		for (int i = 0; i < m_vertexNum - wmSeqLength; i++)
		{
			vecB.push_back(1);
		}
	}

	//将水印序列写入文件中
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

/*为伪随机序列P赋值，并写入到文件中*/
void MeshWatermark::WatermarkSeq::createP()
{
	srand(key);//随机数种子值
	for (int i = 0; i < m*c; i++)
	{
		P.push_back(rand() % 2 == 0 ? -1 : 1);
	}

	//将P写入文件中
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
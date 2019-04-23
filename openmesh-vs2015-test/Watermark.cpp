#include "WaterMark.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <numeric>
#include <algorithm>

#include <sys/timeb.h>
#include <ctime>
#include <climits>

#define AREA_WATERMARK 1
#define CLUSTER_NUM 28

using namespace GeometryProcess;
using namespace OpenMesh;
using namespace std;

unsigned char colorMapSimple[60] =
{
	255, 0, 0//��ɫ
	,0, 255, 0//��ɫ
	,0,0,255//��ɫ
	,255, 255,0//��ɫ
	,0,255,255//ˮ��
	,255,0,255//Ʒ��
	,255,128,64//��ɫ
	,153 ,0 ,51//����
	,102 ,51 ,102//����
	,255 ,153 , 102//��ɫ
	,102 ,0 , 51//�ƺ�
	,153 ,204 , 255//����
	,51 , 51  ,102//��ɫ
	,0,191,255//����ɫ
	,64,128,128//ī��
	,130,0,0//����
	,195,226,204
	,208,172,191
	,158,79,30
	,98,90,5
};

WaterMark::WaterMark()
{
}

WaterMark::~WaterMark()
{
}
//ԭʼˮӡλ��
void WaterMark::setM(int m_m)
{
	m = m_m;
}

int WaterMark::getM()
{
	return m;
}

//��Ƭ����
void WaterMark::setC(int m_c)
{
	c = m_c;
}

int WaterMark::getC()
{
	return c;
}

void WaterMark::setAlpha(double a)
{
	alpha = a;
}

void WaterMark::setKey(int k)
{
	key = k;
}

void WaterMark::Init(PolygonMesh::Mesh * _mesh)
{
	myMesh = _mesh;
	num_vtx = _mesh->n_vertices();
}

//void WaterMark::calCurvature()
//{
//	int vidx = 0;
//	for (auto v_it = myMesh->vertices_begin(); v_it != myMesh->vertices_end(); ++v_it, vidx++)
//	{
//		double Kmin = 0.0;
//		double Kmax = 0.0;
//		double total_area = 0.0;
//		Matrix33 tangMat;
//		tangMat.setZero();
//		for (auto oh_it = myMesh->voh_begin(v_it); oh_it != myMesh->voh_end(v_it); ++oh_it)
//		{
//			Matrix33 tempMat;
//			tempMat.setZero();
//			double area = 0.0;
//			double ang = 0.0;
//			calEdgeVec(oh_it, tempMat);//����|e| * (e * eT)
//			calAngleTris(oh_it, ang, area);//�����߶�Ӧ�����������εļнǺ����
//			total_area += area;
//			tempMat *= ang;
//			tangMat += tempMat;
//		}
//		tangMat /= (0.5 * total_area);//����������֮�͡�����������Ρ�
//		calculateCurvature(tangMat, Kmin, Kmax);
//
//		m_curvatureKmin.push_back(Kmin);
//		m_curvatureKmax.push_back(Kmax);
//	}
//	//visualCoordiate(m_curvatureKmin);//���ӻ���С����
//}

//void WaterMark::calEdgeVec(OpenMesh::HalfedgeHandle he, Matrix33 & estiMat)
//{
//	Vector3 eVector;
//	eVector.setZero();
//	OpenMesh::Vec3d bPt, ePt;
//	double lenVec = 0.0;
//	auto fH = myMesh->from_vertex_handle(he);
//	auto tH = myMesh->to_vertex_handle(he);
//	bPt = myMesh->point(fH);//���
//	ePt = myMesh->point(tH);//�յ�
//	eVector[0] = ePt[0] - bPt[0]; eVector[1] = ePt[1] - bPt[1]; eVector[2] = ePt[2] - bPt[2];
//	lenVec = eVector.norm();//��߳�
//	eVector.normalize();
//	estiMat = eVector * eVector.transpose();
//	estiMat *= lenVec;
//}

//void WaterMark::calAngleTris(OpenMesh::HalfedgeHandle he, double & _angle, double &_area)
//{
//	double area_f = 0.0;
//	double area_of = 0.0;
//	auto ohe = myMesh->opposite_halfedge_handle(he);//��öԱߵİ��handle
//													//Ҫ�����жϣ�����б߽������
//	auto face_he = myMesh->face_handle(he);//��߶�Ӧ����handle
//	auto face_ohe = myMesh->face_handle(ohe);
//	if (myMesh->is_valid_handle(face_he) && myMesh->is_valid_handle(face_ohe))//�򵥵�ʡȥ����
//	{
//		OpenMesh::Vec3d fheNorm = myMesh->normal(face_he);
//		fheNorm.normalize();
//		OpenMesh::Vec3d foheNorm = myMesh->normal(face_ohe);
//		foheNorm.normalize();
//		double inP = (fheNorm[0] * foheNorm[0]) + (fheNorm[1] * foheNorm[1]) + (fheNorm[2] * foheNorm[2]);
//		_angle = acos(inP);//�淨��н�
//		area_f = calTriArea(face_he);
//		area_of = calTriArea(face_ohe);
//		_area = area_f + area_of;//�����������
//	}
//}

//double WaterMark::calTriArea(OpenMesh::FaceHandle faceH)
//{
//	double areaFace = 0.0;
//	RowVector3 point[3];
//	OpenMesh::Vec3d po;
//	int i = 0;
//	for (auto f_vit = myMesh->fv_begin(faceH); f_vit != myMesh->fv_end(faceH); ++f_vit, ++i)
//	{
//		po = myMesh->point(f_vit.handle());
//		point[i] << po[0], po[1], po[2];
//	}
//	point[1] -= point[0];
//	point[2] -= point[0];
//	point[0] = point[1].cross(point[2]);
//	areaFace = abs(0.5 * point[0].norm());
//	return areaFace;
//}

PolygonMesh::Mesh * WaterMark::my_embed_wm(PolygonMesh::Mesh * _mesh)
{
	std::cout << "Enter my_embed_wm() Method." << std::endl;
	
	EigenDeformation eigenDef;
	eigenDef.Init(_mesh);//��ʼ��������ز���

						 /////*
						 ////��Ҫ�ָ�ʱ��ȥ��RenderPGMeshEntity.cpp340�к�362�е�ע��
						 ////*/
						 ////segmentaionMesh(_mesh);

	eigenDef.calLap_Matrix(); // @firejq: ���������������˹���󲢵���matlab��������ֵ�ֽ�
	eigenDef.normVec();//������������λ��
	eigenDef.calR_Matrix();//@firejq: ����Ƶ��ϵ������
	eigenDef.embedWatermark();//way1
							  //eigenDef.embedByL();//way2

	std::cout << "Leave my_embed_wm() Method." << std::endl;
	return _mesh;
}

void WaterMark::my_extract_wm(PolygonMesh::Mesh * _mesh)
{
	EigenDeformation eigenDef;
	eigenDef.Init(_mesh);//��ʼ��������ز���
						 //���ļ��ж�ȡE����ʱ���������¼���
						 //eigenDef.calLap_Matrix();
						 //eigenDef.normVec();//������������λ��
	eigenDef.extractWatermark();//way1
								//eigenDef.extractByL();//way2
}

//void WaterMark::segmentationMesh(PolygonMesh::Mesh * _mesh)
//{
//	num_vtx = _mesh->n_vertices();
//	num_of_cluster = CLUSTER_NUM;//���÷ָ����
//	int dim = 6;//��Ϊ����������������ά��
//	MatrixXX	feature_vector(num_vtx, dim);
//
//	all_non_bdy = new double[num_of_cluster];//��fonud_Boundary��ʹ��
//
//											 //��������
//	VectorXd xCoordinate(num_vtx);
//	VectorXd yCoordinate(num_vtx);
//	VectorXd zCoordinate(num_vtx);
//
//	//vnx=vertex normal x axis,���㷨��
//	VectorXd Vnx(num_vtx);
//	VectorXd Vny(num_vtx);
//	VectorXd Vnz(num_vtx);
//
//	//���������긳ֵ��VectorXd
//	int j = 0;
//	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit, j++)
//	{
//		auto t_p = _mesh->point(vit.handle());
//
//		xCoordinate[j] = t_p[0];
//		yCoordinate[j] = t_p[1];
//		zCoordinate[j] = t_p[2];
//	}
//
//	/*��ȡ���ж���ķ���add by shiqun*/
//	_mesh->request_face_normals();
//	_mesh->request_vertex_normals();
//	_mesh->update_normals();
//	int vnidx = 0;
//	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit, vnidx++)
//	{
//		auto vertex = vit.handle();
//		OpenMesh::Vec3d v_normal;
//		v_normal = _mesh->normal(vertex);
//		Vnx[vnidx] = v_normal.data()[0];
//		Vny[vnidx] = v_normal.data()[1];
//		Vnz[vnidx] = v_normal.data()[2];
//	}
//	/*һ���ֲο�YQ��һ���ֲο�OpenMesh���ų������*/
//
//	feature_vector.col(0) = xCoordinate;
//	feature_vector.col(1) = yCoordinate;
//	feature_vector.col(2) = zCoordinate;
//	feature_vector.col(3) = Vnx;
//	feature_vector.col(4) = Vny;;
//	feature_vector.col(5) = Vnz;
//	//feature_vector.col(6) =v_curvatureKave;
//
//	Engine*	m_engine;
//	m_engine = NULL;//��ʼ��matlab����
//	if ((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab��laplacian��������ֵ
//	{
//		return;
//	}
//
//	int result_buffer_size = 1024 * 1000;
//	char* result_buffer = new char[result_buffer_size];
//	engOutputBuffer(m_engine, result_buffer, result_buffer_size);
//
//	engSetVisible(m_engine, 1);
//	engEvalString(m_engine, "cd('D:\\code\\GeometryProcessing-1\\mat_code')");//��matalab��������·��
//
//	mxArray*	mx_feature_vec = NULL;
//	double*	fv_buffer;
//	fv_buffer = new double[feature_vector.rows()*feature_vector.cols()];
//	int k = 0;
//	for (int i = 0; i<feature_vector.cols(); i++)
//	{
//		for (int j = 0; j<feature_vector.rows(); j++)
//		{
//			fv_buffer[k++] = feature_vector(j, i);
//		}
//	}
//
//	mx_feature_vec = mxCreateDoubleMatrix(feature_vector.rows(),
//		feature_vector.cols(), mxREAL);
//	memcpy((char*)mxGetPr(mx_feature_vec), (char*)fv_buffer,
//		feature_vector.rows()*feature_vector.cols() * sizeof(float));
//
//	engPutVariable(m_engine, "data", mx_feature_vec);
//	char	cmd_buf[128];
//	engEvalString(m_engine, cmd_buf);
//	///use default sigma
//	//sprintf(cmd_buf,"[labels]=sc(data,%d,400,0,%d);",num_of_neighbours, num_of_cluster);
//	///just kmeans here
//	sprintf(cmd_buf, "[labels]=k_means(data,\'random\',%d);", num_of_cluster);
//	engEvalString(m_engine, cmd_buf);
//
//	///Display output information
//	printf("%s", result_buffer);//��matlab�������ʱ���������������Ϣ
//
//	mxArray* mx_labels = NULL;
//	mx_labels = engGetVariable(m_engine, "labels");
//
//	double* m_labels = NULL;
//	m_labels = (double *)mxGetData(mx_labels);
//	for (int i = 0; i <num_vtx; i++)
//	{
//		all_labels.push_back(m_labels[i]);
//	}
//
//	if (mx_labels == NULL)
//	{
//		mxDestroyArray(mx_feature_vec);
//		mxDestroyArray(mx_labels);
//		engClose(m_engine);
//		delete[] fv_buffer;
//		delete[] result_buffer;
//		cout << "Get Data from Matlab error. End Clustering." << std::endl;
//		return;
//	}
//
//	double* labels = mxGetPr(mx_labels);
//	assert((mxGetNumberOfElements(mx_labels)) == num_vtx);
//
//	///** Step 3: Store the results **/
//	mxDestroyArray(mx_feature_vec);
//	mxDestroyArray(mx_labels);
//
//	engClose(m_engine);
//	delete[] fv_buffer;
//	delete[] result_buffer;
//
//	visualSegmentations(_mesh);
//
//	Init(_mesh);
//	saveSegmentation();
//
//	//find_Boundary(_mesh, 1);
//}

//void WaterMark::saveSegmentation()
//{
//	for (int i = 0; i < num_vtx; i++)
//	{
//		Area_map[all_labels[i]].push_back(i);
//	}
//
//	//д���ļ�
//	char Area_name[100];
//	for (int area_idx = 1; area_idx < Area_map.size() + 1; area_idx++)
//	{
//		sprintf(Area_name, "D:\\code\\GeometryProcessing-1\\Txt\\Seg\\Area%d.txt", area_idx);
//		ofstream Areafile;
//		Areafile.open(Area_name, ios_base::out);
//		if (Areafile)
//		{
//			for (int i = 0; i < Area_map[area_idx].size(); i++)
//			{
//				Areafile << Area_map[area_idx][i] << " ";
//			}
//		}
//		Areafile.close();
//	}
//}

//void WaterMark::calculateCurvature(const Matrix33 & tanMat, double & Kmin, double & Kmax)
//{
//	Eigen::Vector3cd eigvals;
//	Vector3 eigvalReal;
//	eigvals = tanMat.eigenvalues();
//	eigvalReal[0] = eigvals[0].real(); eigvalReal[1] = eigvals[1].real(); eigvalReal[2] = eigvals[2].real();
//	Kmax = max(max(eigvalReal[0], eigvalReal[1]), eigvalReal[2]);
//	double temp = max(eigvalReal[0], eigvalReal[1]);
//	if (temp <= eigvalReal[2])
//	{
//		Kmin = temp;
//	}
//	else
//	{
//		Kmin = min(eigvalReal[0], eigvalReal[1]);
//	}
//}

//void WaterMark::visualCoordiate(vector<double> & diffCoor)
//{
//	double dif = 0.0;
//	int iter = 0;
//	double _min = 0.0; //7.63377e-006;
//	double _max = 0.0134309;
//	double dfc = 0.0;
//	double _mean = 0.0;
//	guassNormaize(diffCoor, _mean);//��������һ����˹normalization��
//	for (int i = 0; i < num_vtx; i++)
//	{
//		dfc = diffCoor[i];
//		_min = min(dfc, _min);
//		_max = max(dfc, _max);
//	}
//	//cout<<"min="<<_min<<endl;
//	//cout<<"max="<<_max<<endl;
//	//��������һ��
//	double dist = _max - _min;
//	for (int i = 0; i < num_vtx; i++)
//	{
//		diffCoor[i] = (diffCoor[i] - _min) / dist;
//	}
//	//�������С���ƽ��ģ
//	for (auto v_it = myMesh->vertices_begin(); v_it != myMesh->vertices_end(); ++v_it)
//	{
//		dif = diffCoor[iter];
//		Mesh::Color Cl;
//		double norV = /*(dif - _min) / (_max - _min+ 1e-5f)*/ dif * 100.0f;
//		double clampV = norV < 99.0f ? norV : 99.0f;
//		clampV = clampV > 0.0f ? clampV : 0.0f;
//		int colorPercent = (int)clampV;
//		unsigned char* colorPtr = colorMapSimple + colorPercent * 3;
//
//		Cl.values_[0] = colorPtr[0];
//		Cl.values_[1] = colorPtr[1];
//		Cl.values_[2] = colorPtr[2];
//		myMesh->set_color(v_it.handle(), Cl);
//		iter++;
//	}
//	((PolygonMesh::PGMeshEntity *)myMesh)->update_rendering();
//}

//void WaterMark::guassNormaize(vector<double> & _cur, double & mean_)
//{
//	if (_cur.empty())
//	{
//		cout << "����Ϊ��" << endl;
//		return;
//	}
//	int size_ = _cur.size();
//
//	//������ŷ�Ͼ��������Χ�жԽ��ߵı�ֵ
//	PGMeshEntity * meshE = new PGMeshEntity;
//	meshE->set_mesh(myMesh);
//	meshE->cal_bounding_box();
//	double dig = meshE->get_bb_radius() * 2;
//	assert(dig >0.000001);
//
//	for (int i = 0; i < size_; i++)
//	{
//		_cur[i] /= dig;
//	}
//
//}

//void WaterMark::visualSegmentations(PolygonMesh::Mesh * _mesh)
//{
//	int setIndex = 0;
//	int counter = 0;
//	Mesh::Color Cl;
//	for (auto v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it, counter++)
//	{
//		setIndex = all_labels[counter];//��ȡ����ķָ�����
//									   //setIndex *= 10;//��������
//									   //int coIndex = setIndex > 9 ? setIndex % 9 :setIndex ;
//									   //unsigned char * colorPtr = colorMapSimple + coIndex;//�����ɫָ��
//		unsigned char * colorPtr = colorMapSimple + 3 * (setIndex - 1);//���޸ĺ�Ĵ���
//		Cl.values_[0] = colorPtr[0];
//		Cl.values_[1] = colorPtr[1];
//		Cl.values_[2] = colorPtr[2];
//		_mesh->set_color(v_it.handle(), Cl);
//	}
//	//((PolygonMesh::PGMeshEntity *)_mesh)->update_rendering();
//}

/*ȷ�����Ϊarea�������еı߽綥��
����area�зǱ߽�����
*/
//void WaterMark::find_Boundary(PolygonMesh::Mesh * _mesh, int area)
//{
//	int v_idx = 0;
//	cnt_non_bdy = 0;//�Ǳ߽綥�����
//	for (auto v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it, v_idx++)
//	{
//		if (all_labels[v_idx] == area)
//		{
//			//��Openmesh�У�������Ҳ�Ǵ�0��ʼ��
//			for (auto vv_it = _mesh->vv_begin(v_it); vv_it != _mesh->vv_end(v_it); ++vv_it)
//			{
//				OpenMesh::VertexHandle vertex2 = vv_it.handle();
//				int v2_idx = vertex2.idx();
//
//				//��ĳ�������1��������һ�������ź�����ͬʱ���������Ϊ�߽綥��
//				if (all_labels[v2_idx] != area)
//				{
//					is_boundray[v_idx] = true;
//					break;
//				}
//			}
//		}
//		//����1�еķǱ߽綥��
//		if (all_labels[v_idx] == area && !is_boundray[v_idx])
//		{
//			cnt_non_bdy++;
//		}
//	}
//	all_non_bdy[area - 1] = cnt_non_bdy;
//}

//void WaterMark::cal_area_Lap(PolygonMesh::Mesh * _mesh, int area)
//{
//	/*������Ϊarea�������У��Ǳ߽綥���γɵ�������˹����*/
//	MatrixXX D_matrix;
//	MatrixXX A_matrix;
//	MatrixXX Lap_matrix;
//	D_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//	A_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//	Lap_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//
//	map_vec.clear();
//	for (int j = 0; j < _mesh->n_vertices(); j++)
//	{
//		if (all_labels[j] == area && !is_boundray[j])
//		{
//			//��¼ԭ�����е�������area�еĶ�����
//			map_vec.push_back(j);
//		}
//	}
//
//	//��map_vec�еļ�ֵ��ֵ�ߵ�
//	map<int, int> map_vertex;
//	map_vertex.clear();
//	for (int i = 0; i< map_vec.size(); i++)
//	{
//		map_vertex.insert(pair<int, int>(map_vec[i], i));
//	}
//
//
//	int v_idx = 0;
//	int n_vidx = 0;//��ͬһ�����еĶ������±��
//	VectorX byd_neibor;
//	byd_neibor.setZero(cnt_non_bdy);
//	for (auto v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it, v_idx++)
//	{
//		if (all_labels[v_idx] == area && !is_boundray[v_idx])
//		{
//			//������߽���ھӸ���
//			for (auto vv_it = _mesh->vv_begin(v_it); vv_it != _mesh->vv_end(v_it); ++vv_it)
//			{
//				OpenMesh::VertexHandle vertex2 = vv_it.handle();
//				int v2_idx = vertex2.idx();
//
//				//�ж��ھӶ����Ƿ�Ϊ�߽綥��,���ھӵĶ���������Ϊv2_idx,��is_boundray�еı��Ҫ��ǰһλ
//				if (is_boundray[v2_idx])
//				{
//					byd_neibor[n_vidx]++;
//				}
//			}
//
//			//��ö���i�Ķ�,Ҫ��ȥ�߽��
//			int val = _mesh->valence(v_it.handle());
//			D_matrix(n_vidx, n_vidx) = val - byd_neibor[n_vidx];
//
//			//��ö����1����
//			OpenMesh::VertexHandle vertex1 = v_it.handle();
//			int v1_idx = vertex1.idx();
//			for (auto vv_it = _mesh->vv_begin(v_it); vv_it != _mesh->vv_end(v_it); ++vv_it)
//			{
//				OpenMesh::VertexHandle vertex2 = vv_it.handle();
//				int v2_idx = vertex2.idx();
//				if (!is_boundray[v2_idx])
//				{
//					A_matrix(n_vidx, map_vertex[v2_idx]) = 1;//����Ҫ��vertex1��vertex2���±��
//				}
//			}
//
//			n_vidx++;//��ͬһ�����еĶ�������µı��
//		}
//	}
//
//	Lap_matrix = D_matrix - A_matrix;
//
//	/*��������˹�����������ֵ�ֽ�*/
//	/*Step1. ��������˹����ֵ��һά����*/
//	int l_index = 0;
//	double *PL = new double[cnt_non_bdy*cnt_non_bdy];
//	for (int l_col = 0; l_col < cnt_non_bdy; l_col++)
//	{
//		for (int l_row = 0; l_row < cnt_non_bdy; l_row++)
//		{
//			PL[l_index++] = Lap_matrix(l_row, l_col);
//		}
//	}
//
//	/*Step2. ��matlab���棬�������*/
//	Engine*	m_engine;
//	m_engine = NULL;//��ʼ��matlab����
//	if ((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab��laplacian��������ֵ
//	{
//		return;
//	}
//	engSetVisible(m_engine, 1);
//
//	mxArray *pL = mxCreateDoubleMatrix(cnt_non_bdy, cnt_non_bdy, mxREAL);
//	memcpy((void *)mxGetPr(pL), (void *)PL, cnt_non_bdy*cnt_non_bdy * sizeof(double));
//
//	engPutVariable(m_engine, "L", pL);
//	engEvalString(m_engine, "cd('D:\\code\\GeometryProcessing-1\\mat_code')");//����matalab����
//
//	char buffer[255];
//	buffer[254] = '\0';
//
//	engOutputBuffer(m_engine, buffer, 255);
//	engEvalString(m_engine, "[eigVector,eigValue] = calcLaplacian(L);");
//
//	printf("%s", buffer);//��matlab�������ʱ���������������Ϣ
//
//	mxArray * eigenVectorObj = NULL;
//	eigenVectorObj = engGetVariable(m_engine, "eigVector");
//	double  * eigenVector = NULL;
//	eigenVector = (double*)mxGetData(eigenVectorObj);
//
//	RowVectorX eVector_iter;
//	eVector_iter.resize(cnt_non_bdy);
//
//	m_eigenVector.clear();
//	for (int eig_index = 0; eig_index < cnt_non_bdy; eig_index++)//��eigenvector����m_eigenVector
//	{
//		eVector_iter.setZero();
//		for (int ver_pos = 0; ver_pos < cnt_non_bdy; ver_pos++)
//		{
//			eVector_iter[ver_pos] = eigenVector[eig_index * cnt_non_bdy + ver_pos];
//		}
//		m_eigenVector.push_back(eVector_iter);
//	}
//
//	mxDestroyArray(pL);
//	mxDestroyArray(eigenVectorObj);
//}
//
//void WaterMark::cal_area_Rmatrix(PolygonMesh::Mesh * _mesh, int area)
//{
//	//��������
//	VectorXd xCoordinate(cnt_non_bdy);
//	VectorXd yCoordinate(cnt_non_bdy);
//	VectorXd zCoordinate(cnt_non_bdy);
//
//
//	//���������긳ֵ��VectorXd
//	int j = 0;
//	int k = 0;
//	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit, j++)
//	{
//		if (all_labels[j] == area && !is_boundray[j])
//		{
//			auto t_p = _mesh->point(vit.handle());
//
//			xCoordinate[k] = t_p[0];
//			yCoordinate[k] = t_p[1];
//			zCoordinate[k] = t_p[2];
//
//			////��¼ԭ�����е�������1�еĶ�����
//			//map_vec.push_back(j);
//			k++;
//		}
//	}
//
//	MatrixX3 V_matrix;
//	V_matrix.setZero(cnt_non_bdy, 3);
//
//	V_matrix.col(0) = xCoordinate;//V_matrix�Ƕ����������nx3
//	V_matrix.col(1) = yCoordinate;
//	V_matrix.col(2) = zCoordinate;
//
//	int vetex_index = 0;
//	double *vetex_coor = new double[cnt_non_bdy * 3];
//	for (int col_index = 0; col_index < 3; col_index++)
//	{
//		for (int row_index = 0; row_index < cnt_non_bdy; row_index++)
//		{
//			vetex_coor[vetex_index++] = V_matrix(row_index, col_index);
//		}
//	}
//
//	E_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//	for (int e_index = 0; e_index< cnt_non_bdy; e_index++)
//	{
//		E_matrix.col(e_index) = m_eigenVector[e_index];//E_matrix��������������nxn
//	}
//
//	int eig_index = 0;
//	double *eigVec = new double[cnt_non_bdy*cnt_non_bdy];
//	for (int axi_col = 0; axi_col < cnt_non_bdy; axi_col++)
//	{
//		for (int axi_row = 0; axi_row < cnt_non_bdy; axi_row++)
//		{
//			eigVec[eig_index++] = E_matrix(axi_row, axi_col);
//		}
//	}
//
//	Engine*	m_engine;
//	m_engine = NULL;//��ʼ��matlab����
//	if ((!m_engine && !(m_engine = engOpen(NULL))))// ʹ��matlab��laplacian��������ֵ
//	{
//		return;
//	}
//	engSetVisible(m_engine, 1);
//
//	mxArray *Vetex_Coor = mxCreateDoubleMatrix(cnt_non_bdy, 3, mxREAL);//�����񶥵����������matlab
//	mxArray *Eigen_vetor = mxCreateDoubleMatrix(cnt_non_bdy, cnt_non_bdy, mxREAL);
//
//	memcpy((void *)mxGetPr(Eigen_vetor), (void *)eigVec, (cnt_non_bdy * cnt_non_bdy) * sizeof(double));
//	memcpy((void *)mxGetPr(Vetex_Coor), (void *)vetex_coor, (cnt_non_bdy * 3) * sizeof(double));
//
//	engPutVariable(m_engine, "A", Eigen_vetor);//���󷽳�ΪAX = B,����A = E_matrix�� X = R_matrix�� B = V_matrix
//	engPutVariable(m_engine, "B", Vetex_Coor);
//
//	//buffer�������յ�����Ϣ����matlab�����д�ʱ���������buffer�鿴������Ϣ
//	char buffer[255];
//	buffer[254] = '\0';
//	engOutputBuffer(m_engine, buffer, 255);
//
//	engEvalString(m_engine, "cd('D:\\code\\GeometryProcessing-1\\mat_code')");//shenzhen
//	engEvalString(m_engine, "R = calcR(A, B);");
//
//	printf("%s", buffer);//��matlab�������ʱ���������������Ϣ
//
//	mxArray * r_matrix = NULL;
//	r_matrix = engGetVariable(m_engine, "R");
//	double  * Rmatrix = NULL;
//	Rmatrix = (double*)mxGetData(r_matrix);
//
//	R_matrix.setZero(cnt_non_bdy, 3);
//	int r_index = 0;
//	for (int r_col = 0; r_col < 3; r_col++)
//	{
//		for (int r_row = 0; r_row < cnt_non_bdy; r_row++)
//		{
//			R_matrix(r_row, r_col) = Rmatrix[r_index++];
//		}
//	}
//
//	mxDestroyArray(Vetex_Coor);
//	mxDestroyArray(Eigen_vetor);
//	mxDestroyArray(r_matrix);
//	// ����ʱ���ùر�
//	/*if (m_engine)
//	{
//	engClose(m_engine);
//	m_engine = NULL;
//	}*/
//
//}
//
//void WaterMark::embed_by_area(PolygonMesh::Mesh * _mesh)
//{
//#ifdef AREA_WATERMARK
//
//	//����ʱ��
//	time_t start = 0, end = 0;
//	time(&start);
//
//	segmentationMesh(_mesh);
//
//	num_vtx = _mesh->n_vertices();
//	is_boundray.clear();
//	for (int i = 0; i < num_vtx; i++)
//	{
//		//Ĭ�����ж��㶼���Ǳ߽統ĳ�������1��������һ�������ź�����ͬʱ���������Ϊ�߽綥��
//		is_boundray.push_back(false);
//	}
//
//	for (int area_num = 1; area_num < num_of_cluster + 1; area_num++)
//	{
//		find_Boundary(_mesh, area_num);
//	}
//
//	sort(all_non_bdy, all_non_bdy + num_of_cluster);
//	//for (int i = 0; i<num_of_cluster; i++)
//	//{
//	//	cout<<all_non_bdy[i]<<endl;
//	//}
//	//
//	chip_rate = 5;
//	setM((int)ceil((double)all_non_bdy[0] / chip_rate));
//	cout << "M is: " << getM() << endl;
//	//setM( (int)ceil((double)cnt_non_bdy/chip_rate) );//����ԭʼˮӡλ��
//	setC(chip_rate);//������Ƭ����
//	setAlpha(0.005);
//	setKey(7);//���ò���ˮӡ�����������
//
//	createA();
//	createWB();
//	createP();
//
//	char Efile_name[100];
//	for (int area_num = 1; area_num < num_of_cluster + 1; area_num++)
//	{
//		find_Boundary(_mesh, area_num);
//		cal_area_Lap(_mesh, area_num);
//		cal_area_Rmatrix(_mesh, area_num);
//
//		//��matlab�����У��ѽ�����ֵ��С�������У���С������ֵ��Ӧ����������Ϊ��Ƶ
//		//ÿ������ˮӡλ����ͬ������С����Ķ���������
//		for (int i = 0; i < getM()*chip_rate - 5; i++)
//		{
//			E_matrix(i, i) = E_matrix(i, i) + vecB[i] * P[i] * alpha;
//			/*E_matrix(i,i+1) = E_matrix(i,i+1) + vecB[i] * P[i] * alpha;
//			E_matrix(i+1,i) = E_matrix(i+1,i) + vecB[i] * P[i] * alpha;*/
//			//E_matrix(i-1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
//			//E_matrix(i+1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
//		}
//
//		for (int e_index = 0; e_index< cnt_non_bdy; e_index++)
//		{
//			m_eigenVector[e_index] = E_matrix.col(e_index);//E_matrix��������������nxn
//		}
//
//		//���޸ĺ�ĵ�λ������������д���ļ���
//		ofstream WEfile;
//		sprintf(Efile_name, "D:\\code\\GeometryProcessing-1\\Txt\\A%dWNe.txt", area_num);
//		WEfile.open(Efile_name, ios_base::out);
//		if (WEfile)
//		{
//			for (int i = 0; i < cnt_non_bdy; i++)
//			{
//				for (int j = 0; j < cnt_non_bdy; j++)
//				{
//					WEfile << E_matrix(i, j) << " ";
//				}
//			}
//		}
//		WEfile.close();
//
//		//Ƶ��ϵ������
//		VectorXd Rs(cnt_non_bdy);
//		VectorXd Rt(cnt_non_bdy);
//		VectorXd Ru(cnt_non_bdy);
//
//		Rs = R_matrix.col(0);
//		Rt = R_matrix.col(1);
//		Ru = R_matrix.col(2);
//
//		//���任��Ķ�������
//		VectorXd newxCoordinate(cnt_non_bdy);
//		VectorXd newyCoordinate(cnt_non_bdy);
//		VectorXd newzCoordinate(cnt_non_bdy);
//
//		newxCoordinate.setZero();
//		newyCoordinate.setZero();
//		newzCoordinate.setZero();
//
//		VectorXd  normV(cnt_non_bdy);//Eigen���͵����飬������ʱ����
//		for (int i_for_set = 0; i_for_set< cnt_non_bdy; i_for_set++)
//		{
//			normV = m_eigenVector[i_for_set];
//
//			newxCoordinate = newxCoordinate + Rs[i_for_set] * normV;//X = ��Rs[i] . ei
//			newyCoordinate = newyCoordinate + Rt[i_for_set] * normV;
//			newzCoordinate = newzCoordinate + Ru[i_for_set] * normV;
//		}
//
//		//���޸ĺ�Ķ������긳�ظ������ϵĶ���
//		int vi = 0;
//		int nvi = 0;
//		for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit, vi++)
//		{
//			if (nvi < cnt_non_bdy && map_vec[nvi] == vi)
//			{
//				auto& t_p = _mesh->point(vit.handle());
//
//				t_p[0] = newxCoordinate[nvi];
//				t_p[1] = newyCoordinate[nvi];
//				t_p[2] = newzCoordinate[nvi];
//
//				nvi++;
//			}
//		}
//	}
//
//	time(&end);
//	cout << "ˮӡǶ��ʱ��" << (end - start) << "��" << endl;
//
//	_mesh->update_face_normals();
//	//��������Ƕ�����
//
//#else
//	//ֻ����һ������
//	segmentationMesh(_mesh);
//
//	num_vtx = _mesh->n_vertices();
//	is_boundray.clear();
//	for (int i = 0; i < num_vtx; i++)
//	{
//		//Ĭ�����ж��㶼���Ǳ߽統ĳ�������1��������һ�������ź�����ͬʱ���������Ϊ�߽綥��
//		is_boundray.push_back(false);
//	}
//
//	//ȷ��ÿ������ı߽�
//	find_Boundary(_mesh, 1);
//	//����ÿ�������������˹���󲢽�������ֵ�ֽ�
//	cal_area_Lap(_mesh, 1);
//	cal_area_Rmatrix(_mesh, 1);
//
//	chip_rate = 7;
//	setM((int)ceil((double)cnt_non_bdy / chip_rate));//����ԭʼˮӡλ��
//	setC(chip_rate);//������Ƭ����
//	setAlpha(0.005);
//	setKey(7);
//
//	createA();
//	createWB();
//	createP();
//
//
//	//��matlab�����У��ѽ�����ֵ��С�������У���С������ֵ��Ӧ����������Ϊ��Ƶ
//	for (int i = 0; i < cnt_non_bdy - 5; i++)
//	{
//		E_matrix(i, i) = E_matrix(i, i) + vecB[i] * P[i] * alpha;
//		/*	E_matrix(i,i+1) = E_matrix(i,i+1) + vecB[i] * P[i] * alpha;
//		E_matrix(i+1,i) = E_matrix(i+1,i) + vecB[i] * P[i] * alpha;*/
//		//E_matrix(i-1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
//		//E_matrix(i+1,i) = E_matrix(i,i+1) + wm.vecB[i] * wm.P[i] * wm.alpha;
//	}
//
//	for (int e_index = 0; e_index< cnt_non_bdy; e_index++)
//	{
//		m_eigenVector[e_index] = E_matrix.col(e_index);//E_matrix��������������nxn
//	}
//
//	char Efile_name[100];
//	sprintf(Efile_name, "D:\\code\\GeometryProcessing-1\\Txt\\A%dWE.txt", 4);
//	//���޸ĺ�ĵ�λ������������д���ļ���
//	ofstream WEfile;
//	WEfile.open(Efile_name, ios_base::out);
//	if (WEfile)
//	{
//		for (int i = 0; i < cnt_non_bdy; i++)
//		{
//			for (int j = 0; j < cnt_non_bdy; j++)
//			{
//				WEfile << E_matrix(i, j) << " ";
//			}
//		}
//	}
//	WEfile.close();
//
//	//Ƶ��ϵ������
//	VectorXd Rs(cnt_non_bdy);
//	VectorXd Rt(cnt_non_bdy);
//	VectorXd Ru(cnt_non_bdy);
//
//	Rs = R_matrix.col(0);
//	Rt = R_matrix.col(1);
//	Ru = R_matrix.col(2);
//
//	//���任��Ķ�������
//	VectorXd newxCoordinate(cnt_non_bdy);
//	VectorXd newyCoordinate(cnt_non_bdy);
//	VectorXd newzCoordinate(cnt_non_bdy);
//
//	newxCoordinate.setZero();
//	newyCoordinate.setZero();
//	newzCoordinate.setZero();
//
//	VectorXd  normV(cnt_non_bdy);//Eigen���͵����飬������ʱ����
//	for (int i_for_set = 0; i_for_set< cnt_non_bdy; i_for_set++)
//	{
//		normV = m_eigenVector[i_for_set];
//
//		newxCoordinate = newxCoordinate + Rs[i_for_set] * normV;//X = ��Rs[i] . ei
//		newyCoordinate = newyCoordinate + Rt[i_for_set] * normV;
//		newzCoordinate = newzCoordinate + Ru[i_for_set] * normV;
//	}
//
//
//	//���޸ĺ�Ķ������긳�ظ������ϵĶ���
//	int vi = 0;
//	int nvi = 0;
//	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit, vi++)
//	{
//		if (nvi < cnt_non_bdy && map_vec[nvi] == vi)
//		{
//			auto& t_p = _mesh->point(vit.handle());
//
//			t_p[0] = newxCoordinate[nvi];
//			t_p[1] = newyCoordinate[nvi];
//			t_p[2] = newzCoordinate[nvi];
//
//			nvi++;
//		}
//	}
//
//	_mesh->update_face_normals();
//	std::cout << "Hello world  !" << std::endl;
//	//��������Ƕ�����
//#endif
//}
//
//void WaterMark::extract_by_area(PolygonMesh::Mesh * _mesh)
//{
//	Init(_mesh);
//	num_of_cluster = CLUSTER_NUM;
//
//	//����ʱ��
//	time_t start = 0, end = 0;
//	time(&start);
//
//	//��ȡ�ָ���
//	char Area_name[100];
//	vector<int> readA;
//	for (int area_idx = 1; area_idx < num_of_cluster + 1; area_idx++)
//	{
//		readA.clear();
//		sprintf(Area_name, "D:\\code\\GeometryProcessing-1\\Txt\\Seg\\Area%d.txt", area_idx);
//		ifstream AreaRfile;
//		AreaRfile.open(Area_name, ios_base::in);
//		if (AreaRfile)
//		{
//			// read into memory
//			AreaRfile.seekg(0, AreaRfile.end);
//			int length = AreaRfile.tellg();
//			AreaRfile.seekg(0, AreaRfile.beg);
//
//			char *buffer = new char[length];
//			AreaRfile.read(buffer, length);
//			AreaRfile.close();
//
//			// parse into array
//			std::istringstream iss(buffer);
//			int temp;
//			while (iss >> temp)
//			{
//				readA.push_back(temp);
//			}
//
//			readArea.push_back(readA);
//			delete[] buffer;
//			// print or use it.
//		}
//		AreaRfile.close();
//	}
//
//	//��all_labels��ʼ��
//	all_labels.clear();
//	for (int i = 0; i < num_vtx; i++)
//	{
//		all_labels.push_back(0);
//	}
//
//	//ȷ��ÿ����ķָ���
//	for (int i = 0; i < num_of_cluster; i++)
//	{
//		for (int j = 0; j < readArea[i].size(); j++)
//		{
//			all_labels[readArea[i][j]] = i + 1;
//		}
//	}
//
//	is_boundray.clear();
//	for (int i = 0; i < num_vtx; i++)
//	{
//		//Ĭ�����ж��㶼���Ǳ߽統ĳ�������1��������һ�������ź�����ͬʱ���������Ϊ�߽綥��
//		is_boundray.push_back(false);
//	}
//
//	all_non_bdy = new double[num_of_cluster];
//	for (int area_num = 1; area_num < num_of_cluster + 1; area_num++)
//	{
//		find_Boundary(_mesh, area_num);
//	}
//
//	sort(all_non_bdy, all_non_bdy + num_of_cluster);
//	for (int i = 0; i<num_of_cluster; i++)
//	{
//		cout << all_non_bdy[i] << endl;
//	}
//
//	chip_rate = 5;
//	setM((int)ceil((double)all_non_bdy[0] / chip_rate));
//	setC(chip_rate);
//	cout << "M is: " << getM() << endl;
//	setKey(7);
//	createP();
//	cout << "Test in extract: " << getM()*getC() - chip_rate << endl;
//
//	//����ԭʼ��ˮӡ���ж���
//	double *RB = new double[getM()*getC()];
//	ifstream Wbfile;
//	Wbfile.open("D:\\code\\GeometryProcessing-1\\Txt\\Wb.txt", ios_base::in);
//	if (Wbfile)
//	{
//		// read into memory
//		Wbfile.seekg(0, Wbfile.end);
//		int length = Wbfile.tellg();
//		Wbfile.seekg(0, Wbfile.beg);
//
//		char *buffer = new char[length];
//		Wbfile.read(buffer, length);
//		Wbfile.close();
//
//		// parse into array
//		std::istringstream iss(buffer);
//		int i = 0;
//		while (i<getM()*getC())
//		{
//			iss >> RB[i++];
//		}
//		delete[] buffer;
//		// print or use it.
//	}
//	Wbfile.close();
//
//	char Efile_name[100];
//	for (int area_num = 1; area_num < num_of_cluster + 1; area_num++)
//	{
//		find_Boundary(_mesh, area_num);
//		cal_area_Lap(_mesh, area_num);
//
//		E_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//		for (int e_index = 0; e_index< cnt_non_bdy; e_index++)
//		{
//			E_matrix.col(e_index) = m_eigenVector[e_index];//E_matrix��������������nxn
//		}
//
//		//��ȡ�޸ĺ����������
//		double *pE = new double[cnt_non_bdy * cnt_non_bdy];
//		ifstream RWEfile;
//		sprintf(Efile_name, "D:\\code\\GeometryProcessing-1\\Txt\\A%dWNe.txt", area_num);
//		RWEfile.open(Efile_name, ios_base::out);
//		if (RWEfile)
//		{
//			// read into memory
//			RWEfile.seekg(0, RWEfile.end);
//			int length = RWEfile.tellg();
//			RWEfile.seekg(0, RWEfile.beg);
//
//			char *buffer = new char[length];
//			RWEfile.read(buffer, length);
//			RWEfile.close();
//
//			// parse into array
//			std::istringstream iss(buffer);
//			int i = 0;
//			while (iss >> pE[i++]);
//			delete[] buffer;
//			// print or use it.
//		}
//		RWEfile.close();
//
//		MatrixXX	RE_matrix;
//		RE_matrix.setZero(cnt_non_bdy, cnt_non_bdy);
//		int ecount = 0;
//		for (int i = 0; i < cnt_non_bdy; i++)
//		{
//			for (int j = 0; j < cnt_non_bdy; j++)
//			{
//				RE_matrix(i, j) = pE[ecount++];
//			}
//		}
//
//		////��P��ȡ����
//		//double *P = new double[cnt_non_bdy];
//		//ifstream Pfile;
//		//Pfile.open("D:\\code\\GeometryProcessing-1\\Txt\\P.txt",ios_base::in );
//		//if (Pfile)
//		//{
//		//	// read into memory
//		//	Pfile.seekg (0, Pfile.end);
//		//	int length = Pfile.tellg();
//		//	Pfile.seekg (0, Pfile.beg);
//
//		//	char *buffer = new char[length];
//		//	Pfile.read(buffer, length);
//		//	Pfile.close();
//
//		//	// parse into array
//		//	std::istringstream iss(buffer);
//		//	int i = 0;
//		//	while ( i < cnt_non_bdy)
//		//	{
//		//		iss >> P[i++];
//		//	}
//		//	delete [] buffer;
//		//	// print or use it.
//		//}
//		//Pfile.close();
//
//		double * result = new double[getM()*getC() - chip_rate];
//		for (int ei = 0; ei<getM()*getC() - chip_rate; ei++)
//		{
//			result[ei] = chip_rate*(RE_matrix(ei, ei) - E_matrix(ei, ei))* P[ei];
//		}
//
//		double * wB = new double[getM()*getC() - chip_rate];
//		for (int i = 0; i < getM()*getC() - chip_rate; i++)
//		{
//			wB[i] = _copysign(1.0, result[i]);
//		}
//
//		//�Ƚ϶�ȡ����Wb����ȡ����wB
//		double Bcounter = 0.0;
//		for (int bi = 0; bi < getM()*getC() - chip_rate; bi++)
//		{
//			if (RB[bi] == wB[bi])
//			{
//				Bcounter++;
//			}
//		}
//
//		time(&end);
//		cout << "ˮӡ��ȡʱ��" << (end - start) << "��" << endl;
//
//		double corr = Bcounter / (getM()*getC() - chip_rate);
//		cout << "The result is: " << corr << endl;
//
//		delete[] pE;
//		delete[] result;
//		delete[] wB;
//	}
//
//	delete[] all_non_bdy;
//	delete[] RB;
//
//	cout << "This is the end of extract_by_area." << endl;
//}

void WaterMark::createA()
{
	srand(time(NULL));//���������ֵ
	for (int i = 0; i < m; i++)
	{
		vecA.push_back(rand() % 2);
	}
}

void WaterMark::createWB()
{
	/*for( int i = 0; i < vecA.size() *c; i ++)
	{
	for( int j = 0; j <  vecA.size(); j ++)
	{
	if( j*c <= i && i < (j+1)*c)
	vecB.push_back( vecA[j]);
	}
	}*/
	//����ֵͼ��ת����ԭʼ��ˮӡ���ж���
	double *B = new double[1188];
	ifstream Wfile;
	Wfile.open("D:\\firejq\\����\\Watermarking\\Txt\\AMesh\\flower.txt", ios_base::in);
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
		while (i<1188)
		{
			iss >> B[i++];
		}
		delete[] buffer;
		// print or use it.
	}
	Wfile.close();

	for (int i = 0; i < 1188/*vecA.size() * c*/; i++)
	{
		vecB.push_back(B[i]);
		if (vecB[i] == 0)
			vecB[i] = -1;
	}

	//��ˮӡ����д���ļ���
	ofstream Wbfile;
	Wbfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Wb.txt", ios_base::out);
	if (Wbfile)
	{
		for (int i = 0; i < 1188/*vecA.size() * c*/; i++)
		{
			Wbfile << vecB[i] << " ";
		}
	}
	Wbfile.close();
}

void WaterMark::createP()
{
	srand(key);//���������ֵ
	for (int i = 0; i < m*c; i++)
	{
		P.push_back(rand() % 2);
	}


	for (int i = 0; i < P.size(); i++)
	{
		if (P[i] == 0)
			P[i] = -1;
	}

	//��Pд���ļ���
	ofstream Pfile;
	Pfile.open("D:\\firejq\\����\\Watermarking\\Txt\\P.txt", ios_base::out);
	if (Pfile)
	{
		for (int i = 0; i < m*c; i++)
		{
			Pfile << P[i] << " ";
		}
	}
	Pfile.close();
}
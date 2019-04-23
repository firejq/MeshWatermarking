#include <iostream>

// −−−−−−−−−−−−−−−−−−−− OpenMesh
//#include "PGMeshTypes.h"
// --------------------

// −−−−−−−−−−−−−−−−−−−− Eigen
//#include <Eigen/Sparse>
//#include <Eigen/Dense>
// --------------------

//using namespace Eigen;
//using namespace GeometryProcess;

///////////////////////////////////
//#include <GeometryProcess/SceneInterface/SceneInterface.h>
//#include <GeometryProcess/PolygonMesh/PGMeshTypes.h>
//using namespace GeometryProcess::PolygonMesh;
///////////////////////////////////


#include "Watermark.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

//嵌入水印入口
void scene_embed_wm(string mesh_model, string wm_result)
{
	PolygonMesh::Mesh _mesh;
	//_mesh = ((PolygonMesh::PGMeshEntity *)active_entity_)->get_mesh();

	if (!OpenMesh::IO::read_mesh(_mesh, mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read successfully" << std::endl;
	std::cout << "vertice number: " << _mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;


	//BaseEntity * _meshEntity;//调用另一个init函数时，可以这样初始化
	//_meshEntity = active_entity_;

	WaterMark mytest;
	PolygonMesh::Mesh * wm_mesh = mytest.my_embed_wm(&_mesh);
	
	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(*wm_mesh, wm_result)) {
			std::cout << "Cannot write mesh to file \'" << wm_result << "\'" << std::endl;
			return;
		}
		else {
			std::cout << "Write mesh successfully!" << std::endl;
		}
	}
	catch (std::exception& x) {
		std::cout << x.what() << std::endl;
		return;
	}

	////修改完点坐标后，重新绘制模型
	//((PolygonMesh::PGMeshEntity *)active_entity_)->update_rendering();
}

//提取水印入口
void scene_extract_wm(string wm_mesh_model)
{
	//PolygonMesh::Mesh * _ori_mesh;
	//PolygonMesh::Mesh * _wm_mesh;
	//entity_list_ =entity_manager_->get_entities();
	//_ori_mesh =(PolygonMesh::Mesh *)  entity_list_[0];
	//_wm_mesh =(PolygonMesh::Mesh *)  entity_list_[1];

	//EigenDeformation eigenDef;
	//eigenDef. Init( _ori_mesh);//初始化网格相关参数
	PolygonMesh::Mesh _wm_mesh;
	//_wm_mesh = ((PolygonMesh::PGMeshEntity *)active_entity_)->get_mesh();

	if (!OpenMesh::IO::read_mesh(_wm_mesh, wm_mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read successfully" << std::endl;
	std::cout << "vertice number: " << _wm_mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _wm_mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;

	//BaseEntity * _meshEntity;//调用另一个init函数时，可以这样初始化
	//_meshEntity = active_entity_;

	WaterMark mytest;
	mytest.my_extract_wm(&_wm_mesh);

}


int main() {
	std::cout << "Hello world!" << std::endl;

	//scene_embed_wm("D:\\firejq\\repo\\MeshModels\\lowpolycow\\cow.obj", "result.obj");

	scene_extract_wm("result.obj");

	return 0;


	MyMesh mesh;
	/*MyMesh::VertexHandle vhandle[8];
	vhandle[0] = mesh.add_vertex(MyMesh::Point(-2, -2, 2));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(1, -1, 1));D:\MATLAB\matlabR2012a\extern\lib\win32\microsoft
	vhandle[2] = mesh.add_vertex(MyMesh::Point(1, 1, 1));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(-1, 1, 1));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(1, -1, -1));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(1, 1, -1));
	vhandle[7] = mesh.add_vertex(MyMesh::Point(-1, 1, -1));

	std::vector<MyMesh::VertexHandle> face_vhandles;
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[7]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[5]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[7]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[7]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "output.obj")) {
			std::cout << "Cannot write mesh to file 'output.off'" << std::endl;
			return 1;
		}
		else {
			std::cout << "Write mesh successfully!" << std::endl;
		}
	}
	catch (std::exception& x) {
		std::cout << x.what() << std::endl;
		return 1;
	}
	*/

	if (!OpenMesh::IO::read_mesh(mesh, "D:\\firejq\\毕设\\MeshModels\\lowpolydeer\\deer.obj")) {
		std::cout << "read error\n";
		return 1;
	}
	else {
		std::cout << "read successfully" << std::endl;

		std::cout << "vertice number: " << mesh.n_vertices() << std::endl;
		std::cout << "face number: " << mesh.n_faces() << std::endl;

	}







	return 0;

}
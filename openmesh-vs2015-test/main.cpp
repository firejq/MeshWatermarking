#include "MeshWatermark.h"

//嵌入水印入口
void scene_embed_wm(string mesh_model, string wm_result)
{
	PolygonMesh::Mesh _mesh;

	if (!OpenMesh::IO::read_mesh(_mesh, mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read successfully" << std::endl;
	std::cout << "vertice number: " << _mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;

	//WaterMark mytest;
	//PolygonMesh::Mesh * wm_mesh = mytest.my_embed_wm(&_mesh);

	MeshWatermark meshWm;
	meshWm.Init(&_mesh);//初始化网格相关参数
	meshWm.embedWatermark();//way1
							  //meshWm.embedByL();//way2

	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(_mesh, wm_result)) {
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

}

//提取水印入口
void scene_extract_wm(string wm_mesh_model)
{
	PolygonMesh::Mesh _wm_mesh;

	if (!OpenMesh::IO::read_mesh(_wm_mesh, wm_mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read successfully" << std::endl;
	std::cout << "vertice number: " << _wm_mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _wm_mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;

	//WaterMark mytest;
	//mytest.my_extract_wm(&_wm_mesh);

	MeshWatermark meshWm;
	meshWm.Init(&_wm_mesh);//初始化网格相关参数
						 //从文件中读取E矩阵时，不需重新计算
						 //meshWm.calLap_Matrix();
						 //meshWm.normVec();//将特征向量单位化
	meshWm.extractWatermark();//way1
								//meshWm.extractByL();//way2

}


int main() {
	std::cout << "Hello world!" << std::endl;

	//string cur_path = _getcwd(NULL, 0);
	//std::cout << cur_path << std::endl;
	//return 0;

	string ori_mesh = "D:\\firejq\\repo\\openmesh-vs2015-test\\mesh_models\\original\\lowpolycow\\cow.obj";
	string wm_mesh = "D:\\firejq\\repo\\openmesh-vs2015-test\\mesh_models\\watermarked\\cowresult.obj";

	//scene_embed_wm(ori_mesh, wm_mesh);

	scene_extract_wm(wm_mesh);

	return 0;

}
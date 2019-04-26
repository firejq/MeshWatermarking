#include "MeshWatermark.h"

//嵌入水印入口
void scene_embed_wm(string mesh_model, string watermark_path, string wm_result, int c = 7, double a = 0.005)
{
	std::cout << "-----------------------------" << std::endl;

	PolygonMesh::Mesh _mesh;

	if (!OpenMesh::IO::read_mesh(_mesh, mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read original mesh successfully" << std::endl;
	std::cout << "vertice number: " << _mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Embed process starts:" << std::endl;

	MeshWatermark meshWm;
	meshWm.Init(&_mesh, watermark_path, c, a);//初始化网格相关参数
	meshWm.embedWatermark();//嵌入水印

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

	std::cout << "Embed process finished successfully!" << std::endl;
	std::cout << "-----------------------------" << std::endl;

}

//提取水印入口
void scene_extract_wm(string wm_mesh_model, string ori_mesh_model, string ori_watermark, string extr_watermark)
{
	std::cout << "-----------------------------" << std::endl;

	PolygonMesh::Mesh _wm_mesh;

	if (!OpenMesh::IO::read_mesh(_wm_mesh, wm_mesh_model)) {
		std::cout << "read error\n";
		return;
	}

	std::cout << "read watermarked mesh successfully" << std::endl;
	std::cout << "vertice number: " << _wm_mesh.n_vertices() << std::endl;
	std::cout << "face number: " << _wm_mesh.n_faces() << std::endl;
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Extract process starts:" << std::endl;

	MeshWatermark meshWm;
	meshWm.Init(&_wm_mesh, ori_watermark);//初始化网格相关参数
						 //meshWm.calLap_Matrix();
						 //meshWm.normVec();//将特征向量单位化
	meshWm.extractWatermark(extr_watermark);//提取水印
	
	std::cout << "Extract process finished successfully!" << std::endl;
	std::cout << "-----------------------------" << std::endl;

}


int main() {
	std::cout << "Hello world!" << std::endl;
	
	/*
		obj: cow389/deer817/cat1136/wolf1368
		wm:  flower/heart
	*/
	string obj = "cat";
	string wm = "flower";

	string ori_mesh = ROOT_PATH + "file\\mesh_models\\original\\lowpoly" + obj + "\\" + obj + ".obj";
	string wm_mesh = ROOT_PATH + "file\\mesh_models\\watermarked\\" + obj + "-result.obj";
	string emb_watermark = ROOT_PATH + "file\\wmBinaryImage\\Embed\\" + wm + ".txt";
	string extr_watermark = ROOT_PATH + "file\\wmBinaryImage\\Extract\\extr_" + wm + "_from_" + obj + ".txt";
	const int chip_rate = 7;
	const double alpha = 0.005;
	//执行嵌入
	scene_embed_wm(ori_mesh, emb_watermark, wm_mesh, chip_rate, alpha);
	//执行提取
	scene_extract_wm(wm_mesh, ori_mesh, emb_watermark, extr_watermark);

	return 0;

}
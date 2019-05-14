#pragma once
#ifndef _INIT_H_
#define _INIT_H_  

/********c++ headers********/
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios_base;
/**************************/


/******OpenMesh Init*********/
#include <OpenMesh/Core/io/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
namespace PolygonMesh
{
	/** this structure defines the traits of the mesh
	*/
	struct CustomTraits : OpenMesh::DefaultTraits
	{
		// let Point and Normal be a vector made from doubles
		typedef OpenMesh::Vec3d Point;
		typedef OpenMesh::Vec3d Normal;

		// add normal property to vertices and faces
		VertexAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
		FaceAttributes(OpenMesh::Attributes::Normal);
		HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
	};

	//typedef OpenMesh::TriMesh_ArrayKernelT<CustomTraits>  TriangularMesh;
	typedef OpenMesh::PolyMesh_ArrayKernelT<CustomTraits> Mesh;

}
/**************************/


/********Eigen Init**********/
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3>              MatrixX3;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectorX;
/**************************/


/*********matlab**********/
#include "engine.h"
/**************************/


/******constant configuration*********/
const string ROOT_PATH = "D:\\firejq\\repo\\openmesh-vs2015-test\\";
/**************************/


#endif
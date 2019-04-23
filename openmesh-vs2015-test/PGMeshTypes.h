#pragma once
#ifndef _GEOMETRY_MESH_TYPES_
#define _GEOMETRY_MESH_TYPES_  

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
#endif

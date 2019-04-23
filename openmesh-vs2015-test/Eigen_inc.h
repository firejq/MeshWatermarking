#pragma once
#ifndef types_h
#define types_h

//#define EIGEN_USE_MKL_ALL

//#define EIGEN_DONT_VECTORIZE
//#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

////#include <Eigen/Sparse>
//// #include <Eigen/SparseLU>
//// #include <Eigen/SparseCholesky>
//#include <Eigen/Dense> 
//#include <Eigen/LU>
//#include <Eigen/Core>
//#include <Eigen/Geometry>
//#include <Eigen/Sparse>
//#include <Eigen/Cholesky>
//#include <Eigen/SVD>
////#include <GeometryProcess/EigenDeformation/Vector2.h>
////#include <GeometryProcess/EigenDeformation/Vector3.h>
//#include <Eigen/Geometry> 

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace GeometryProcess
{
#define USE_DOUBLE

	typedef int IndexType;

#ifdef USE_DOUBLE
	typedef double ScalarType;
#define TW_TYPE_SCALARTYPE TW_TYPE_DOUBLE
#else
	typedef float ScalarType;
#define TW_TYPE_SCALARTYPE TW_TYPE_FLOAT
#endif
	//
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 3, Eigen::RowMajor> PointMatrixType;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 3, Eigen::RowMajor> FaceMatrixType;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 4, Eigen::RowMajor> QuadMatrixType;
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 2, Eigen::RowMajor> UVMatrixType;
	//
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 2>              MatrixX2;
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 3>              MatrixX3;
	typedef Eigen::Matrix<ScalarType, 2, 2>              Matrix22;
	typedef Eigen::Matrix<ScalarType, 3, 3>              Matrix33;
	typedef Eigen::Matrix<ScalarType, 4, 4>              Matrix44;
	typedef Eigen::Matrix<ScalarType, 7, 7>              Matrix77;
	typedef Eigen::Matrix<ScalarType, 6, 6>              Matrix66;
	typedef Eigen::Matrix<ScalarType, 10, 10>              Matrix1010;
	typedef Eigen::Matrix<ScalarType, 2, 3>              Matrix23;
	//
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorX;
	typedef Eigen::Matrix<ScalarType, 2, 1> Vector2;
	typedef Eigen::Matrix<ScalarType, 3, 1> Vector3;
	typedef Eigen::Matrix<ScalarType, 4, 1> Vector4;
	typedef Eigen::Matrix<ScalarType, 5, 1> Vector5;
	typedef Eigen::Matrix<ScalarType, 7, 1> Vector7;
	typedef Eigen::Matrix<ScalarType, 6, 1> Vector6;
	typedef Eigen::Matrix<ScalarType, 10, 1> Vector10;
	//
	typedef Eigen::Matrix<ScalarType, 1, Eigen::Dynamic> RowVectorX;
	typedef Eigen::Matrix<ScalarType, 1, 2>              RowVector2;
	typedef Eigen::Matrix<ScalarType, 1, 3>              RowVector3;
	typedef Eigen::Matrix<ScalarType, 1, 4>              RowVector4;
	typedef Eigen::Matrix<IndexType, 1, 2>              RowVector2i;
	//
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic> MatrixXXi;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 3>              MatrixX3i;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 2>              MatrixX2i;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 12>             MatrixX12i;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 8>              MatrixX8i;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 6>              MatrixX6i;
	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 4>              MatrixX4i;

	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 26>             MatrixX26i;
	//-----------
	typedef Eigen::Matrix<IndexType, 3, 6>              Matrix36i;

	typedef Eigen::Matrix<IndexType, Eigen::Dynamic, 1> VectorXi;
	typedef Eigen::Matrix<IndexType, 3, 1> Vector3i;
	typedef Eigen::Matrix<IndexType, 1, 3> RowVector3i;
	typedef Eigen::Matrix<IndexType, 1, Eigen::Dynamic> RowVectorXi;
	/*	typedef Eigen::Matrix<IndexType,  3,              3> Matrix33i;*/

	typedef Eigen::AngleAxis<ScalarType> AngleAxis;
	typedef Eigen::Quaternion<ScalarType> Quaternion;
	//
	typedef Eigen::SparseMatrix<ScalarType> SparseMatrixd;
	//typedef Eigen::Triplet<ScalarType,IndexType> SparseMatrixTriplet;

	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 3, Eigen::RowMajor> RowMatrixX3;
	typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> RowMatrixX1;
}


#endif


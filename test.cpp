#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "matrix.hpp"
#include <Eigen/Dense>
#include <time.h>
#include <stdio.h>

using namespace Eigen;

bool compareMatrixEigen(matrix<float> mat, MatrixXf eig, float eps = 1e-4)
{
	if (mat.rows() != eig.rows()) return false;
	if (mat.cols() != eig.cols()) return false;
	for (unsigned int i = 0; i < mat.rows(); ++i)
		for (unsigned int j = 0; j < mat.cols(); ++j)
			if (std::abs(mat(i,j) - eig(i,j)) > eps) return false;
	return true;
}

TEST_CASE( "constructors, properties, delete", "[matrix]") {

	SECTION( "default constructor (empty matrix), rows(), cols(), length(), data()") {

		// create empty matrix object	
		matrix<float> mA;

		// test metadata
		REQUIRE( mA.rows() == 0);
		REQUIRE( mA.cols() == 0);
		REQUIRE( mA.length() == 0);
		REQUIRE( mA.data() == 0);

	}
	SECTION( "zero-filled n x m matrix, element(), access-operators() and []") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;
	
		// create zero-filled matrix
		matrix<float> mA(nRows, nCols);

		// test metadata
		REQUIRE( mA.rows() == nRows);
		REQUIRE( mA.cols() == nCols);
		REQUIRE( mA.length() == nRows * nCols);
		REQUIRE( mA.data() != 0);
	
		// check each element 
		for (int i = 0; i < nRows * nCols; ++i){
			REQUIRE( mA.element(i) == 0);
			REQUIRE( mA(i) == 0);
			REQUIRE( *(mA[i]) == 0);
		}
	}
	

	SECTION( "initialize matrix with eigen matrix, check comparison function") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;
	
		// create Eigen Matrix, fill with random values
		Matrix<float, Dynamic, Dynamic, RowMajor> eA;
		eA.setRandom(nRows, nCols);

		// create matrix, fill with values of Eigen matrix
		matrix<float> mA(nRows, nCols, eA.data()); 

		// test metadata
		REQUIRE( mA.rows() == nRows);
		REQUIRE( mA.cols() == nCols);

		// compare matrices (Eigen and matrix)
		REQUIRE( compareMatrixEigen(mA, eA) == true);	
	
		// change one element to check comparison function
		mA(0) += 1;
		REQUIRE( compareMatrixEigen(mA, eA) == false);	
		
	

	}
	
	SECTION( "initialize eigen matrix with matrix") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;
	
		// create matrix, fill with random values
		matrix<float> mA(nRows, nCols);
		mA.setRandom();
		
		// create Eigen matrix, fill with matrix data
		Map<Matrix<float, Dynamic, Dynamic, RowMajor> > eA(mA.data(),nRows, nCols);
		
		// test metadata
		REQUIRE( mA.rows() == nRows);
		REQUIRE( mA.cols() == nCols);
	
		// compare matrices
		REQUIRE( compareMatrixEigen(mA, eA) == true);	
		

	}

	SECTION( "initialize matrix with data from csv file, save to csv, == operator") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;

		// set equality epsilon
		float eps = 1e-5;
	
		// create matrix, fill with random values
		matrix<float> mA(nRows, nCols);
		mA.setRandom();
	
		// save to file
		mA.save("./testFile.csv");
		
		// create second matrix, read data from testfile
		matrix<float> mB("./testFile.csv");
			
		// test metadata
		REQUIRE( mB.rows() == nRows);
		REQUIRE( mB.cols() == nCols);
	
		// test elements
		for (int i = 0; i < nRows * nCols; ++i){
			REQUIRE( std::abs(mA.element(i) - mB.element(i)) < eps);
			REQUIRE( std::abs(mA(i) - mB(i)) < eps);
			REQUIRE( std::abs(*(mA[i]) - *(mB[i])) < eps);

			// mA and mB should not use the same storage
			REQUIRE( mA[i] != mB[i]);
		}

		// compare matrices with == operator
		REQUIRE( mA == mB);

		// change one element to check == operator
		mA(0) += 1;

		// compare matrices with == operator again (should not be equal)
		REQUIRE( !( mA == mB));
		
		remove("./testFile.csv");

	}

	SECTION( "copy constructor") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;

		// create matrix, fill with random values
		matrix<float> mA(nRows, nCols);
		mA.setRandom();

		// create copy
		matrix<float> mB(mA);
	
		// compare matrices with == operator
		REQUIRE( mA == mB);
		
		// mA and mB should use different sorages
		REQUIRE( mA.data() != mB.data());

	}
	
	SECTION( "delete matrix (clear())") {

		// set dimensions
		unsigned int nRows = 6;
		unsigned int nCols = 9;
	
		// create matrix And fill with random values
		matrix<float> mA(nRows, nCols);
		mA.setRandom();

		// test metadata
		REQUIRE( mA.rows() == nRows);
		REQUIRE( mA.cols() == nCols);
		REQUIRE( mA.length() == nRows * nCols);
		REQUIRE( mA.data() != 0);
	
		// clear matrix
		mA.clear();
	
		// test metadata
		REQUIRE( mA.rows() == 0);
		REQUIRE( mA.cols() == 0);
		REQUIRE( mA.length() == 0);
		REQUIRE( mA.data() == 0);

	}



}

TEST_CASE( "matrix initializations", "[matrix]") {

	// set dimensions
	unsigned int nRows = 5;
	unsigned int nCols = 8;

	// set equality epsilon
	float eps = 1e-5;

	// create matrix, fill with random values between 1 and 2
	matrix<float> mA(nRows, nCols);
	mA.setRandom(1.,2.);

	// check elements
	for (int i = 0; i < nRows * nCols; ++i)
		REQUIRE( ((mA.element(i) >= 1.) & (mA.element(i) <= 2.)));

	SECTION( "setZero"){

		mA.setZero();
	
		// check elements
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( mA.element(i) == 0);
	}

	SECTION( "set"){

		mA.set(11.3);
	
		// check elements
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( std::abs(mA.element(i) - 11.3) <= eps);
	}

	SECTION( "set to matrix"){

		// create second matrix, fill with random values between 3 and 4
		matrix<float> mB(nRows, nCols);
		mB.setRandom(3.,4.);
	
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( ((mB.element(i) >= 3.) & (mB.element(i) <= 4.)));
	
		// fill matrix mA with values of mB
		mA.set(mB);

		REQUIRE(mA == mB);

		// if dimensions do not agree using set() matrix should not be changed
		matrix<float> mC(nRows, nCols+1);
		mC.setRandom();
		mA.set(mC);

		// mA still equal to mB
		REQUIRE(mA == mB);

	}

	SECTION( "setRandom"){

		// with no arguments, fill with random numbers between 0 and 1
		mA.setRandom();
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( ((mA.element(i) >= 0.) & (mA.element(i) <= 1.)));
		
		// two random matrices should not be equal
		matrix<float> mB(nRows, nCols);
		mB.setRandom();
		REQUIRE( !(mA == mB));
				
		// with two arguments, fill with random numbers between argument 1 and 2
		mA.setRandom(3., 12.);
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( ((mA.element(i) >= 3.) & (mA.element(i) <= 12.)));

		// incorrect order of min and max should not matter
		mA.setRandom(12., -3.);
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( ((mA.element(i) >= -3.) & (mA.element(i) <= 12.)));
	
		// same argument twice sould lead to matrix filled with that value without error
		mA.setRandom(-0.3, -0.3);
		for (int i = 0; i < nRows * nCols; ++i)
			REQUIRE( std::abs(mA.element(i) - -0.3) <= eps);
	}
}


TEST_CASE( "matrix construction", "[matrix]") {

	// create two random matrices A and B as base for calculations.
	// A and B are equal as matrix and EigenMatrix 

	unsigned int nRows = 7;
	unsigned int nCols = 11;
	clock_t start, end;


	// matrix A: matrix (mA) and Eigen (eA)
	matrix<float> mA(nRows, nCols);
	mA.setRandom();
	Map<Matrix<float, Dynamic, Dynamic, RowMajor> > eA(mA.data(),nRows, nCols);
	
	// matrix B: matrix (mB) and Eigen (eB)
	matrix<float> mB(nRows, nCols);
	mB.setRandom();
	Map<Matrix<float, Dynamic, Dynamic, RowMajor> > eB(mB.data(),nRows, nCols);

	REQUIRE( compareMatrixEigen(mA, eA) == true);	
	REQUIRE( compareMatrixEigen(mB, eB) == true);	

	REQUIRE( mA.rows() == nRows);
	REQUIRE( mA.cols() == nCols);
	REQUIRE( mB.rows() == nRows);
	REQUIRE( mB.cols() == nCols);

	SECTION( " == operator" ) {

		REQUIRE( mA == mA);
		REQUIRE( !(mB == mA)); 

	}


	SECTION( "basic operations, out of place." ) {
	
		matrix<float> mResult;
		MatrixXf eResult;

		mResult = mA + mB;
		eResult = eA + eB;
		REQUIRE( compareMatrixEigen(mResult, eResult) == true);	

		mResult = mA - mB;
		eResult = eA - eB;
		REQUIRE( compareMatrixEigen(mResult, eResult) == true);	

		mResult = mA * mB;
		eResult = eA.array() * eB.array();
		REQUIRE( compareMatrixEigen(mResult, eResult) == true);	

		mResult = mA / mB;
		eResult = eA.array() / eB.array();
		REQUIRE( compareMatrixEigen(mResult, eResult) == true);	

		mResult = matrix<float>::multiply(mA, mB, 0, 1);
		eResult = eA * eB.transpose();
		REQUIRE( compareMatrixEigen(mResult, eResult) == true);	

	}
}

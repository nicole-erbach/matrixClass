#define _USE_MATH_DEFINES
#include "matrix.hpp"
#include <cstring> // memcpy
#include <stdlib.h> // rand
#include <cmath> // isfinite
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "blas_api.h"

typedef unsigned int uint;

/** 
* @brief default constructor, construct matrix with no memory and zero size
*/

template<class T>
matrix<T>::matrix(){

	nRows = 0;
	nCols = 0;
	storage = 0;
	postResize();
}


/** 
* @brief construct zero-filled n x m matrix
*/

template<class T>
matrix<T>::matrix(uint n, uint m)
	: nRows(n)
	, nCols(m)
{
	
	storage = new T[nRows * nCols];
	setZero();
	postResize();
}


/** 
* @brief construct n x m matrix and fill with data at storageInit
*/

template<class T>
matrix<T>::matrix(uint n, uint m, T* storageInit)
	: nRows(n)
	, nCols(m)
{

	storage = new T[nRows * nCols];
	memcpy(storage, storageInit, nRows * nCols * sizeof(T));
	postResize();

}

/** 
* @brief construct zero-filled n-element vector (1 x n matrix)
*/

template<class T>
matrix<T>::matrix(uint n)
	: nRows(1)
	, nCols(n)
{

	storage = new T[nCols];
	setZero();
	postResize();
}


/** 
* @brief construct n-element vector (1 x n matrix) and fill with data at storageInit
*/

template<class T>
matrix<T>::matrix(uint n, T* storageInit)
	: nRows(1)
	, nCols(n)
{

	storage = new T[nCols];
	memcpy(storage, storageInit, nCols * sizeof(T));
	postResize();

}

/** 
* @brief load matrix from file
*/

template<class T>
matrix<T>::matrix(std::string pathFilename)
{
	storage = 0;
	readFromFile(pathFilename);
	postResize();

}

/** 
* @brief read file to get number of rows and cols and resize matrix
*/

template<class T>
void matrix<T>::prepareForFileread(std::string pathFilename)
{
	// reset matrix
	nCols = 0;
	nRows = 0;
	delete[] storage;
	
	// count number of rows and cols in file, count cols in each row to check consistency 
	std::string line, inLine;
	std::ifstream file;
	uint colsInLine = 0;	

	file.open(pathFilename.c_str());
	if (file.is_open())
	{
		// count cols in first line
		if (std::getline(file, line))
		{
			std::istringstream stream(line);
			while(std::getline(stream, inLine, ','))
				nCols++;
			nRows++;
		}
		else 
		{
			std::cout << pathFilename << " is empty" << std::endl;
			nRows = 0;
			nCols = 0;
			storage = 0;
			file.close();
			return;
		}	

		while (std::getline(file, line))
		{
			colsInLine = 0;
			nRows++;
			std::istringstream stream(line);
			while(std::getline(stream, inLine, ','))
				colsInLine++;
			if (colsInLine != nCols)
			{
				std::cout << colsInLine << ", " << nCols << std::endl;
				std::cout << "colums in line " << nRows << " do not agree with columns in line 1" << std::endl;
				file.close();
				nRows = 0;
				nCols = 0;
				storage = 0;
				return;
			}		
		}
	}
	else
	{
		//printf("could not open file %s\n", pathFilename.c_str());
		std::cout << "could not open file " << pathFilename << std::endl;//printf("could not open file %s\n", pathFilename.c_str());
		nRows = 0;
		nCols = 0;
		storage = 0;
		return;
	}

	file.close();
}

/** 
* @brief read file and store data in matrix
*/

template<class T>
void matrix<T>::readFromFile(std::string pathFilename)
{

	prepareForFileread(pathFilename);
	storage = new T[nCols * nRows];

	std::string line, inLine;
	std::ifstream file;
	
	uint n = 0;

	// read data from file
	file.open(pathFilename.c_str());
	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			std::istringstream stream(line);
			while(std::getline(stream, inLine, ','))
			{
				storage[n] = atof(inLine.c_str());
				n++;
			}
		}

	}
	else 
	{
		nCols = 0;
		nRows = 0;
		delete[] storage;
		storage = 0;
		std::cout << "problems reading file " << pathFilename << std::endl;
	}

	file.close();

}

template<class T>
matrix<T>::~matrix()
{

	delete[] storage;
	storage = 0;
}

/** 
* @brief copy constructor
*/

template<class T>
matrix<T>::matrix(const matrix<T>& matrixInit)
{
	
	nRows = matrixInit.rows();
	nCols = matrixInit.cols();
	storage = new T[nRows * nCols];
	memcpy(storage, matrixInit.data(),nRows * nCols * sizeof(T));
	postResize();
		
}

/** 
* @brief delete matrix
*/

template<class T>
void matrix<T>::clear()
{
	nRows = 0;
	nCols = 0;
	delete[] storage;
	storage = 0;
	postResize();
}

/** 
* @brief method to call after resizing to handle metadata. To be implemented in inherited classes. 
*/
template<class T>
void matrix<T>::postResize()
{
	return;
}

/** 
* @brief return number of rows
*/

template<class T>
uint matrix<T>::rows() const
{
	return(nRows);
}

/** 
* @brief return number of columns
*/

template<class T>
uint matrix<T>::cols() const
{
	return nCols;
}

/** 
* @brief return number of elements (rows * columns)
*/

template<class T>
uint matrix<T>::length() const
{
	return nCols * nRows;
}


/** 
* @brief return pointer to first element
*/

template<class T>
T* matrix<T>::data() const
{
	return storage;
}

/** 
* @brief return value of n-th element
*/

template<class T>
T matrix<T>::element(uint n) const
{
	return storage[n];
}

/** 
* @brief true if this is a vector
*/

template<class T>
bool matrix<T>::isVector() const
{
	if (nRows == 1 && nCols > 1) return true;
	if (nRows > 1 && nCols == 1) return true;
	return false;
}


/** 
* @brief set all values to zero
*/

template<class T>
void matrix<T>::setZero()
{
	memset(storage, 0, nCols * nRows * sizeof(T));
}

/** 
* @brief set all values to value
*/

template<class T>
void matrix<T>::set(T value)
{
	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] = value;
}

/** 
* @brief as = operator, avoiding memory allocation when copying
*/

template<class T>
void matrix<T>::set(const matrix<T>& source)
{
	if (source.rows() == nRows && source.cols() == nCols)
		memcpy(storage, source.data(), nRows * nCols * sizeof(T));
	else
		printf("set matrix: dimension do not agree. Nothing done\n");

}
/** 
* @brief fill matrix with random values of range min:max (default: 0:1)
*/

template<class T>
void matrix<T>::setRandom(T min, T max)
{
	if (min > max){
		T tmp = min;
		min = max;
		max = tmp;
	}
	//srand(time(NULL));
	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] = ((T)(rand()) / RAND_MAX) * (max - min) + min;
}

/** 
* @brief print matrix to console
*/

template<class T>
void matrix<T>::print() const
{

	std::cout << std::endl;

	std::cout << "size: " << nRows << " x " << nCols << std::endl;	
	for (uint i = 0; i < nRows; ++i)
	{
		for (uint j = 0; j < nCols; ++j)
			printf("%4.2f ", storage[i*nCols + j]);
			//std::cout << storage[i*nCols + j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/** 
* @brief print size of matrix to console
*/

template<class T>
void matrix<T>::size() const
{
	std::cout << "size: " << nRows << " x " << nCols << std::endl;	
}

/** 
* @brief resize matrix, existing values are copied and keep remained
*/

template<class T>
void matrix<T>::resize(uint newRows, uint newCols)
{

	if((newRows == nRows) && (newCols == nCols)) return;

	#if N_DEBUG
	if ((newRows < nRows) || (newCols < nCols))
		printf("warning: possible data lost due to resizing to smaller matrix\n");
	#endif 

	T* newData = new T[newRows * newCols];
	memset(newData, 0, newCols * newRows * sizeof(T));

	uint colsToCopy, rowsToCopy;
	
	newRows < nRows? rowsToCopy = newRows : rowsToCopy = nRows;
	newCols < nCols? colsToCopy = newCols : colsToCopy = nCols;

	for (uint i = 0; i < rowsToCopy; ++i)
		for (uint j = 0; j < colsToCopy; ++j)
			newData[i * newCols + j] = storage[i * nCols + j];

	nRows = newRows;
	nCols = newCols;

	delete[] storage;
	storage = newData;	
	postResize();
}


/** 
* @brief resize vector, keep existing values
*/

template<class T>
void matrix<T>::resize(uint newElements)
{
	if((nRows > 1) && (nCols > 1)) 
	{
		printf("matrix resizing needs two dimensions, no resizing done\n");
		return;
	}		

	if (newElements == 0)
	{
		this->clear();
		return;
	}

	uint l, elementsToCopy;
	nRows < nCols ? l = nCols : l = nRows;
	l < newElements ? elementsToCopy = l : elementsToCopy = newElements;

	#if N_DEBUG
	if (newElements < l)
		printf("warning: possible data lost due to resizing to smaller vector\n");
	#endif 

	T* newData = new T[newElements];
	memset(newData, 0, newElements * sizeof(T));
	
	memcpy(newData, storage, elementsToCopy * sizeof(T));

	nRows <= 1 ? nCols = newElements : nRows = newElements;
	if (nRows == 0) nRows = 1;
	if (nCols == 0) nCols = 1;
	

	delete[] storage;
	storage = newData;	
	postResize();
}

/** 
* @brief append row at the bottom of matrix
*/

template<class T>
void matrix<T>::addRow(const matrix<T>& newRow)
{

	if (nCols != newRow.length()) 
	{
		printf("dimensions do not agree, row not added\n");
		return;
	}

	T* newData = new T[(nRows+1) * nCols];
	memcpy(newData, storage, nRows * nCols * sizeof(T));
	memcpy(&newData[nRows * nCols], newRow.data(), nCols * sizeof(T));

	nRows++;
	
	delete[] storage;
	storage = newData;
	postResize();

}

/** 
* @brief append column at the right end of matrix
*/

template<class T>
void matrix<T>::addCol(const matrix<T>& newCol)
{

	if (nRows != newCol.length()) 
	{
		printf("dimensions do not agree, col not added\n");
		return;
	}

	T* newData = new T[(nCols+1) * nRows];
	
	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			newData[i*(nCols+1) + j] = storage[i*nCols + j];

	for (uint i = 0; i < nRows; ++i)
		newData[(i+1)*(nCols+1) -1] = newCol.element(i);


	nCols++;
	
	delete[] storage;
	storage = newData;
	postResize();

}

/** 
* @brief append multiple rows at the end of matrix
*/

template<class T>
void matrix<T>::addMatrixBottom(const matrix<T>& newBottom)
{

	if (this->length() == 0)
	{	
		(*this) = newBottom;
		return;
	}

	if (nCols != newBottom.cols()) 
	{
		printf("dimensions do not agree, matrix not added\n");
		return;
	}

	T* newData = new T[(nRows+newBottom.rows()) * nCols];
	
	memcpy(newData, storage, nRows * nCols * sizeof(T));
	memcpy(&newData[nRows * nCols], newBottom.data(),newBottom.length()*sizeof(T));

	nRows += newBottom.rows();

 	delete[] storage;
	storage = newData;
	postResize();
}	

/** 
* @brief append multiple columns at the right end of matrix
*/

template<class T>
void matrix<T>::addMatrixRight(const matrix<T>& newRight)
{

	if (nRows != newRight.rows()) 
	{
		printf("dimensions do not agree, matrix not added\n");
		return;
	}

	T* newData = new T[(nCols+newRight.cols()) * nRows];

	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			newData[i*(nCols+newRight.cols()) + j] = storage[i*nCols + j];

	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < newRight.cols(); ++j)
			newData[i*(nCols+newRight.cols()) + nCols + j] = newRight.element(i*newRight.cols() + j);

	nCols += newRight.cols();
	
	delete[] storage;
	storage = newData;
	postResize();

}

/** 
* @brief append element at the end of vector
*/

template<class T>
void matrix<T>::addElement(T newElement)
{

	if((nRows > 1) && (nCols > 1)) 
	{
		printf("need vector to add a single element, no adding done\n");
		return;
	}		

	T* newData = new T[nCols * nRows + 1];
		
	memcpy(newData, storage, nCols * nRows * sizeof(T));
		
	delete[] storage;
	storage = newData;

	storage[nCols * nRows] = newElement;
	
	nCols > nRows ? nCols++ : nRows ++;	
	
	if (nCols == 0) nCols ++;
	if (nRows == 0) nRows ++;
	postResize();

}	

/** 
* @brief reshape matrix while keeping the order of elements (matlab style)
*/

template<class T>
void matrix<T>::reshape(uint newRows, uint newCols)
{

	if ((newRows * newCols) != (nRows * nCols))
	{
		printf("reshape: dimensions do not fit. No reshaping done\n");
		return;
	}

	nRows = newRows;
	nCols = newCols;
	postResize();
}

/** 
* @brief return n-th row
*/

template<class T>
matrix<T> matrix<T>::getRow(uint n) const
{

	matrix<T> selectedRow(1,nCols,&storage[n*nCols]);
	return selectedRow;
}

/** 
* @brief return n-th column
*/

template<class T>
matrix<T> matrix<T>::getCol(uint n) const
{

	matrix<T> selectedCol(nRows,1);
		
	for (uint i = 0; i < nRows; ++i)
		selectedCol(i) = storage[i*nCols + n];

	return selectedCol;
}

/** 
* @brief return pointer to n-th row
*/

template<class T>
T* matrix<T>::row(uint n) const
{
	return &storage[n*nCols];
}




/** 
* @brief set n-th row to newRow (raw data array)
*/

template<class T>
void matrix<T>::setRow(uint n, T* newRow) 
{
	if (n > nRows)
	{
		std::cout << "row " << n << " out of range "<< nRows << " , nothing copied" << std::endl;
		return;
	}
	memcpy(&storage[n*nCols], newRow, nCols * sizeof(T));
}

/** 
* @brief set n-th row to newRow (matrix vector)
*/

template<class T>
void matrix<T>::setRow(uint n, const matrix<T>& newRow) 
{
	if (n > nRows)
	{
		std::cout << "row " << n << " out of range "<< nRows << " , nothing copied" << std::endl;
		return;
	}

	if (newRow.cols() != nCols)
	{
		std::cout << "no of columns does not match, row not copied" << std::endl;
		return;
	}

	memcpy(&storage[n*nCols], newRow.data(), nCols * sizeof(T));
}

/** 
* @brief set n-th column to newCol (raw data array)
*/

template<class T>
void matrix<T>::setCol(uint n, T* newCol)
{
	if (n > nCols)
	{
		std::cout << "column " << n << " out of range "<< nCols << " , nothing copied" << std::endl;
		return;
	}
	
	for (uint i = 0; i < nRows; ++i)
		storage[n + i*nCols] = newCol[i];
}

/** 
* @brief fill with linear spaced values from from to to.
*/

template<class T>
void matrix<T>::linspace(T from, T to)
{
	if (!std::isfinite(to)) to = (T) (nCols * nRows + from -1);
	
	T diff = (to-from)/(nCols * nRows - 1);

	for (uint i = 0; i < nRows * nCols; ++i)
		storage[i] = from + i*diff;
}

/** 
* @brief fill each row with linear spaced values from from to to .
*/

template<class T>
void matrix<T>::linspacePerRow(T from, T to)
{
	if (!std::isfinite(to)) to = (T) (nCols + from -1);
	
	T diff = (to-from)/(nCols-1);

	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			(*this)(i,j) = from + j * diff;

}

/** 
* @brief fill each column with linear spaced values from from to to .
*/

template<class T>
void matrix<T>::linspacePerCol(T from, T to)
{
	if (!std::isfinite(to)) to = (T) (nRows + from -1);
	
	T diff = (to-from)/(nRows-1);

	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			(*this)(i,j) = from + i * diff;

}

/** 
* @brief set all values below thresh to zero, others to one  
*/

template<class T>
void matrix<T>::binarize(T thresh)
{
	
	for (uint i = 0; i < this->length(); ++i)
		if (storage[i] < thresh) 
			storage[i] = 0.;
		else storage[i] = 1;
}

/** 
* @brief set all values below thresh to zero, leave others  
*/

template<class T>
void matrix<T>::threshold(T thresh)
{
	
	for (uint i = 0; i < this->length(); ++i)
		if (storage[i] < thresh) 
			storage[i] = 0.;
}



/** 
* @brief set n-th column to newCol (MatrixBase vector)
*/

template<class T>
void matrix<T>::setCol(uint n, const matrix<T>& newCol)
{
	if (n > nCols)
	{
		std::cout << "column " << n << " out of range "<< nCols << " , nothing copied" << std::endl;
		return;
	}
	
	if (newCol.rows() != nRows)
	{
		std::cout << "no of rows does not match, column not copied" << std::endl;
		return;
	}
	
	for (uint i = 0; i < nRows; ++i)
		storage[n + i*nCols] = newCol(i);
}

/** 
* @brief assignment operator
*/


template<class T>
matrix<T> matrix<T>::operator=(const matrix<T>& rhs)
{
	if (&rhs == this) return *this;

	if (rhs.rows() == nRows && rhs.cols() == nCols)
	{
		memcpy(storage, rhs.data(), nCols * nRows * sizeof(T));
		return *this;
	}
	else
	{	
		delete[] storage;
		storage = new T[rhs.rows() * rhs.cols() * sizeof(T)];
	
		memcpy(storage, rhs.data(), rhs.cols() * rhs.rows() * sizeof(T));
	
		nCols = rhs.cols();
		nRows = rhs.rows();
		return *this;
		postResize();
	}
}

/** 
* @brief elementwise addition, return result. Automatic broadcasting of rhs if neccessary.
*/

template<class T>
matrix<T> matrix<T>::operator+(const matrix<T>& rhs) const
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting +" << std::endl;
			for (uint i = 0; i < nRows; i++)
				result.setRow(i,this->getRow(i) + rhs);
			return result;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting +" << std::endl;
			for (uint i = 0; i < nCols; i++)
				result.setCol(i,this->getCol(i) + rhs);
			return result;
		}

	}


	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no addition performed\n");
		return *this;
	}

	matrix<T> result(nRows, nCols);

	for (uint i = 0; i < nCols * nRows; ++i)
		result.data()[i] = storage[i] + rhs.data()[i];

	return result;
}

/** 
* @brief elementwise subtraction, return result. Automatic broadcasting of rhs if neccessary.
*/

template<class T>
matrix<T> matrix<T>::operator-(const matrix<T>& rhs) const
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting -" << std::endl;
			for (uint i = 0; i < nRows; i++)
				result.setRow(i,this->getRow(i) - rhs);
			return result;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting -" << std::endl;
			for (uint i = 0; i < nCols; i++)
				result.setCol(i,this->getCol(i) - rhs);
			return result;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no subtraction performed\n");
		return *this;
	}

	matrix<T> result(nRows, nCols);

	for (uint i = 0; i < nCols * nRows; ++i)
		result.data()[i] = storage[i] - rhs.data()[i];

	return result;
}

/** 
* @brief elementwise multiplication, return result
*/

template<class T>
matrix<T> matrix<T>::operator*(const matrix<T>& rhs) const
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting *" << std::endl;
			for (uint i = 0; i < nRows; i++)
				result.setRow(i,this->getRow(i) * rhs);
			return result;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting *" << std::endl;
			for (uint i = 0; i < nCols; i++)
				result.setCol(i,this->getCol(i) * rhs);
			return result;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no multiplication performed\n");
		return *this;
	}

	matrix<T> result(nRows, nCols);

	for (uint i = 0; i < nCols * nRows; ++i)
		result.data()[i] = storage[i] * rhs.data()[i];

	return result;
}

/** 
* @brief elementwise division, return result
*/

template<class T>
matrix<T> matrix<T>::operator/(const matrix<T>& rhs) const
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting /" << std::endl;
			for (uint i = 0; i < nRows; i++)
				result.setRow(i,this->getRow(i) / rhs);
			return result;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			matrix<T> result(nRows, nCols);
			std::cout << "warning: automatic broadcasting /" << std::endl;
			for (uint i = 0; i < nCols; i++)
				result.setCol(i,this->getCol(i) / rhs);
			return result;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no division performed\n");
		return *this;
	}

	matrix<T> result(nRows, nCols);

	for (uint i = 0; i < nCols * nRows; ++i)
		result.data()[i] = storage[i] / rhs.data()[i];

	return result;
}

/** 
* @brief add rhs to matrix (in place addition)
*/

template<class T>
void matrix<T>::operator+=(const matrix<T>& rhs)
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			std::cout << "warning: automatic broadcasting +=" << std::endl;
			for (uint i = 0; i < nRows; i++)
				this->setRow(i,this->getRow(i) + rhs);
			return;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			std::cout << "warning: automatic broadcasting +=" << std::endl;
			for (uint i = 0; i < nCols; i++)
				this->setCol(i,this->getCol(i) + rhs);
			return;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no addition performed\n");
		return ;
	}

	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] += rhs.data()[i];
}

/** 
* @brief multiply matrix with rhs elementwise (in place multiplication)
*/

template<class T>
void matrix<T>::operator*=(const matrix<T>& rhs)
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			std::cout << "warning: automatic broadcasting *=" << std::endl;
			for (uint i = 0; i < nRows; i++)
				this->setRow(i,this->getRow(i) * rhs);
			return;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			std::cout << "warning: automatic broadcasting *=" << std::endl;
			for (uint i = 0; i < nCols; i++)
				this->setCol(i,this->getCol(i) * rhs);
			return;
		}

	}

	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no multiplication performed\n");
		return;
	}

	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] *= rhs.data()[i];
}

/** 
* @brief multiply each element in matrix with value (in place multiplication)
*/

template<class T>
void matrix<T>::operator*=(T value)
{
	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] *= value;
}

/** 
* @brief subtract rhs from matrix (in place subtraction)
*/

template<class T>
void matrix<T>::operator-=(const matrix<T>& rhs)
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			std::cout << "warning: automatic broadcasting -=" << std::endl;
			for (uint i = 0; i < nRows; i++)
				this->setRow(i,this->getRow(i) - rhs);
			return;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			std::cout << "warning: automatic broadcasting -=" << std::endl;
			for (uint i = 0; i < nCols; i++)
				this->setCol(i,this->getCol(i) - rhs);
			return;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no subtraction performed\n");
		return;
	}

	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] -= rhs.data()[i];
}

/** 
* @brief divide matrix by rhs elementwise (in place division)
*/

template<class T>
void matrix<T>::operator/=(const matrix<T>& rhs)
{
	// check for automatic broadcasting first
	if (rhs.isVector() && !this->isVector())
	{
		if (rhs.rows() == 1 && rhs.cols() == nCols)
		{
			std::cout << "warning: automatic broadcasting /=" << std::endl;
			for (uint i = 0; i < nRows; i++)
				this->setRow(i,this->getRow(i) / rhs);
			return;
		}
		if (rhs.rows() == nRows && rhs.cols() == 1)
		{
			std::cout << "warning: automatic broadcasting /=" << std::endl;
			for (uint i = 0; i < nCols; i++)
				this->setCol(i,this->getCol(i) / rhs);
			return;
		}

	}
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
	{
		printf("dimensions do not agree, no division performed\n");
		return;
	}

	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] /= rhs.data()[i];
}

/** 
* @brief divide each element of matrix by value (in place division)
*/

template<class T>
void matrix<T>::operator/=(T value)
{
	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] /= value;
}




/** 
* @brief comparision operator, true if no element pairs differ more than EPS (1e-10)
*/

template<class T>
bool matrix<T>::operator==(const matrix<T>& rhs) const
{
	if ((nCols != rhs.cols()) || nRows != rhs.rows())
		return 0;

	for (uint i = 0; i < nCols * nRows; ++i)
		if (std::abs(storage[i] - rhs.data()[i]) > EPS) return 0;

	return 1;
}

/** 
* @brief access operator [] to return pointer to element i
*/

template<class T>
T* matrix<T>::operator[](const uint i) const
{
	if (i >= nRows * nCols) 
	{
		printf("index out of range\n");
		return 0;
	}

	return &storage[i];
}

/** 
* @brief access operator () return value of i-th element
*/

template<class T>
T matrix<T>::operator()(const uint i) const
{
	if (i >= nRows * nCols) 
	{
		printf("index out of range\n");
		return 0;
	}

	return storage[i];
}

/** 
* @brief return reference of i-th element (for writing)
*/

template<class T>
T& matrix<T>::operator()(const uint i)
{
	if (i >= nRows * nCols) 
	{
		printf("index out of range\n");
		return storage[nCols * nRows -1];
	}

	return storage[i];
}

/** 
* @brief return element(i,j)
*/

template<class T>
T matrix<T>::operator()(const uint i, const uint j) const
{
	if ((i >= nRows) || (j >= nCols))
	{
		printf("index out of range\n");
		return 0;
	}

	return storage[i*nCols + j];
}

/** 
* @brief return element(i,j) for writing
*/

template<class T>
T& matrix<T>::operator()(const uint i, const uint j)
{
	if ((i >= nRows) || (j >= nCols))
	{
		printf("index out of range\n caution: last element manipulated\n");
		return storage[nCols * nRows-1];
	}

	return storage[i*nCols + j];
}

/** 
* @brief set each element to abs(element)
*/

template<class T>
void matrix<T>::abs()
{
	for (uint i = 0; i < nCols * nRows; ++i)
		storage[i] = std::abs(storage[i]);
}

/** 
* @brief return transposed matrix
*/

template<class T>
matrix<T> matrix<T>::transpose() const
{

	matrix<T> transposed(nCols, nRows);
	
	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			transposed.data()[j*nRows + i] = storage[i*nCols + j];

	return transposed;

}

/** 
* @brief perform in place transposition
*/

template<class T>
void matrix<T>::transposeInPlace()
{
	if (nCols == 1 || nRows == 1)
	{
		std::swap(nCols, nRows);
		return;
	}
	matrix<T> trans = this->transpose();
	
	memcpy(storage, trans.data(), nCols * nRows * sizeof(T));
	
	nCols = nRows;
	nRows = trans.rows();
	postResize();
}


/** 
* @brief return true if all elements are nonzero
*/

template<class T>
bool matrix<T>::all() const
{
	for (uint i = 0; i < nCols * nRows; ++i)
		if (storage[i] == 0.) return 0;

	return 1;
}

/** 
* @brief return true if any element is nonzero
*/

template<class T>
bool matrix<T>::any() const
{
	for (uint i = 0; i < nCols * nRows; ++i)
		if (storage[i] != 0.) return 1;

	return 0;
}

/** 
* @brief return true if matrix contains no NANs and no INFs
*/

template<class T>
bool matrix<T>::allFinite() const
{
	for (uint i = 0; i < nCols * nRows; ++i)
		if (!std::isfinite(storage[i])) return 0;

	return 1;
}

/** 
* @brief return mean of all elements
*/

template<class T>
T matrix<T>::mean() const
{
	T sum = 0;
	for (uint i = 0; i < nCols * nRows; ++i)
		sum += storage[i];

	return (sum/(nCols * nRows));
}

/** 
* @brief return sum of all elements
*/

template<class T>
T matrix<T>::sum() const
{
	T sum = 0;
	for (uint i = 0; i < nCols * nRows; ++i)
		sum += storage[i];

	return sum;
}

/** 
* @brief return max value of all elements
*/

template<class T>
T matrix<T>::max() const
{
	T max = storage[0];
	for (uint i = 1; i < nCols * nRows; ++i)
		if (storage[i] > max) max = storage[i];

	return max;
}

/** 
* @brief return index of max value of all elements
*/

template<class T>
T matrix<T>::maxIndex() const
{
	T max = storage[0];
	uint maxInd = 0;
	for (uint i = 1; i < nCols * nRows; ++i)
		if (storage[i] > max) 
		{
			max = storage[i];
			maxInd = i;
		}	

	return maxInd;
}


/** 
* @brief return min value of all elements
*/

template<class T>
T matrix<T>::min() const
{
	T min = storage[0];
	for (uint i = 1; i < nCols * nRows; ++i)
		if (storage[i] < min) min = storage[i];

	return min;
}

/** 
* @brief return index of min value of all elements
*/

template<class T>
T matrix<T>::minIndex() const
{
	T min = storage[0];
	uint minInd = 0;
	for (uint i = 1; i < nCols * nRows; ++i)
		if (storage[i] < min) 
		{
			min = storage[i];
			minInd = i;
		}	

	return minInd;
}


/** 
* @brief load matrix from file, discard existing data
*/

template<class T>
void matrix<T>::load(std::string pathFilename) 
{
	readFromFile(pathFilename);
	postResize();
}

/** 
* @brief save matrix to file
*/

template<class T>
void matrix<T>::save(std::string pathFilename) const
{
	if ((nRows == 0) || (nCols == 0)) 
	{
		std::cout << "no data to write to file " << pathFilename << std::endl;
		return;
	}

	std::ofstream out(pathFilename.c_str(), std::ios::out | std::ios::trunc);

	uint n = 0;

	if (out.is_open())
	{
		for (uint i = 0; i < nRows; ++i)
		{
			for (uint j = 0; j < nCols; ++j)
			{
				if (j < nCols-1)
					out << storage[n] << ",";
				else out << storage[n];
				n++;
			}
			out << "\n";
		}	
		out.close();
	}
	else std::cout << "could not open file " << pathFilename << " to write" << std::endl;
}

/** 
* @brief get vector with sum of each row
*/

template<class T>
matrix<T> matrix<T>::sumPerRow() const
{
	matrix<T> result(nRows, 1);
	T sum = 0;
	for (uint i = 0; i < nRows; ++i)
	{
		for (uint j = 0; j < nCols; ++j)
			sum += storage[i*nCols + j];
		result.data()[i] = sum;
		sum = 0;
	}
	
	return result;
}

/** 
* @brief get vector with mean of each row
*/

template<class T>
matrix<T> matrix<T>::meanPerRow() const
{
	matrix<T> result(nRows, 1);
	T sum = 0;
	for (uint i = 0; i < nRows; ++i)
	{
		for (uint j = 0; j < nCols; ++j)
			sum += storage[i*nCols + j];
		result.data()[i] = sum / nCols;
		sum = 0;
	}
	
	return result;
}

/** 
* @brief get vector with sum of each column
*/

template<class T>
matrix<T> matrix<T>::sumPerCol() const
{
	matrix<T> result(1, nCols);
	T sum = 0;
	for (uint i = 0; i < nCols; ++i)
	{
		for (uint j = 0; j < nRows; ++j)
			sum += storage[j*nCols + i];
		result.data()[i] = sum;
		sum = 0;
	}
	
	return result;
}

/** 
* @brief get vector with mean of each column
*/

template<class T>
matrix<T> matrix<T>::meanPerCol() const
{
	matrix<T> result(1, nCols);
	T sum = 0;
	for (uint i = 0; i < nCols; ++i)
	{
		for (uint j = 0; j < nRows; ++j)
			sum += storage[j*nCols + i];
		result.data()[i] = sum / nRows;
		sum = 0;
	}
	
	return result;
}

/** 
* @brief get vector with max value of each row, second col of return vector contains indices
*/

template<class T>
matrix<T> matrix<T>::maxPerRow() const
{
	matrix<T> result(nRows, 2);
	T max = 0;
	uint maxInd = 0;
	for (uint i = 0; i < nRows; ++i)
	{
		max = storage[i*nCols];
		maxInd = 0;
		for (uint j = 0; j < nCols; ++j)
			if (storage[i*nCols + j] > max) 
			{
				max = storage[i*nCols + j];
				maxInd = j;
			}
		result(i,0) = max;
		result(i,1) = maxInd;
	}
	
	return result;
}

/** 
* @brief get vector with min value of each row, second col of return vector contains indices
*/

template<class T>
matrix<T> matrix<T>::minPerRow() const
{
	matrix<T> result(nRows, 2);
	T min = 0;
	uint minInd = 0;
	for (uint i = 0; i < nRows; ++i)
	{
		min = storage[i*nCols];
		minInd = 0;
		for (uint j = 0; j < nCols; ++j)
			if (storage[i*nCols + j] < min) 
			{
				min = storage[i*nCols + j];
				minInd = j;
			}
		result(i,0) = min;
		result(i,1) = minInd;
	}
	
	return result;
}

/** 
* @brief get vector with min value of each column, second row contains indices
*/

template<class T>
matrix<T> matrix<T>::minPerCol() const
{
	matrix<T> result(2, nCols);
	T min = 0;
	uint minInd = 0;
	for (uint i = 0; i < nCols; ++i)
	{
		min = storage[i];
		minInd = 0;
		for (uint j = 0; j < nRows; ++j)
			if (storage[j*nCols + i] < min) 
			{
				min = storage[j*nCols + i];
				minInd = j;
			}
		result(0,i) = min;
		result(1,i) = minInd;
	}
	
	return result;
}

/** 
* @brief get vector with max value of each column, second row contains indices
*/

template<class T>
matrix<T> matrix<T>::maxPerCol() const
{
	matrix<T> result(2, nCols);
	T max = 0;
	uint maxInd = 0;
	for (uint i = 0; i < nCols; ++i)
	{
		max = storage[i];
		maxInd = 0;
		for (uint j = 0; j < nRows; ++j)
			if (storage[j*nCols + i] > max) 
			{
				max = storage[j*nCols + i];
				maxInd = j;
			}
		result(0,i) = max;
		result(1,i) = maxInd;
	}
	
	return result;
}


/** 
* @brief normalize matrix to max value = 1
*/

template<class T>
void matrix<T>::normalize()
{
	T maxValue = max();
	(*this)/=maxValue;
}

/** 
* @brief normalize each row to max value = 1
*/

template<class T>
void matrix<T>::normalizePerRow()
{
	matrix<T> maxValues = maxPerRow();
	
	for (uint i = 0; i < nRows; ++i)
		for (uint j = 0; j < nCols; ++j)
			storage[i*nCols + j] /= maxValues(i);
}

/** 
* @brief normalize each column to max value = 1
*/

template<class T>
void matrix<T>::normalizePerCol()
{
	matrix<T> maxValues = maxPerCol();

	for (uint i = 0; i < nCols; ++i)
		for (uint j = 0; j < nRows; ++j)
			storage[j*nCols + i] /= maxValues(i);
}


/**
* @brief run convolution, allocates memory, not yet real time capable
*
*/

template<class T>
void matrix<T>::convBase(T* a, uint lengthA, T* b, uint lengthB, T* result)
{

	uint l = lengthA + lengthB -1;
//	T* result = new T[l];
	memset(result, 0, l*sizeof(T));

	for (uint i = 0; i < l; ++i)
		for (uint j = 0; j < std::max(lengthA, lengthB); ++j)
			result[i] += (j < lengthA ? a[j] : 0) * (i-j < lengthB? b[i-j] : 0); 


}

/**
* @brief smooth in place with gaussian kernel of length len. Allocates new memory, not yet real time capable
*
*/

template<class T>
void matrix<T>::smoothInPlace(uint len)
{
	
	matrix<T> kernel(len);
	
	// calculate hamming coefficients
	for (uint i = 0; i < len; ++i)
		kernel(i) = (0.54-0.46*cos((i*2*M_PI)/len));
	
	matrix<T> result(len + nCols -1);
	
	for (uint i = 0; i < nRows; ++i)
	{
		convBase(storage, nCols, kernel.data(), len, result.data());
		this->setRow(i, result[(uint)round(len/2)]);	
		//this->setRow(i, result.data()[(uint)round(len/2)]);	
	}	

	(*this)/=kernel.sum();
}

/**
* @brief return smoothed copy (with gaussian kernel of length len). Allocates new memory, not yet real time capable
*
*/

template<class T>
matrix<T> matrix<T>::smooth(uint len)
{
	
	matrix<T> kernel(len);
	
	// calculate hamming coefficients
	for (uint i = 0; i < len; ++i)
		kernel(i) = (0.54-0.46*cos((i*2*M_PI)/len));
	
	matrix<T> resultRow(len + nCols -1);
	matrix<T> result(nRows, nCols);
	
	for (uint i = 0; i < nRows; ++i)
	{
		convBase(storage, nCols, kernel.data(), len, resultRow.data());
		result.setRow(i, resultRow[(uint)round(len/2)]);	
		//this->setRow(i, result.data()[(uint)round(len/2)]);	
	}	

	result/=kernel.sum();
	return result;
}

/**
* @brief matrix multiplication result = a*b using blas, memory for result matrix is allocated*
*/

template<class T>
matrix<T> matrix<T>::multiply(const matrix<T>& a, const matrix<T>& b, bool transA, bool transB)
{
	
	// Problem: blas needs columnwise storage, these matrices are stored rowwise.
	// solution: (A*B)' = B' * A' 
	
	uint m,n,k,l;
	m = a.rows();
	n = a.cols();
	k = b.rows();
	l = b.cols();
	if (transA) std::swap(m,n);
	if (transB) std::swap(k,l);

	// storae correction:
	//if (!transA) std::swap(m,n);
	//if (!transB) std::swap(k,l);


	if (n != k)
	{
		printf("matrix multiplication: dimensions do not fit\n");
		matrix<T> a;
		return a;
	}

	matrix<T> result(m, l);
	multiplyMatrix(transB, transA, l, m, k, b.data(), a.data(), result.data()); 
	
	return result;

}

/**
* @brief matrix multiplication c = a*b using blas, taking result matrix to avoid memory allocation
*/

template<class T>
void matrix<T>::multiply(const matrix<T>& a, const matrix<T>& b, matrix<T>& c, bool transA, bool transB)
{
	
	// Problem: blas needs columnwise storage, these matrices are stored rowwise.
	// solution: (A*B)' = B' * A' 
	
	uint m,n,k,l;
	m = a.rows();
	n = a.cols();
	k = b.rows();
	l = b.cols();
	if (transA) std::swap(m,n);
	if (transB) std::swap(k,l);

	if (n != k)
		printf("matrix multiplication: factors dimensions do not fit\n");

	if (c.rows() != m || c.cols()!=l)
		printf("matrix multiplication: result dimensions do not fit\n");
	
	multiplyMatrix(transB, transA, l, m, k, b.data(), a.data(), c.data()); 
}

/**
* @brief return linear indices of nonzero elements. 
*/

template<class T>
matrix<T> matrix<T>::find() const
{
	std::cout << "deprecated function: find()" << std::endl;
	matrix<T> result;
	for (uint i = 0; i < nCols * nRows; ++i)
		if (storage[i] != 0) result.addElement((T) i);

	return result;
}

template class matrix<float>;
template class matrix<double>;

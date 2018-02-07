#ifndef _MATRICES_H
#define _MATRICES_H 1

#include <string>

typedef unsigned int uint;

/**
* @brief matrix template class for double and float matrices and vectors
+ @author nicole
*/

template<class T = double>
class matrix{

	protected:
	
	/** @brief number of rows*/
	uint nRows;
	/** @brief number of columns*/
	uint nCols;
	/** @brief pointer to contiguous memory storing matrixdata*/
	T* storage;
	/** @brief amount to increase storage if one elemen is added. Increases if single elements are often added*/
	uint toResize = 1;
	/** @brief allocated memory (elements)*/
	uint allocated = 0;
	
	/** @brief epsilon. (x < EPS -> x == 0)*/
	const T EPS = 1e-5;

	void prepareForFileread(std::string pathFilename);
	void readFromFile(std::string pathFilename);
	void convBase(T* a, uint lengthA, T* b, uint lengthB, T* result);

	virtual void postResize();

    public:

	// constructors
	matrix();
	matrix(uint n, uint m);
	matrix(uint n, uint m, T* dataInit);
	matrix(uint n);
	matrix(uint n, T* dataInit);
	matrix(std::string pathFilename);

	// copyConstructor
	matrix(const matrix<T>& matrixInit);

	// destructor
	~matrix();

	// delete
	void clear();

	// get properties	
	uint rows() const;
	uint cols() const;
	uint length() const;
	T* data() const;
	T element(uint n) const;
	bool isVector() const;

	// initialize
	void setZero();
	void set(T value);
	void set(const matrix<T>& source);
	void setRandom(T min = 0., T max = 1.);

	// resizing and concenating
	void resize(uint newRows, uint newCols);
	void resize(uint newElements);
	void addRow(const matrix<T>& newRow);
	void addCol(const matrix<T>& newCol);
	void addElement(T newElement);
	void addMatrixBottom(const matrix<T>& newBottom);
	void addMatrixRight(const matrix<T>& newRight);
	void reshape(uint newRows, uint newCols);

	void print() const;
	void size() const;

	matrix<T> getRow(uint n) const;	
	matrix<T> getCol(uint n) const;	
	T* row(uint n) const;
	void setRow(uint n, T* data);	
	void setRow(uint n, const matrix<T>& data);	
	void setCol(uint n, T* data);	
	void setCol(uint n, const matrix<T>& data);	

	void linspace(T from=1., T to=(0.0/0.0));
	void linspacePerRow(T from=1., T to=(0.0/0.0));
	void linspacePerCol(T from=1., T to=(0.0/0.0));
	void binarize(T thresh);
	void threshold(T thresh);


	// assignment operator
	matrix<T> operator=(const matrix<T>& rhs);

	// elementwise operators
	matrix<T> operator+(const matrix<T>& rhs) const;
	matrix<T> operator-(const matrix<T>& rhs) const;
	matrix<T> operator*(const matrix<T>& rhs) const;
	matrix<T> operator/(const matrix<T>& rhs) const;
	void operator+=(const matrix<T>& rhs);
	void operator*=(const matrix<T>& rhs);
	void operator*=(T value);
	void operator-=(const matrix<T>& rhs);
	void operator/=(const matrix<T>& rhs);
	void operator/=(T value);
	bool operator==(const matrix<T>& rhs) const;	
	T* operator[](const uint i) const;
	T operator()(const uint i) const;
	T& operator()(const uint i);
	T operator()(const uint i, const uint j) const;
	T& operator()(const uint i, const uint j);
	void abs();
	static matrix<T> multiply(const matrix<T>& a, const matrix<T>& b, bool transA = 0, bool transB = 0);
	static void multiply(const matrix<T>& a, const matrix<T>& b, matrix<T>& c, bool transA = 0, bool transB = 0);


	matrix<T> transpose() const;
	void transposeInPlace();

	bool all() const;
	bool any() const;
	bool allFinite() const;
	T mean() const;
	T sum() const;
	T max() const;
	T maxIndex() const;
	T min() const;
	T minIndex() const;
	matrix<T> sumPerRow() const;
	matrix<T> meanPerRow() const;
	matrix<T> sumPerCol() const;
	matrix<T> meanPerCol() const;
	matrix<T> maxPerRow() const;
	matrix<T> minPerRow() const;
	matrix<T> maxPerCol() const;
	matrix<T> minPerCol() const;


	// file i/o
	void load(std::string pathFilename);
	void save(std::string pathFilename) const;


	void normalize();
	void normalizePerRow();
	void normalizePerCol();
	
	void smoothInPlace(uint len);
	matrix<T> smooth(uint len);

	matrix<T> find() const;

};




#endif


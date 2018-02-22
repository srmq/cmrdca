/*
 * GenDissimMatrix.h
 *
 *  Created on: Feb 21, 2018
 *      Author: srmq
 */

#ifndef SRC_GENDISSIMMATRIX_H_
#define SRC_GENDISSIMMATRIX_H_

namespace util {

class IDissimMatrix {
public:
	virtual double getDissim(const int i, const int j) const = 0;
	virtual void putDissim(const int i, const int j, const double dissim) = 0;
	virtual int length() const = 0;
	virtual ~IDissimMatrix() = 0;
};

inline IDissimMatrix::~IDissimMatrix() { }

template<class T = double>
class GenDissimMatrix : public IDissimMatrix {
private:
	T **matrix;
	const int nElems;
	T maxDissim;

public:
	GenDissimMatrix(const int nElements, T maxDissim) : nElems(nElements), maxDissim(maxDissim) {
		this->matrix = new T*[nElems];
		for (int i = 0; i < nElements; i++) {
			this->matrix[i] = new T[i+1];
		}
	}
	GenDissimMatrix(const int nElements) : GenDissimMatrix(nElements, -1.0){}
	virtual ~GenDissimMatrix() {
		for (int i = 0; i < this->nElems; i++) {
			delete[] this->matrix[i];
		}
		delete[] this->matrix;
	}

	double getDissim(const int i, const int j) const {
		const T storedDissim = (j > i) ? this->matrix[j][i] : this->matrix[i][j];
		if (storedDissim > 0)
			return (double)storedDissim;
		else
			return (double)maxDissim;
	}

	void putDissim(const int i, const int j, const double dissim) {
		(j > i) ? this->matrix[j][i] = (T)dissim : this->matrix[i][j] = (T) dissim;
		if (dissim > (double)maxDissim)
			maxDissim = (T)dissim;
	}
	int length() const {return this->nElems;}

};

}



#endif /* SRC_GENDISSIMMATRIX_H_ */

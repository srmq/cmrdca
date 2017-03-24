/*
 * DissimMatrix.h
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#ifndef DISSIMMATRIX_H_
#define DISSIMMATRIX_H_

namespace util {

class DissimMatrix {
private:
	double **matrix;
	const int nElems;

public:
	DissimMatrix(const int nElements);
	virtual ~DissimMatrix();

	double getDissim(const int i, const int j) const {
		return (j > i) ? this->matrix[j][i] : this->matrix[i][j];
	}

	void putDissim(const int i, const int j, const double dissim) {
		(j > i) ? this->matrix[j][i] = dissim : this->matrix[i][j] = dissim;
	}
	int length() const {return this->nElems;}
};

} /* namespace util */
#endif /* DISSIMMATRIX_H_ */

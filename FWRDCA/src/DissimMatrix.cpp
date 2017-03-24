/*
 * DissimMatrix.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#include "DissimMatrix.h"

namespace util {

DissimMatrix::DissimMatrix(const int nElements) : nElems(nElements) {
	this->matrix = new double*[nElems];
	for (int i = 0; i < nElements; i++) {
		this->matrix[i] = new double[i+1];
	}

}

DissimMatrix::~DissimMatrix() {
	for (int i = 0; i < this->nElems; i++) {
		delete[] this->matrix[i];
	}
	delete[] this->matrix;
}

} /* namespace util */

/*
 * FuzzyCluster.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#include "FuzzyCluster.h"
#include <map>
#include <memory>
#include <cstring>
#include <cassert>
#include <set>

namespace util {

FuzzyCluster::FuzzyCluster(const unsigned int nelems, const int p) : nelems(nelems), p(p), center(-1), lambdaWeights(new double[p], std::default_delete<double[]>()), elements(new std::set<int>){
	std::fill_n(lambdaWeights.get(), p, 1.0/(double)p);
}

FuzzyCluster::FuzzyCluster(const FuzzyCluster& copyFrom) : nelems(copyFrom.nelems), p(copyFrom.p), center(copyFrom.center), lambdaWeights(new double[copyFrom.p], std::default_delete<double[]>()), membershipDegrees(copyFrom.membershipDegrees), elements(new std::set<int>(*(copyFrom.elements.get()))) {
	memcpy(this->lambdaWeights.get(), copyFrom.lambdaWeights.get(), (this->p)*sizeof(double));
}

FuzzyCluster::~FuzzyCluster() {

}

std::shared_ptr<double> FuzzyCluster::getWeights(int *p) {
	if (p != NULL) (*p) = this->p;
	return this->lambdaWeights;
}

void FuzzyCluster::setWeights(std::shared_ptr<double> weightPtr) {
	this->lambdaWeights = weightPtr;
}

double FuzzyCluster::getMembershipDegree(const int elem) const {
	double result = -1.0;
	const std::map<int,double>::const_iterator it = this->membershipDegrees.find(elem);
	if (it != this->membershipDegrees.end()) {
		result = it->second;
	}
	return result;
}

std::shared_ptr<std::set<int> > FuzzyCluster::getElements() const {
	assert(this->elements->size() == nelems);
	return this->elements;
}


} /* namespace util */


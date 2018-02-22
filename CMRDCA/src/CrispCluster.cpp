/*
 * CrispCluster.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#include "CrispCluster.h"
#include <map>
#include <memory>
#include <cstring>
#include <cassert>
#include <set>
#include <utility>

namespace util {

CrispCluster::CrispCluster(const int p) : p(p), medoids(new std::set<int>()), lambdaWeights(new double[p], std::default_delete<double[]>()), elements(new std::set<int>){
	std::fill_n(lambdaWeights.get(), p, 1.0 /*/(double)p */); //product of weights = 1
}

CrispCluster::CrispCluster(const CrispCluster& copyFrom) : p(copyFrom.p), medoids(new std::set<int>(*(copyFrom.medoids.get()))), lambdaWeights(new double[copyFrom.p], std::default_delete<double[]>()), elements(new std::set<int>(*(copyFrom.elements.get()))) {
	memcpy(this->lambdaWeights.get(), copyFrom.lambdaWeights.get(), (this->p)*sizeof(double));
}

CrispCluster::~CrispCluster() {

}

std::shared_ptr<double> CrispCluster::getWeights(int *p) {
	if (p != NULL) (*p) = this->p;
	return this->lambdaWeights;
}

void CrispCluster::setWeights(std::shared_ptr<double> weightPtr) {
	this->lambdaWeights = weightPtr;
}

std::shared_ptr<std::set<int> > CrispCluster::getElements() const {
	return this->elements;
}

bool CrispCluster::insert(const int element) {
	std::pair<std::set<int>::iterator, bool> result = this->elements.get()->insert(element);
	return result.second;
}

bool CrispCluster::remove(const int element) {
	bool found;
	std::set<int>::iterator pos = this->elements.get()->find(element);
	found = (pos != this->elements.get()->end());
	if (found) {
		this->elements.get()->erase(pos);
	}
	return found;
}


} /* namespace util */


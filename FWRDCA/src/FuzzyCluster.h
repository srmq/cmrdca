/*
 * FuzzyCluster.h
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#ifndef FUZZYCLUSTER_H_
#define FUZZYCLUSTER_H_

#include <map>
#include <memory>
#include <cassert>
#include <utility>
#include <set>

namespace util {

class FuzzyCluster {
public:
	FuzzyCluster(const unsigned int nelems, const int p);
	FuzzyCluster(const util::FuzzyCluster& copyFrom);
	virtual ~FuzzyCluster();
	int getCenter() const {
		return center;
	}
	void setCenter(int center) {
		this->center = center;
	}
	std::shared_ptr<double> getWeights(int *p);
	void setWeights(std::shared_ptr<double> weightPtr);

	void updateMemberhipDegree(const int elem, const double value) {
		this->membershipDegrees[elem] = value;
		this->elements->insert(elem);
	}
	double getMembershipDegree(int elem) const;
	double weightOf(int criterion) const {
		assert(criterion >= 0 && criterion < this->p);
		return this->lambdaWeights.get()[criterion];
	}
	std::shared_ptr<std::set<int>> getElements() const;

private:
	unsigned int nelems;
	const int p;
	int center;
	std::shared_ptr<double> lambdaWeights;
	std::map<int, double> membershipDegrees;
	std::shared_ptr<std::set<int>> elements;


};

} /* namespace util */
#endif /* FUZZYCLUSTER_H_ */

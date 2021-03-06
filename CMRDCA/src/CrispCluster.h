/*
 * CrispCluster.h
 *
 *  Created on: Jun 19, 2013
 *      Author: srmq
 */

#ifndef CRISPCLUSTER_H_
#define CRISPCLUSTER_H_

#include <map>
#include <memory>
#include <cassert>
#include <utility>
#include <set>

namespace util {

class CrispCluster {
public:
	CrispCluster(const int p);
	CrispCluster(const util::CrispCluster& copyFrom);
	virtual ~CrispCluster();
	std::shared_ptr<std::set<int> > getMedoids()  const {
		return medoids;
	}
	void setMedoids(std::shared_ptr<std::set<int> > newMedoids) {
		this->medoids = newMedoids;
	}
	std::shared_ptr<double> getWeights(int *p);
	void setWeights(std::shared_ptr<double> weightPtr);

	bool insert(int element);
	bool remove(int element);
	void setElements(std::shared_ptr<std::set<int> > elements) {
		this->elements = elements;
	}

	double weightOf(int criterion) const {
		assert(criterion >= 0 && criterion < this->p);
		return this->lambdaWeights.get()[criterion];
	}
	std::shared_ptr<std::set<int>> getElements() const;

private:
	const int p;
	std::shared_ptr<std::set<int> > medoids;
	std::shared_ptr<double> lambdaWeights;
	std::shared_ptr<std::set<int>> elements;


};

} /* namespace util */
#endif /* CRISPCLUSTER_H_ */

/*
 * CMRDCA.h
 *
 *  Created on: May 27, 2013
 *      Author: srmq
 */

#ifndef CMRDCA_H_
#define CMRDCA_H_

#include "DissimMatrix.h"
#include "CrispCluster.h"
#include <vector>
#include <memory>
#include <random>
#include <limits>
#include <ctime>
#include <iostream>
#include <sstream>


namespace clustering {

class CMRDCA {
public:
	CMRDCA(const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices);
	virtual ~CMRDCA();
	virtual void cluster(int Kclusters);

	static const int BIG_CONSTANT = std::numeric_limits<int>::max()/100;
	static const long TOTALTIMELIMITSECONDS = 86400;
	static const int TIMELIMIT = 1800;
	bool timeLimitAchieved = false;
	double calcJ(const std::shared_ptr<std::vector<util::CrispCluster> > &clusters) const;
	double calcJ(const util::CrispCluster &cluster) const;
	std::shared_ptr<std::vector<util::CrispCluster> > getClusters() { return clusters; 	 }
	std::shared_ptr<std::vector<util::CrispCluster> > getClustersCopy() const;
	static void seed_random_engine(unsigned seed);

private:
	const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices;
	time_t initialTime = time(NULL);
	static std::default_random_engine generator;
	bool clusterAssign(std::vector<util::CrispCluster> &clusters);
	int isCenterOf(int index);

	double weightedAvgDissim(int element, int medoid, const util::CrispCluster &cluster) const;
	double minimizeRegret(const util::CrispCluster &c, double currentRegret) const;
	double minimizeRegret(const util::CrispCluster &c, int center, double currentRegret) const;
protected:
	const int nElems;
	std::uniform_int_distribution<int> distribution;
	const int nCriteria;
	double maxWeightAbsoluteDifferenceGlobal = 1.0;
	std::vector<int> clusterIndexForElement;
	void bestPrototypes();

	int K;
	std::shared_ptr<std::vector<util::CrispCluster> > clusters;
	int currentIteration;
	double lastJ;
	double epsilon = 1E-6;
	int iterationLimit = 1000;
	void initialize();
	double updateWeights(util::CrispCluster &cluster, double maxValue, int clusterNum);
	bool timeIsUp() const;

};

} /* namespace clustering */
#endif /* CMRDCA_H_ */

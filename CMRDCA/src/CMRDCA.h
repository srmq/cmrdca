/*
 * CMRDCA.h
 *
 *  Created on: May 27, 2013
 *      Author: srmq
 */

#ifndef CMRDCA_H_
#define CMRDCA_H_

#include "DissimMatrix.h"
#include "FuzzyCluster.h"
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

	int getM() const {
		return m;
	}

	void setM(int m = DEFAULT_M) {
		this->m = m;
	}

	static const int DEFAULT_M = 2;
	static const int BIG_CONSTANT = std::numeric_limits<int>::max()/100;
	static const long TOTALTIMELIMITSECONDS = 86400;
	static const int TIMELIMIT = 1800;
	bool timeLimitAchieved = false;
	double calcJ(const std::shared_ptr<std::vector<util::FuzzyCluster> > &clusters) const;
	double calcJ(const util::FuzzyCluster &cluster) const;
	std::shared_ptr<std::vector<util::FuzzyCluster> > getClusters() { return clusters; 	 }
	std::shared_ptr<std::vector<util::FuzzyCluster> > getClustersCopy() const;
	static int getBestClusterIndex(const std::shared_ptr<std::vector<util::FuzzyCluster> >& clusters, int i);
	static void seed_random_engine(unsigned seed);

	bool isPossibilisticMode() const {
		return possibilisticMode;
	}

	void setPossibilisticMode(bool possibilisticMode = false) {
		this->possibilisticMode = possibilisticMode;
	}

private:
	const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices;
	int m = DEFAULT_M;
	bool possibilisticMode = false;
	double oneOverMMinusOne = 1.0/(m - 1.0);
	time_t initialTime = time(NULL);
	static std::default_random_engine generator;


	void updateUik(int i, util::FuzzyCluster &fc, int K);
	double maxRegret(int i, int gk, const util::FuzzyCluster &cluster) const;
	double calcRegret(const util::FuzzyCluster &c, double currentRegret) const;
	double calcRegret(const util::FuzzyCluster &c, int center, double currentRegret) const;
protected:
	const int nElems;
	std::uniform_int_distribution<int> distribution;
	const int nCriteria;
	double maxWeightAbsoluteDifferenceGlobal = 1.0;

	int K;
	std::shared_ptr<std::vector<util::FuzzyCluster> > clusters;
	int currentIteration;
	double lastJ;
	double epsilon = 1E-4;
	int iterationLimit = 1000;
	void initialize();
	void updateMembershipDegrees(util::FuzzyCluster &fc, int K);
	void bestPrototypes();
	double updateWeights(util::FuzzyCluster &cluster, double maxValue, int clusterNum);
	bool timeIsUp() const;
	// void addEquations(int el, double uik, int ck, glp_prob *lp, int pvarIndex);

};

} /* namespace clustering */
#endif /* CMRDCA_H_ */

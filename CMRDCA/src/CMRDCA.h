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
#include <set>

namespace clustering {

struct MedoidDistance {
	int medoidIndex;
	double totalDistance;

	bool operator<(const MedoidDistance& other) const {
		return totalDistance < other.totalDistance;
	}

};


class CMRDCA {
public:
	CMRDCA(const std::vector<std::shared_ptr<util::IDissimMatrix>>& dissimMatrices);
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
	int getIterationsToConverge() const {return this->currentIteration; }
	void setNumbOfMedoids(int number) {this->numbOfMedoids = number;}
	int getNumbOfMedoids() const { return this->numbOfMedoids; }
	int getIterationLimit() const {return this->iterationLimit;}
	void setIterationLimit(int numIter) { this->iterationLimit = numIter; }
	std::shared_ptr<std::set<int>> blackListElementsForMedoids(double percentOfMeanVariance);
	inline void setUseLocalMedoids(bool localMedoids) { this->useLocalMedoids = localMedoids; }
	inline void setBlackListPercentOfMeanVariance(double lessThanPercentage) { this->blackListPercentOfMeanVariance = lessThanPercentage; }

private:
	time_t initialTime = time(NULL);
	static std::default_random_engine generator;
	bool useLocalMedoids = false;
	double blackListPercentOfMeanVariance = 0.3;

	double weightedAvgDissim(int element, const util::CrispCluster &cluster) const;
	double weightedAvgDissim(int element, const util::CrispCluster &cluster, const std::set<int> &medoids) const;
	double updateWeights(util::CrispCluster &cluster, double maxValue, int clusterNum);
	std::pair<typename std::set<int>,  double>
		findBestState(std::pair<typename std::set<int>,  double> &initialState,
				      std::vector<int> &initPos,
				      std::vector<MedoidDistance> &candidatMedoidsDistances,
				      const std::set<std::set<int> > &existingMedoids);

protected:
	const std::vector<std::shared_ptr<util::IDissimMatrix>>& dissimMatrices;
	const int nElems;
	std::uniform_int_distribution<int> distribution;
	const int nCriteria;
	double maxWeightAbsoluteDifferenceGlobal;
	std::vector<int> clusterIndexForElement;
	void bestPrototypes();

	int K;
	std::shared_ptr<std::vector<util::CrispCluster> > clusters;
	int currentIteration;
	double lastJ;
	double epsilon;
	int iterationLimit;
	int numbOfMedoids;
	std::shared_ptr<std::set<int> > blackListedMedoids;
	void initialize();
	bool timeIsUp() const;
	bool clusterAssign(std::vector<util::CrispCluster> &clusters);
	double distanceToMedoids(int elemIndex, const util::IDissimMatrix& dissimMatrix, const std::set<int> &medoidSet) const;
	void computeBestKMedoidSet(util::CrispCluster &cluster, int k, std::set<int> &result, const std::set<std::set<int> > &existingMedoids);
	void computeBestClusterLocalKMedoidSet(util::CrispCluster &cluster, const int k, std::set<int> &result, const std::set<std::set<int> > &existingMedoids);
	double calcJ(const util::CrispCluster &cluster, const std::set<int> &medoids) const;
};

} /* namespace clustering */
#endif /* CMRDCA_H_ */

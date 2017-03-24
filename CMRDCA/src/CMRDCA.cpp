/*
 * CMRDCA.cpp
 *
 *  Created on: May 27, 2013
 *      Author: srmq
 */

#include "CMRDCA.h"
#include "DissimMatrix.h"
#include <vector>
#include <set>
#include "FuzzyCluster.h"
#include <memory>
#include <random>
#include <cmath>
#include <iostream>
#include <sstream>
#include <ctime>
#include <algorithm>

namespace clustering {

std::default_random_engine CMRDCA::generator(1u);

CMRDCA::CMRDCA(const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices) :
		timeLimitAchieved(false),
		dissimMatrices(dissimMatrices),
		m(CMRDCA::DEFAULT_M),
		possibilisticMode(false),
		oneOverMMinusOne(1.0/(m - 1.0)),
		initialTime(time(NULL)),
		nElems(dissimMatrices.front()->length()),
		distribution(0,this->nElems - 1),
		nCriteria(dissimMatrices.size()),
		maxWeightAbsoluteDifferenceGlobal(1.0),
		currentIteration(0),
		lastJ(std::numeric_limits<double>::max()),
		epsilon(1E-4),
		iterationLimit(1000)
		{
	this->K = 0;
}

CMRDCA::~CMRDCA() {
	// TODO Auto-generated destructor stub
}

void CMRDCA::cluster(int Kclusters) {
	this->K = Kclusters;
	this->initialize();
	bool stop;
	do {
		// Step 1: computation of the best prototypes
		bestPrototypes();

		for (int k = 0; k < this->K; k++) {
			util::FuzzyCluster &currentCluster = this->clusters.get()->at(k);
			const double maxValue = calcJ(currentCluster);
			const double regret = updateWeights(currentCluster, maxValue, k);
			if (regret == -1) {
				std::clog << "WARNING: Could not optimize weights for cluster " + k;
			}
		}

		// Step 3: definition of the best partition
		for (int k = 0; k < this->K; k++) {
			util::FuzzyCluster &fuzzyCluster = this->clusters.get()->at(k);
			updateMembershipDegrees(fuzzyCluster, K);
		}

		// Stop criterion
		const double newJ = this->calcJ(this->clusters);
		stop = abs(newJ - this->lastJ) <= this->epsilon
				|| this->currentIteration > this->iterationLimit;
		this->lastJ = newJ;
		if (timeIsUp()) {
			this->timeLimitAchieved = true;
			std::stringstream ssmsg;
			ssmsg << "WARNING: Returning because time limit of ";
			ssmsg << TOTALTIMELIMITSECONDS;
			ssmsg << " seconds was achieved after ";
			ssmsg << this->currentIteration;
			ssmsg << " iterations!";
			std::clog << ssmsg.str() << std::endl;
		}
	} while (!stop && !this->timeLimitAchieved);

}

bool CMRDCA::timeIsUp() const {
	time_t now;
	time(&now);
	const double timeDif = difftime(now, this->initialTime);
	assert(timeDif >= 0);

	return (timeDif > CMRDCA::TOTALTIMELIMITSECONDS);
}

double CMRDCA::updateWeights(util::FuzzyCluster &cluster, double maxValue, int clusterNum) {
	double regret = -1;
	//FIXME todo

	return regret;
}

void CMRDCA::bestPrototypes() {
	this->currentIteration++;
	std::vector<util::FuzzyCluster> * const clustvecpoint = this->clusters.get();
	for (int k = 0; k < K; k++) {
		double bestResultForThisK = BIG_CONSTANT;
		int newGk = 0;
		util::FuzzyCluster &currentCluster = clustvecpoint->at(k);
		for (int h = 0; h < this->nElems; h++) {
			double resultForThisK = 0.0;
			for(int i = 0; i < this->nElems; i++) {
				resultForThisK +=
						pow(currentCluster.getMembershipDegree(i), this->m)
						*this->maxRegret(i, h, currentCluster);
			}
			if (resultForThisK < bestResultForThisK) {
				bestResultForThisK = resultForThisK;
				newGk = h;
			}
		}
		currentCluster.setCenter(newGk);
	}
}

void CMRDCA::initialize() {
	this->timeLimitAchieved = false;
	time(&this->initialTime);

	this->clusters.reset(new std::vector<util::FuzzyCluster>());
	std::vector<util::FuzzyCluster> * const clustvecpoint = this->clusters.get();
	clustvecpoint->reserve(this->K);
	{
		std::set<int> centers;
		for (int i = 0; i < this->K; i++) {
			util::FuzzyCluster fc = util::FuzzyCluster(this->nElems, this->nCriteria);
			int nextCenter;
			do {
				nextCenter = distribution(CMRDCA::generator);
			} while (centers.find(nextCenter) != centers.end());
			centers.insert(nextCenter);
			fc.setCenter(nextCenter);

			clustvecpoint->push_back(fc);
		}
	}
	for (unsigned int clust = 0; clust < clustvecpoint->size(); ++clust) {
		updateMembershipDegrees(clustvecpoint->at(clust), this->K);
	}
	this->lastJ = calcJ(this->clusters);

}

double CMRDCA::calcJ(const std::shared_ptr<std::vector<util::FuzzyCluster> > &clusters) const {
	double J = 0.0;

	std::vector<util::FuzzyCluster> * const clustvecpoint = clusters.get();
	for (unsigned int clust = 0; clust < clustvecpoint->size(); ++clust) {
		J += calcJ(clustvecpoint->at(clust));
	}
	return J;
}

double CMRDCA::calcJ(const util::FuzzyCluster &cluster) const {
	return calcRegret(cluster, BIG_CONSTANT);
}

double CMRDCA::calcRegret(const util::FuzzyCluster &c, double currentRegret) const {
	return calcRegret(c, c.getCenter(), currentRegret);
}

double CMRDCA::calcRegret(const util::FuzzyCluster &c, int center, double currentRegret) const {
	double sumRegret = 0.0;
	for (int el = 0; el < this->nElems; el++) {
		sumRegret += pow(c.getMembershipDegree(el), this->m)*maxRegret(el, center, c);
		if (sumRegret > currentRegret) return BIG_CONSTANT;
	}

	return sumRegret;
}

void CMRDCA::updateMembershipDegrees(util::FuzzyCluster &fc, int K) {
	for (int i = 0; i < this->nElems; ++i) {
		updateUik(i, fc, K);
	}
}

void CMRDCA::updateUik(int i, util::FuzzyCluster &fc, int K) {
	double uik = 0.0;
	const double sumMultiplier = pow(maxRegret(i, fc.getCenter(), fc), oneOverMMinusOne);
	for (int h = 0; h < K; h++) {
		std::vector<util::FuzzyCluster> * const clustvecpoint = this->clusters.get();
		const util::FuzzyCluster &theHcluster= clustvecpoint->at(h);
		const double denominator = pow(maxRegret(i, theHcluster.getCenter(), theHcluster), oneOverMMinusOne);
		uik += sumMultiplier/denominator;
	}
	fc.updateMemberhipDegree(i, 1.0/uik);


}

double CMRDCA::maxRegret(int i, int gk, const util::FuzzyCluster &cluster) const {
	double maxRegret = std::numeric_limits<double>::min();
	double myRegret;
	for (int j = 0; j < nCriteria; ++j) {
		myRegret = this->dissimMatrices[j]->getDissim(i, gk) * cluster.weightOf(j);
		if (myRegret > maxRegret) {
			maxRegret = myRegret;
		}
	}
	return maxRegret;
}


std::shared_ptr<std::vector<util::FuzzyCluster> > CMRDCA::getClustersCopy() const {
	std::shared_ptr<std::vector<util::FuzzyCluster> >
	result(new std::vector<util::FuzzyCluster>(this->clusters->begin(), this->clusters->end()
					));
	return result;
}

int CMRDCA::getBestClusterIndex(
		const std::shared_ptr<std::vector<util::FuzzyCluster> >& clusters,
		int i) {

	double bestMemberShipDegree = BIG_CONSTANT*-1;
	int bestIndex = -1;
	const int clustersSize = (int) clusters->size();
	for (int k = 0; k < clustersSize; k++) {
		util::FuzzyCluster const &fc = clusters->at(k);
		const double membershipDegree = fc.getMembershipDegree(i);
		if (membershipDegree > bestMemberShipDegree) {
			bestMemberShipDegree = membershipDegree;
			bestIndex = k;
		}
	}
	return bestIndex;
}

void CMRDCA::seed_random_engine(unsigned seed) {
	CMRDCA::generator.seed(seed);
}

} /* namespace clustering */


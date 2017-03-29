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
#include "CrispCluster.h"
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
		initialTime(time(NULL)),
		dissimMatrices(dissimMatrices),
		nElems(dissimMatrices.front()->length()),
		distribution(0,this->nElems - 1),
		nCriteria(dissimMatrices.size()),
		maxWeightAbsoluteDifferenceGlobal(1.0),
		clusterIndexForElement(this->nElems, -1),
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
	bool changed;
	do {
		// Step 1: computation of the best prototypes
		bestPrototypes();

		if (this->dissimMatrices.size() > 1) {
			// Step 2: update weights
			for (int k = 0; k < this->K; k++) {
				util::CrispCluster &currentCluster = this->clusters.get()->at(k);
				const double maxValue = calcJ(currentCluster);
				const double regret = updateWeights(currentCluster, maxValue, k);
				if (regret == -1) {
					std::clog << "WARNING: Could not optimize weights for cluster " + k;
				}
			}
		}

		// Step 3: definition of the best partition
		changed = clusterAssign(*(this->clusters.get()));

		// Stop criterion
		const double newJ = this->calcJ(this->clusters);
		stop = !changed || abs(newJ - this->lastJ) <= this->epsilon
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

double CMRDCA::updateWeights(util::CrispCluster &cluster, double maxValue, int clusterNum) {
	int p;
	cluster.getWeights(&p);

	double num = 1.0;
	for (int h = 0; h < p; h++) { //productory
		double sumNum = 0;
		for (std::set<int>::const_iterator elI = cluster.getElements().get()->begin(); elI != cluster.getElements().get()->end(); elI++) {
			sumNum += this->dissimMatrices[h]->getDissim(*elI, cluster.getCenter());
		}
		num *= sumNum;
	}
	double powNum = pow(num, 1.0/p);

	std::shared_ptr<double> newWeightsShrPtr(new double[p], std::default_delete<double[]>());
	double *newWeights = newWeightsShrPtr.get();

	bool allDenNonZero = true;

	for (int j = 0; j < p; j++) {
		double denominator = 0;
		for (std::set<int>::const_iterator elI = cluster.getElements().get()->begin(); elI != cluster.getElements().get()->end(); elI++) {
			denominator += this->dissimMatrices[j]->getDissim(*elI, cluster.getCenter());
		}
		if (denominator > epsilon) {
			newWeights[j] = powNum / denominator;
		} else {
			allDenNonZero = false;
			break;
		}
	}

	double regret;
	if (allDenNonZero) {
		cluster.setWeights(newWeightsShrPtr);
		regret = calcJ(cluster);
	} else {
		regret = -1;
	}


	return regret;
}

void CMRDCA::bestPrototypes() {
	this->currentIteration++;
	std::vector<util::CrispCluster> * const clustvecpoint = this->clusters.get();
	for (int k = 0; k < K; k++) {
		double bestResultForThisK = BIG_CONSTANT;
		int newGk = 0;
		util::CrispCluster &currentCluster = clustvecpoint->at(k);
		for (int candidateMedoid = 0; candidateMedoid < this->nElems; candidateMedoid++) {
			double totDissimForThisMedoid = 0.0;
			for(std::set<int>::const_iterator cit = currentCluster.getElements().get()->begin(); cit != currentCluster.getElements().get()->end(); cit++) {
				totDissimForThisMedoid +=
						this->weightedAvgDissim(*cit, candidateMedoid, currentCluster);
				if (totDissimForThisMedoid >= bestResultForThisK)
					break;
			}
			if (totDissimForThisMedoid < bestResultForThisK) {
				bestResultForThisK = totDissimForThisMedoid;
				newGk = candidateMedoid;
			}
		}
		currentCluster.setCenter(newGk);
	}
}

void CMRDCA::initialize() {
	this->timeLimitAchieved = false;
	time(&this->initialTime);

	this->clusters.reset(new std::vector<util::CrispCluster>());
	std::vector<util::CrispCluster> * const clustvecpoint = this->clusters.get();
	clustvecpoint->reserve(this->K);
	{
		std::set<int> centers;
		for (int i = 0; i < this->K; i++) {
			util::CrispCluster fc = util::CrispCluster(this->nCriteria);
			int nextCenter;
			do {
				nextCenter = distribution(CMRDCA::generator);
			} while (centers.find(nextCenter) != centers.end());
			centers.insert(nextCenter);
			fc.setCenter(nextCenter);

			clustvecpoint->push_back(fc);
		}
	}

	clusterAssign(*(this->clusters.get()));

	this->lastJ = calcJ(this->clusters);

}

double CMRDCA::calcJ(const std::shared_ptr<std::vector<util::CrispCluster> > &clusters) const {
	double J = 0.0;

	std::vector<util::CrispCluster> * const clustvecpoint = clusters.get();
	for (unsigned int clust = 0; clust < clustvecpoint->size(); ++clust) {
		J += calcJ(clustvecpoint->at(clust));
	}
	return J;
}

double CMRDCA::calcJ(const util::CrispCluster &cluster) const {
	double clusterJ = 0.0;
	for (std::set<int>::const_iterator cit = cluster.getElements().get()->begin(); cit != cluster.getElements().get()->end(); cit++) {
		clusterJ += weightedAvgDissim(*cit, cluster.getCenter(), cluster);
	}

	return clusterJ;

}

double CMRDCA::minimizeRegret(const util::CrispCluster &c, double currentRegret) const {
	return minimizeRegret(c, c.getCenter(), currentRegret);
}

double CMRDCA::minimizeRegret(const util::CrispCluster &c, int center, double currentRegret) const {
	double sumRegret = 0.0;
	for (std::set<int>::const_iterator cit = c.getElements().get()->begin(); cit != c.getElements().get()->end(); cit++) {
		sumRegret += weightedAvgDissim(*cit, center, c);
		if (sumRegret > currentRegret) return BIG_CONSTANT;
	}

	return sumRegret;
}

double CMRDCA::weightedAvgDissim(int element, int medoid, const util::CrispCluster &cluster) const {
	double avgDissim = 0;
	double sumWeights = 0;
	for (int j = 0; j < nCriteria; ++j) {
		const double critJWeight = cluster.weightOf(j);
		avgDissim += this->dissimMatrices[j]->getDissim(element, medoid) * critJWeight;
		sumWeights += critJWeight;
	}
	avgDissim /= sumWeights;

	return avgDissim;
}


std::shared_ptr<std::vector<util::CrispCluster> > CMRDCA::getClustersCopy() const {
	std::shared_ptr<std::vector<util::CrispCluster> >
	result(new std::vector<util::CrispCluster>(this->clusters->begin(), this->clusters->end()
					));
	return result;
}

void CMRDCA::seed_random_engine(unsigned seed) {
	CMRDCA::generator.seed(seed);
}

bool CMRDCA::clusterAssign(std::vector<util::CrispCluster> &clusters) {
	bool hasChanged = false;
	for (int i = 0; i < this->nElems; i++) {
		int centerOf;
		if ((centerOf = isCenterOf(i)) == -1) {
			//element i is not the medoid of any cluster
			int clusterIndexMinDist = -1;
			double minDist = this->BIG_CONSTANT;
			if (this->clusterIndexForElement[i] >= 0 && this->clusterIndexForElement[i] < (int)clusters.size()) {
				clusterIndexMinDist = clusterIndexForElement[i];
				minDist = weightedAvgDissim(i, clusters[clusterIndexMinDist].getCenter(), clusters[clusterIndexMinDist]);
			}

			for (size_t k = 0; k < clusters.size(); k++) {
				const double myDist = weightedAvgDissim(i, clusters[k].getCenter(), clusters[k]);
				if ((myDist + epsilon) < minDist) {
					minDist = myDist;
					clusterIndexMinDist = k;
				}
			}

			if (this->clusterIndexForElement[i] != clusterIndexMinDist) {
				if (this->clusterIndexForElement[i] >= 0 && this->clusterIndexForElement[i] < (int)clusters.size()) {
					clusters[this->clusterIndexForElement[i]].remove(i);
				}
				this->clusterIndexForElement[i] = clusterIndexMinDist;
				clusters[this->clusterIndexForElement[i]].insert(i);
				hasChanged = true;
			}
		} else {
			//i is medoid of a cluster
			if (this->clusterIndexForElement[i] != centerOf) {
				if (this->clusterIndexForElement[i] >= 0 && this->clusterIndexForElement[i] < (int)clusters.size()) {
					clusters[this->clusterIndexForElement[i]].remove(i);
				}
				this->clusterIndexForElement[i] = centerOf;
				clusters[centerOf].insert(i);
				hasChanged = true;
			}
		}
	}

	return hasChanged;
}

// Return the index of the first cluster in which the index passed as parameter is the center
int CMRDCA::isCenterOf(int index) {
	std::vector<util::CrispCluster> &cp = *(this->clusters.get());
	for (size_t i = 0; i < cp.size(); i++) {
		if (cp[i].getCenter() == index) {
			return i;
		}
	}
	return -1;
}

} /* namespace clustering */


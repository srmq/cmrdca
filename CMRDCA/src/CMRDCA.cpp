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

CMRDCA::CMRDCA(const std::vector<std::shared_ptr<util::IDissimMatrix>>& dissimMatrices) :
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
		epsilon(1E-6),
		iterationLimit(1000),
		numbOfMedoids(1)
		{
	this->K = 0;
}

CMRDCA::~CMRDCA() {
	// TODO Auto-generated destructor stub
}

void CMRDCA::cluster(int Kclusters) {
	this->K = Kclusters;
	std::cout << "INFO: initializing" << std::endl;
	this->initialize();
	bool stop;
	bool changed;
	do {
		// Step 1: computation of the best prototypes
		std::cout << "INFO: step 1 - computation of the best prototypes" << std::endl;
		bestPrototypes();

		if (this->dissimMatrices.size() > 1) {
			std::cout << "INFO: step 2 - updating weights" << std::endl;
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
		std::cout << "INFO: step 3 - assigning elements to clusters" << std::endl;
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

double CMRDCA::distanceToMedoids(int elemIndex, const util::IDissimMatrix& dissimMatrix, const std::set<int> &medoidSet) const {
	double result = 0;
	for (std::set<int>::const_iterator medoidIt = medoidSet.begin(); medoidIt != medoidSet.end(); medoidIt++) {
		result += dissimMatrix.getDissim(elemIndex, *medoidIt);
	}
	return result;
}


double CMRDCA::updateWeights(util::CrispCluster &cluster, double maxValue, int clusterNum) {
	int p;
	cluster.getWeights(&p);

	double num = 1.0;
	for (int h = 0; h < p; h++) { //productory
		double sumNum = 0;
		for (std::set<int>::const_iterator elI = cluster.getElements().get()->begin(); elI != cluster.getElements().get()->end(); elI++) {
			sumNum += distanceToMedoids(*elI, *(this->dissimMatrices[h].get()), *cluster.getMedoids().get());
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
			denominator += distanceToMedoids(*elI, *(this->dissimMatrices[j].get()), *cluster.getMedoids().get());
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

struct FrontierElement {
	std::pair<typename std::set<int>,  double> state;
	std::vector<int> currPos;
	bool operator<(const FrontierElement& other) const {
		return state.second < other.state.second;
	}
};

static FrontierElement avanceOne(int posIndex, const FrontierElement& currState,
		std::vector<MedoidDistance> &candidatMedoidsDistances) {
	const int oldPos = currState.currPos[posIndex];
	const int newPos = oldPos + 1;
	std::vector<int> newPosVector = currState.currPos;
	newPosVector[posIndex] += 1;

	const double oldElemCost = candidatMedoidsDistances.at(oldPos).totalDistance;
	const double newElemCost = candidatMedoidsDistances.at(newPos).totalDistance;

	double setCost = currState.state.second;
	setCost -= oldElemCost;
	setCost += newElemCost;
	std::set<int> medoids = currState.state.first;
	medoids.erase(candidatMedoidsDistances.at(oldPos).medoidIndex);
	medoids.insert(candidatMedoidsDistances.at(newPos).medoidIndex);

	return FrontierElement( { std::make_pair(medoids, setCost), newPosVector });

}

std::pair<typename std::set<int>,  double>
	CMRDCA::findBestState(std::pair<typename std::set<int>,  double> &initialState,
			      	      std::vector<int> &initPos,
			      	      std::vector<MedoidDistance> &candidatMedoidsDistances,
			      	      const std::set<std::set<int> > &existingMedoids) {

	std::set<FrontierElement> frontier({FrontierElement({initialState, initPos})});
	while(frontier.size() > 0) {
		std::set<FrontierElement>::iterator firstElementIt = frontier.begin();
		FrontierElement elem = *firstElementIt;
		frontier.erase(firstElementIt);
		if (existingMedoids.find(elem.state.first) == existingMedoids.end() ) {
			//is a new medoid set
			return elem.state;
		} else {
			const int lastSetElementAtIndex = elem.currPos[elem.currPos.size() - 1];
			if (lastSetElementAtIndex < (int)(candidatMedoidsDistances.size()-1)) {
				FrontierElement el = avanceOne(elem.currPos.size() - 1, elem, candidatMedoidsDistances);
				frontier.insert(el);
			} else {
				std::cerr << "WARNING: last element in set went to last position" << std::endl;
			}
			for (int i = initPos.size()-2; i >= 0; i--) {
				const int myIndex = elem.currPos[i];
				const int rightNeighborIndex = elem.currPos[i+1];
				assert(rightNeighborIndex > myIndex);
				if (((rightNeighborIndex - myIndex) > 1) && (myIndex < (int)(candidatMedoidsDistances.size()-1))) {
					FrontierElement el = avanceOne(i, elem, candidatMedoidsDistances);
					frontier.insert(el);
				}
			}
		}
	}
	std::cerr << "WARNING: returning empty set for medoids" << std::endl;
	return std::make_pair(std::set<int>(), -1);
}

void CMRDCA::computeBestClusterLocalKMedoidSet(util::CrispCluster &cluster, const int k, std::set<int> &result, const std::set<std::set<int> > &existingMedoids) {
	std::vector<MedoidDistance> medoids;
	medoids.reserve(cluster.getElements()->size());
	for (std::set<int>::const_iterator it = cluster.getElements()->begin(); it != cluster.getElements()->end(); it++) {
		int candidateMedoid = *it;
		if (this->blackListedMedoids->find(candidateMedoid) != this->blackListedMedoids->end())
			continue;
		std::set<int> uniset({candidateMedoid});
		double elemTotalDistance = calcJ(cluster, uniset);
		medoids.push_back(MedoidDistance({candidateMedoid, elemTotalDistance}));
	}
	std::sort(medoids.begin(), medoids.end());
	result.clear();
	if ((int)medoids.size() >= k) {
		double initCost = 0;
		for (int i = 0; i < k; i++) {
			result.insert(medoids[i].medoidIndex);
			initCost += medoids[i].totalDistance;
		}
		double myNewJ = calcJ(cluster, result);
		const double oldJ = calcJ(cluster);
		if (oldJ < myNewJ) {
			result.clear();
			result.insert(cluster.getMedoids()->begin(), cluster.getMedoids()->end());
			std::cerr << "WARNING: unable to find k non blacklisted best medoids inside the cluster. Returning the same set" << std::endl;
		} else {
			//myNewJ is better, now checking if this medoid set already exists elsewhere
			if (existingMedoids.find(result) != existingMedoids.end()) {
				//oh no! repeated medoid set
				result.clear();
				std::pair<typename std::set<int>,  double> initialState = std::make_pair(result, initCost);
				std::vector<int> initPos;
				initPos.reserve(k);
				for (int i = 0; i < k; i++) {
					initPos.push_back(i);
				}
				std::pair<typename std::set<int>,  double> bestState = findBestState(initialState, initPos, medoids, existingMedoids);
				if (!bestState.first.empty()) {
					result.insert(bestState.first.begin(), bestState.first.end());
					myNewJ = calcJ(cluster, result);
					if (oldJ < myNewJ) {
						result.clear();
						result.insert(cluster.getMedoids()->begin(), cluster.getMedoids()->end());
						std::cerr << "WARNING: unable to find k non blacklisted best medoids inside the cluster. Returning the same set" << std::endl;
					}
				} else {
					result.clear();
					result.insert(cluster.getMedoids()->begin(), cluster.getMedoids()->end());
					std::cerr << "WARNING: unable to find k non blacklisted best medoids inside the cluster. Returning the same set" << std::endl;
				}

			}
		}
	} else {
		std::cerr << "WARNING: unable to find k non blacklisted medoids for cluster. Returning the same set" << std::endl;
		result.insert(cluster.getMedoids()->begin(), cluster.getMedoids()->end());
	}
}

void CMRDCA::computeBestKMedoidSet(util::CrispCluster &cluster, const int k, std::set<int> &result, const std::set<std::set<int> > &existingMedoids) {
	std::vector<MedoidDistance> medoids;
	medoids.reserve(this->nElems);
	for (int candidateMedoid = 0; candidateMedoid < this->nElems; candidateMedoid++) {
		if (this->blackListedMedoids->find(candidateMedoid) != this->blackListedMedoids->end())
			continue;
		std::set<int> uniset({candidateMedoid});
		double elemTotalDistance = calcJ(cluster, uniset);
		medoids.push_back(MedoidDistance({candidateMedoid, elemTotalDistance}));
	}
	std::sort(medoids.begin(), medoids.end());
	std::set<int> initialSet;
	double initCost = 0;
	for (int i = 0; i < k; i++) {
		initialSet.insert(medoids[i].medoidIndex);
		initCost += medoids[i].totalDistance;
	}
	result.clear();
	if (existingMedoids.find(initialSet) == existingMedoids.end()) {
		result.insert(initialSet.begin(), initialSet.end());
	} else {
		std::pair<typename std::set<int>,  double> initialState = std::make_pair(initialSet, initCost);
		std::vector<int> initPos;
		initPos.reserve(k);
		for (int i = 0; i < k; i++) {
			initPos.push_back(i);
		}
		std::pair<typename std::set<int>,  double> bestState = findBestState(initialState, initPos, medoids, existingMedoids);
		result.insert(bestState.first.begin(), bestState.first.end());
	}
}

void CMRDCA::bestPrototypes() {
	this->currentIteration++;
	std::set<std::set<int> >  medoidSet;
	for (int k = 0; k < K; k++) {
		std::shared_ptr<std::set<int> > newMedoids(new std::set<int>());
		util::CrispCluster &currentCluster = this->clusters->at(k);
		if (this->useLocalMedoids) {
			computeBestClusterLocalKMedoidSet(currentCluster, this->numbOfMedoids, *(newMedoids.get()), medoidSet);
		} else {
			computeBestKMedoidSet(currentCluster, this->numbOfMedoids, *(newMedoids.get()), medoidSet);
		}
		currentCluster.setMedoids(newMedoids);
		if (newMedoids->size() == 0) {
			std::cerr << "WARNING: empty medoid set" << std::endl;
		}
		medoidSet.insert(*(newMedoids.get()));
	}
}

void CMRDCA::initialize() {
	this->timeLimitAchieved = false;
	time(&this->initialTime);

	this->clusters.reset(new std::vector<util::CrispCluster>());
	std::vector<util::CrispCluster> * const clustvecpoint = this->clusters.get();
	clustvecpoint->reserve(this->K);
	if (this->blackListPercentOfMeanVariance > 0) {
		this->blackListedMedoids = blackListElementsForMedoids(this->blackListPercentOfMeanVariance);
	} else {
		this->blackListedMedoids = std::shared_ptr<std::set<int>>(new std::set<int>());
	}
	{
		std::set<std::set<int>> centers;
		for (int i = 0; i < this->K; i++) {
			util::CrispCluster fc = util::CrispCluster(this->nCriteria);
			std::shared_ptr<std::set<int> > medoidSet(new std::set<int>());
			int nextMedoidCenter;
			do {
				medoidSet->clear();
				while ((int)medoidSet->size() < this->numbOfMedoids) {
					nextMedoidCenter = distribution(CMRDCA::generator);
					if ( (this->blackListedMedoids->find(nextMedoidCenter) != this->blackListedMedoids->end()) || (medoidSet->find(nextMedoidCenter) != medoidSet->end()) )
						continue;
					else {
						medoidSet->insert(nextMedoidCenter);
					}
				}
			} while (centers.find(*(medoidSet.get())) != centers.end());
			centers.insert(*(medoidSet.get()));
			fc.setMedoids(medoidSet);

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
		clusterJ += weightedAvgDissim(*cit, cluster);
	}

	return clusterJ;

}

double CMRDCA::calcJ(const util::CrispCluster &cluster, const std::set<int> &medoids) const {
	double clusterJ = 0.0;
	for (std::set<int>::const_iterator cit = cluster.getElements().get()->begin(); cit != cluster.getElements().get()->end(); cit++) {
		clusterJ += weightedAvgDissim(*cit, cluster, medoids);
	}

	return clusterJ;

}

double CMRDCA::weightedAvgDissim(int element, const util::CrispCluster &cluster) const {
	return weightedAvgDissim(element, cluster, *(cluster.getMedoids().get()));
}

double CMRDCA::weightedAvgDissim(int element, const util::CrispCluster &cluster, const std::set<int> &medoids) const {
	double avgDissim = 0;
	double sumWeights = 0;
	for (int j = 0; j < nCriteria; ++j) {
		const double critJWeight = cluster.weightOf(j);
		avgDissim += distanceToMedoids(element,  *(this->dissimMatrices[j].get()), medoids) * critJWeight;
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
		int clusterIndexMinDist = -1;
		double minDist = this->BIG_CONSTANT;
		if (this->clusterIndexForElement[i] >= 0 && this->clusterIndexForElement[i] < (int)clusters.size()) {
			clusterIndexMinDist = clusterIndexForElement[i];
			minDist = weightedAvgDissim(i, clusters[clusterIndexMinDist]);
		}

		for (size_t k = 0; k < clusters.size(); k++) {
			const double myDist = weightedAvgDissim(i, clusters[k]);
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
	}

	return hasChanged;
}

std::shared_ptr<std::set<int>> CMRDCA::blackListElementsForMedoids(double percentOfMeanVariance) {
	std::shared_ptr<std::set<int>> result(new std::set<int>());

	double varEl[this->nElems];
	double meanVar = 0;

	for (int el = 0; el < this->nElems; el++) {
		double sumSQ = 0;
		double avg = 0;
		for (int elJ = 0; elJ < this->nElems; elJ++) {
		    double dissimToElj = 0;
		    for (int p = 0; p < nCriteria; p++) {
		        dissimToElj += this->dissimMatrices[p]->getDissim(el, elJ)/nCriteria;
		    }
		    sumSQ += dissimToElj*dissimToElj;
		    avg += dissimToElj;
		}
		avg /= this->nElems;
		varEl[el] = sumSQ/this->nElems - avg*avg;
		meanVar += varEl[el];
	}
	meanVar /= this->nElems;

	for (int el = 0; el < this->nElems; el++) {
		if (varEl[el] < percentOfMeanVariance*meanVar) {
			result->insert(el);
		}
	}

	return result;
}


} /* namespace clustering */


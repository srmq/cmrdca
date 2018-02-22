/*
 * CMRDCAGlobal.cpp

 *
 *  Created on: Sep 2, 2013
 *      Author: srmq
 */

#include "CMRDCA.h"
#include "CMRDCAGlobal.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include "CrispCluster.h"

namespace clustering {

CMRDCAGlobal::CMRDCAGlobal(
		const std::vector<std::shared_ptr<util::IDissimMatrix> >& dissimMatrices) :
		CMRDCA(dissimMatrices) {
}

CMRDCAGlobal::~CMRDCAGlobal() {
}

void CMRDCAGlobal::cluster(int Kclusters) {
	this->K = Kclusters;
	std::cout << "INFO: initializing" << std::endl;
	this->initialize();
	{
		std::vector<util::CrispCluster>::iterator clusterIterator = this->clusters.get()->begin();
		std::shared_ptr<double> clusterWeights = (*clusterIterator).getWeights(NULL);
		clusterIterator++;
		while (clusterIterator != this->clusters.get()->end()) {
			(*clusterIterator).setWeights(clusterWeights);
			clusterIterator++;
		}
	}

	bool stop;
	bool changed;
	do {
		std::cout << "INFO: step 1 - computation of the best prototypes" << std::endl;
		// Step 1: computation of the best prototypes
		bestPrototypes();

		const double maxValue = calcJ(this->clusters);

		if (this->dissimMatrices.size() > 1) {
			// Step 2: update weights
			std::cout << "INFO: step 2 - updating weights" << std::endl;
			const double regret = updateWeights(this->clusters, maxValue);
			if (regret == -1) {
				std::clog << "WARNING: Could not optimize weights at iteration " + this->currentIteration;
			}
		}

		std::cout << "INFO: step 3 - assigning elements to clusters" << std::endl;
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

double CMRDCAGlobal::updateWeights(
		std::shared_ptr<std::vector<util::CrispCluster> >& clusters,
		double maxValue) {

	double regret = -1;

	double denominatorByMatrix[nCriteria];

	bool allDenNonZero = true;

	for (int j = 0; j < nCriteria; j++) {
		denominatorByMatrix[j] = 0;
		for (int el = 0; el < nElems; el++) {
			denominatorByMatrix[j] += distanceToMedoids(el, *(this->dissimMatrices[j].get()), *(clusters->at(clusterIndexForElement[el]).getMedoids().get()));
		}
		if (denominatorByMatrix[j] < epsilon) {
			allDenNonZero = false;
		}
	}

	if (allDenNonZero) {
		double numerator = denominatorByMatrix[0];
		for (int h = 1; h < nCriteria; h++) {
			numerator *= denominatorByMatrix[h];
		}
		numerator = pow(numerator, 1.0/(double)nCriteria);
		double *sharedWeights = clusters->at(0).getWeights(NULL).get();
		for (int i = 0; i < nCriteria; i++) {
			sharedWeights[i] = numerator/denominatorByMatrix[i];
		}
		regret = calcJ(clusters);
	}


	return regret;


}

}

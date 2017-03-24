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
#include "FuzzyCluster.h"

namespace clustering {

CMRDCAGlobal::CMRDCAGlobal(
		const std::vector<std::shared_ptr<util::DissimMatrix> >& dissimMatrices) :
		CMRDCA(dissimMatrices) {
}

CMRDCAGlobal::~CMRDCAGlobal() {
}

void CMRDCAGlobal::cluster(int Kclusters) {
	this->K = Kclusters;
	this->initialize();
	{
		std::vector<util::FuzzyCluster>::iterator clusterIterator = this->clusters.get()->begin();
		std::shared_ptr<double> clusterWeights = (*clusterIterator).getWeights(NULL);
		clusterIterator++;
		while (clusterIterator != this->clusters.get()->end()) {
			(*clusterIterator).setWeights(clusterWeights);
			clusterIterator++;
		}
	}

	bool stop;
	do {
		// Step 1: computation of the best prototypes
		bestPrototypes();

		const double maxValue = calcJ(this->clusters);
		const double regret = updateWeights(this->clusters, maxValue);
		if (regret == -1) {
			std::clog << "WARNING: Could not optimize weights at iteration " + this->currentIteration;
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

double CMRDCAGlobal::updateWeights(
		std::shared_ptr<std::vector<util::FuzzyCluster> >& clusters,
		double maxValue) {

	double regret = -1;

	//FIXME TODO


	return regret;


}

}

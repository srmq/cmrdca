/*
 * FWRDCAGlobal.cpp

 *
 *  Created on: Sep 2, 2013
 *      Author: srmq
 */

#include "FWRDCA.h"
#include "FWRDCAGlobal.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include "FuzzyCluster.h"

namespace clustering {

FWRDCAGlobal::FWRDCAGlobal(
		const std::vector<std::shared_ptr<util::DissimMatrix> >& dissimMatrices) :
		FWRDCA(dissimMatrices) {
}

FWRDCAGlobal::~FWRDCAGlobal() {
}

void FWRDCAGlobal::cluster(int Kclusters) {
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

double FWRDCAGlobal::updateWeights(
		std::shared_ptr<std::vector<util::FuzzyCluster> >& clusters,
		double maxValue) {

	double regret = -1;

	glp_prob *lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MIN);

	//pVars
	glp_add_cols(lp, this->nElems);
	for (int i = 1; i <= this->nElems; i++) {
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, BIG_CONSTANT);
		glp_set_obj_coef(lp, i, 1.0);
	}
	//weights
	{
		int i = glp_add_cols(lp, this->nCriteria);
		assert(i == this->nElems + 1);
		for (;i <= (this->nElems + this->nCriteria); i++) {
			glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
			glp_set_obj_coef(lp, i, 0.0);

		}
	}

	for (std::vector<util::FuzzyCluster>::const_iterator it = this->clusters.get()->begin(); it != this->clusters.get()->end(); it++) {
		const util::FuzzyCluster &cluster = (*it);
		for (int i = 0; i < this->nElems; i++) {
			addEquations(i, cluster.getMembershipDegree(i), cluster.getCenter(), lp, (i+1));
		}
	}

	{ // sum of weights = 1.0
		int rowNum = glp_add_rows(lp, 1);
		glp_set_row_bnds(lp, rowNum, GLP_FX, 1.0, 1.0);
		double weightCoefs[this->nCriteria+1];
		std::fill(weightCoefs, weightCoefs + this->nCriteria+1, 1.0);
		int weightIndices[this->nCriteria+1];
		{
			int j = 0;
			for (int i = this->nElems+1; i <= (this->nElems + this->nCriteria); i++, j++) {
				weightIndices[j+1] = i;
			}
		}
		glp_set_mat_row(lp, rowNum, this->nCriteria, weightIndices, weightCoefs);
	}

	{ //maxWeightDifference
		if (this->maxWeightAbsoluteDifferenceGlobal != 1.0) {
			for (int i = this->nElems+1; i < (this->nElems + this->nCriteria); i++) {
				for (int j = i+1; j <= (this->nElems + this->nCriteria); j++) {
					int index = glp_add_cols(lp, 1);
					glp_set_col_bnds(lp, index, GLP_DB, 0.0, this->maxWeightAbsoluteDifferenceGlobal);
					glp_set_obj_coef(lp, index, 0.0);
					int rowNum = glp_add_rows(lp, 3);
					glp_set_row_bnds(lp, rowNum, GLP_DB, 0.0, this->maxWeightAbsoluteDifferenceGlobal);
					{
						double weightCoefs[] = {0, 1.0};
						int weightIndex[] = {0, index};
						glp_set_mat_row(lp, rowNum, 1, weightIndex, weightCoefs);
					}
					{
						rowNum++;
						glp_set_row_bnds(lp, rowNum, GLP_LO, 0.0, 0.0);
						double weightCoefs[] = {0, 1.0, -1.0, 1.0};
						int weightIndices[] = {0, index, i, j};
						glp_set_mat_row(lp, rowNum, 3, weightIndices, weightCoefs);
					}
					{
						rowNum++;
						glp_set_row_bnds(lp, rowNum, GLP_LO, 0.0, 0.0);
						double weightCoefs[] = {0, 1.0, -1.0, 1.0};
						int weightIndices[] = {0, index, j, i};
						glp_set_mat_row(lp, rowNum, 3, weightIndices, weightCoefs);
					}

				}
			}

		}
	} //end maxWeightDifference

	{ // lessthanequal maxvalue
		int rowNum = glp_add_rows(lp, 1);
		glp_set_row_bnds(lp, rowNum, GLP_UP, maxValue, maxValue);
		int elemIndices[this->nElems+1];

		for (int i = 1; i <= this->nElems; i++) {
			elemIndices[i] = i;
		}
		double elemCoefs[this->nElems+1];
		std::fill(elemCoefs, elemCoefs + this->nElems+1, 1.0);
		glp_set_mat_row(lp, rowNum, this->nElems, elemIndices, elemCoefs);
	}

	//glp_simplex(lp, &this->glpParms);
	{
		int lpStatus = GLP_OPT; //glp_get_status(lp);
		if ( lpStatus == GLP_OPT || lpStatus == GLP_FEAS) {
			glp_iocp mip_options;
			glp_init_iocp(&mip_options);
			mip_options.presolve = GLP_ON;
			glp_intopt(lp, &mip_options);
			int mipStatus =  glp_mip_status(lp);
			if (mipStatus == GLP_OPT || mipStatus == GLP_FEAS) {
				regret = glp_mip_obj_val(lp);
				int nWeights;
				std::shared_ptr<double> clusterWeights = this->clusters.get()->at(0).getWeights(&nWeights);
				for (int i = 0; i < nWeights; i++) {
					const double weightValue = glp_mip_col_val(lp, (this->nElems)+1+i);
					clusterWeights.get()[i] = weightValue;
				}
				if (mipStatus == GLP_FEAS) {
					std::stringstream ssmsg;
					ssmsg << "INFO: Using non-optimal values for weights on iteration ";
					ssmsg << this->currentIteration;
					std::clog << ssmsg.str() << std::endl;
				}
			}
		}
		glp_delete_prob(lp);
	}

	return regret;


}

}

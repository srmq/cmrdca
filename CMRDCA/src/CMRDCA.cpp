/*
 * FWRDCA.cpp
 *
 *  Created on: May 27, 2013
 *      Author: srmq
 */

#include "FWRDCA.h"
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
#include <glpk.h>
#include <algorithm>

namespace clustering {

std::default_random_engine FWRDCA::generator(1u);

FWRDCA::FWRDCA(const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices) :
		timeLimitAchieved(false),
		dissimMatrices(dissimMatrices),
		m(FWRDCA::DEFAULT_M),
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

FWRDCA::~FWRDCA() {
	// TODO Auto-generated destructor stub
}

void FWRDCA::cluster(int Kclusters) {
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

bool FWRDCA::timeIsUp() const {
	time_t now;
	time(&now);
	const double timeDif = difftime(now, this->initialTime);
	assert(timeDif >= 0);

	return (timeDif > FWRDCA::TOTALTIMELIMITSECONDS);
}

double FWRDCA::updateWeights(util::FuzzyCluster &cluster, double maxValue, int clusterNum) {
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

	for (int i = 0; i < this->nElems; i++) {
		addEquations(i, cluster.getMembershipDegree(i), cluster.getCenter(), lp, (i+1));
	}


		int rowNum = glp_add_rows(lp, 1);
	if (!this->possibilisticMode){ // sum of weights = 1.0
		glp_set_row_bnds(lp, rowNum, GLP_FX, 1.0, 1.0);
	} else {
		glp_set_row_bnds(lp, rowNum, GLP_LO, 1.0, 1.0);
	}
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

	glp_simplex(lp, &this->glpParms);
	{
		int lpStatus = glp_get_status(lp);
		if ( lpStatus == GLP_OPT || lpStatus == GLP_FEAS) {
			glp_intopt(lp, NULL); //FIXME colocar parametros
			int mipStatus =  glp_mip_status(lp);
			if (mipStatus == GLP_OPT || mipStatus == GLP_FEAS) {
				regret = glp_mip_obj_val(lp);
				int nWeights;
				std::shared_ptr<double> clusterWeights = cluster.getWeights(&nWeights);
				for (int i = 0; i < nWeights; i++) {
					const double weightValue = glp_mip_col_val(lp, (this->nElems)+1+i);
					clusterWeights.get()[i] = weightValue;
				}
				if (mipStatus == GLP_FEAS) {
					std::stringstream ssmsg;
					ssmsg << "INFO: Using non-optimal values on iteration ";
					ssmsg << this->currentIteration;
					ssmsg << " for cluster ";
					ssmsg << clusterNum;
					std::clog << ssmsg.str() << std::endl;
				}
			}
		}
		glp_delete_prob(lp);
	}

	return regret;
}

void FWRDCA::addEquations(int el, double uik, int ck, glp_prob *lp, int myVarIndex) {
	int j = nCriteria - 1;
	int i = j - 1;
	int lmaxIndex = glp_add_cols(lp, 1);
	glp_set_col_bnds(lp, lmaxIndex, GLP_DB, 0.0, BIG_CONSTANT);
	glp_set_obj_coef(lp, lmaxIndex, 0.0);

	int indexWeighti = (this->nElems)+1+i;
	int indexWeightj = (this->nElems)+1+j;

	double coefWeightJ = dissimMatrices[j]->getDissim(el, ck);
	double coefWeightI = dissimMatrices[i]->getDissim(el, ck);

	int rowNum = glp_add_rows(lp, 2);

	{
		glp_set_row_bnds(lp, rowNum, GLP_UP, 0.0, 0.0);
		int indices[] = {0, indexWeightj, lmaxIndex};
		double coefs[] = {0, coefWeightJ, -1.0};
		glp_set_mat_row(lp, rowNum, 2, indices, coefs);
	}

	{
		rowNum++;
		glp_set_row_bnds(lp, rowNum, GLP_UP, 0.0, 0.0);
		int indices[] = {0, indexWeighti, lmaxIndex};
		double coefs[] = {0, coefWeightI, -1.0};
		glp_set_mat_row(lp, rowNum, 2, indices, coefs);
	}

	int cIndex = glp_add_cols(lp, 1);
	glp_set_col_bnds(lp, cIndex, GLP_DB, 0.0, 1.0);
	glp_set_obj_coef(lp, cIndex, 0.0);
	glp_set_col_kind(lp, cIndex, GLP_BV);

	rowNum = glp_add_rows(lp, 2);

	{
		glp_set_row_bnds(lp, rowNum, GLP_UP, 0.0, 0.0);
		int indices[] = {0, lmaxIndex, indexWeightj, cIndex};
		double coefs[] = {0, 1.0, -1.0*dissimMatrices[j]->getDissim(el, ck), -1.0*BIG_CONSTANT};
		glp_set_mat_row(lp, rowNum, 3, indices, coefs);
	}

	{
		rowNum++;
		glp_set_row_bnds(lp, rowNum, GLP_UP, BIG_CONSTANT, BIG_CONSTANT);
		int indices[] = {0, lmaxIndex, indexWeighti, cIndex};
		double coefs[] = {0, 1.0, -1.0*dissimMatrices[i]->getDissim(el, ck), BIG_CONSTANT};
		glp_set_mat_row(lp, rowNum, 3, indices, coefs);

	}

	i--;
	while (i >= 0) {
		indexWeighti = (this->nElems)+1+i;
		int newMaxIndex = glp_add_cols(lp, 1);
		glp_set_col_bnds(lp, newMaxIndex, GLP_DB, 0.0, BIG_CONSTANT);
		glp_set_obj_coef(lp, newMaxIndex, 0.0);

		{
			rowNum = glp_add_rows(lp, 1);
			glp_set_row_bnds(lp, rowNum, GLP_LO, 0.0, 0.0);
			int indices[] = {0, lmaxIndex, newMaxIndex};
			double coefs[] = {0, -1.0, 1.0};
			glp_set_mat_row(lp, rowNum, 2, indices, coefs);
		}

		{
			rowNum = glp_add_rows(lp, 1);
			glp_set_row_bnds(lp, rowNum, GLP_LO, 0.0, 0.0);
			int indices[] = {0, newMaxIndex, indexWeighti};
			double coefs[] = {0, 1.0, -1.0*dissimMatrices[i]->getDissim(el, ck)};
			glp_set_mat_row(lp, rowNum, 2, indices, coefs);
		}

		{ //boolvar
			int dIndex = glp_add_cols(lp, 1);
			glp_set_col_bnds(lp, dIndex, GLP_DB, 0.0, 1.0);
			glp_set_obj_coef(lp, dIndex, 0.0);
			glp_set_col_kind(lp, dIndex, GLP_BV);

			rowNum = glp_add_rows(lp, 1);
			glp_set_row_bnds(lp, rowNum, GLP_UP, 0.0, 0.0);
			{
			int indices[] = {0, newMaxIndex, lmaxIndex, dIndex};
			double coefs[] = {0, 1.0, -1.0, -1.0*BIG_CONSTANT};
			glp_set_mat_row(lp, rowNum, 3, indices, coefs);
			}

			rowNum = glp_add_rows(lp, 1);
			glp_set_row_bnds(lp, rowNum, GLP_UP, BIG_CONSTANT, BIG_CONSTANT);
			int indices[] = {0, newMaxIndex, indexWeighti, dIndex};
			double coefs[] = {0, 1.0, -1.0*dissimMatrices[i]->getDissim(el, ck), BIG_CONSTANT};
			glp_set_mat_row(lp, rowNum, 3, indices, coefs);


		}
		lmaxIndex = newMaxIndex;
		i--;
	} // end while

	{
		rowNum = glp_add_rows(lp, 1);
		glp_set_row_bnds(lp, rowNum, GLP_FX, 0.0, 0.0);
		int indices[] = {0, myVarIndex, lmaxIndex};
		double coefs[] = {0, -1.0, pow(uik, this->m)};
		glp_set_mat_row(lp, rowNum, 2, indices, coefs);
	}

}

void FWRDCA::bestPrototypes() {
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

void FWRDCA::glpkInitialize() {
	glp_init_smcp(&this->glpParms);
	glpParms.tm_lim = FWRDCA::TIMELIMIT*1000;
}

void FWRDCA::initialize() {
	this->timeLimitAchieved = false;
	time(&this->initialTime);
	glpkInitialize();

	this->clusters.reset(new std::vector<util::FuzzyCluster>());
	std::vector<util::FuzzyCluster> * const clustvecpoint = this->clusters.get();
	clustvecpoint->reserve(this->K);
	{
		std::set<int> centers;
		for (int i = 0; i < this->K; i++) {
			util::FuzzyCluster fc = util::FuzzyCluster(this->nElems, this->nCriteria);
			int nextCenter;
			do {
				nextCenter = distribution(FWRDCA::generator);
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

double FWRDCA::calcJ(const std::shared_ptr<std::vector<util::FuzzyCluster> > &clusters) const {
	double J = 0.0;

	std::vector<util::FuzzyCluster> * const clustvecpoint = clusters.get();
	for (unsigned int clust = 0; clust < clustvecpoint->size(); ++clust) {
		J += calcJ(clustvecpoint->at(clust));
	}
	return J;
}

double FWRDCA::calcJ(const util::FuzzyCluster &cluster) const {
	return calcRegret(cluster, BIG_CONSTANT);
}

double FWRDCA::calcRegret(const util::FuzzyCluster &c, double currentRegret) const {
	return calcRegret(c, c.getCenter(), currentRegret);
}

double FWRDCA::calcRegret(const util::FuzzyCluster &c, int center, double currentRegret) const {
	double sumRegret = 0.0;
	for (int el = 0; el < this->nElems; el++) {
		sumRegret += pow(c.getMembershipDegree(el), this->m)*maxRegret(el, center, c);
		if (sumRegret > currentRegret) return BIG_CONSTANT;
	}

	return sumRegret;
}

void FWRDCA::updateMembershipDegrees(util::FuzzyCluster &fc, int K) {
	for (int i = 0; i < this->nElems; ++i) {
		updateUik(i, fc, K);
	}
}

void FWRDCA::updateUik(int i, util::FuzzyCluster &fc, int K) {
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

double FWRDCA::maxRegret(int i, int gk, const util::FuzzyCluster &cluster) const {
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


std::shared_ptr<std::vector<util::FuzzyCluster> > FWRDCA::getClustersCopy() const {
	std::shared_ptr<std::vector<util::FuzzyCluster> >
	result(new std::vector<util::FuzzyCluster>(this->clusters->begin(), this->clusters->end()
					));
	return result;
}

int FWRDCA::getBestClusterIndex(
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

void FWRDCA::seed_random_engine(unsigned seed) {
	FWRDCA::generator.seed(seed);
}

} /* namespace clustering */


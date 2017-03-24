/*
 * FuzzyClusterPartition.h
 *
 *  Created on: Jan 6, 2014
 *      Author: srmq
 */

#ifndef FUZZYINDICES_H_
#define FUZZYINDICES_H_
#include <map>
#include <memory>
#include <vector>
#include <set>
#include <armadillo>
#include <cmath>
#include "FuzzyCluster.h"

namespace util {

class FuzzyIndices {
public:
	FuzzyIndices(const std::shared_ptr<std::vector<util::FuzzyCluster> > &fuzzyClusters);
	void insertPertinenceDegree(std::pair<int, int> elementPartition, double value);
	double getPertinenceDegree(std::pair<int, int> elementPartition) const;
	class AndersonIndices {
		public:
			AndersonIndices(FuzzyIndices &outer);
			double randIndex_4a() const;
			double jmcabIndex_4b() const;
			double hubertIndex_4c() const;
			double wallaceIndex_4d1() const;
			double wallaceIndex_4d2() const;
			double fowlkesMallow_4e() const;
			double jacard_4f() const;
			double adjRandHA_4g() const;
			double adjRandCampello_4h() const;
			double adjRandBrouwer_4i() const;
			double minkowski_4j() const;
			double hGamma_4k() const;
			double yule_4l() const;
			double chisq_4m() const;
			double goodmnKrkl_4n() const;

			double sFRHRIndex() const;

		private:
			FuzzyIndices &outer;
			double indexA, indexB, indexC, indexD;
			double chiSq_4m;
			double gdmnKrkl_4n;
			double sFRHR;
	};

	class CampelloIndices {
		public:
			CampelloIndices(FuzzyIndices &outer);
			double randIndex() const;
			double adjustedRandIndex() const;
			double jaccardCoef() const;
			double fowlkesMallowsIndex() const;
			double minkowskiMeasure() const;

		private:
			FuzzyIndices &outer;
			double a, b, c, d;
			double computeV(const int j1, const int j2) const;
			double computeX(const int j1, const int j2) const;
			double computeY(const int j1, const int j2) const;
			double computeZ(const int j1, const int j2) const;

			/**
			 * a = |V intersection Y|
			 */
			double compute_a(const std::vector<int> &elements) const;

			/**
			 * b = |V intersection Z|
			 */
			double compute_b(const std::vector<int> &elements) const;

			/**
			 * c = |X intersection Y|
			 */
			double compute_c(const std::vector<int> &elements) const;

			/**
			 * d = |X intersection Z|
			 */
			double compute_d(const std::vector<int> &elements) const;



	};
private:
	/**
	 * std::pair is of <elementIndex, partitionIndex>, double is the pertinence degree
	 */
	std::map<std::pair<int, int>, double> aPrioriFuzzyPartition;
	std::set<int> partitionIndices;
	const std::shared_ptr<std::vector<util::FuzzyCluster> >& fuzzyClusters;

	std::shared_ptr<std::vector<int>> elements() const;
	int numberOfClusters() const;
	int numberOfPartitionsReference() const;
};


}

inline void util::FuzzyIndices::insertPertinenceDegree(
		std::pair<int, int> elementPartition, double value) {
	aPrioriFuzzyPartition[elementPartition] = value;
	partitionIndices.insert(elementPartition.second);
}

inline double util::FuzzyIndices::getPertinenceDegree(
		std::pair<int, int> elementPartition) const {
	std::map<std::pair<int, int>, double>::const_iterator it = aPrioriFuzzyPartition.find(elementPartition);
	double result = 0.0;
	if (it != aPrioriFuzzyPartition.end()) {
		result = it->second;
	}
	return result;
}

inline int util::FuzzyIndices::numberOfClusters() const {
	return fuzzyClusters->size();
}

inline int util::FuzzyIndices::numberOfPartitionsReference() const {
	return partitionIndices.size();
}

inline double util::FuzzyIndices::AndersonIndices::randIndex_4a() const {
	const double num = indexA + indexD;
	const double den = indexA + indexB + indexC + indexD;
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::jmcabIndex_4b() const {
	return (indexB + indexC)/(indexA + indexB + indexC + indexD);
}

inline double util::FuzzyIndices::AndersonIndices::hubertIndex_4c() const {
	return ((indexA+indexD)-(indexB+indexC))/(indexA + indexB + indexC + indexD);
}

inline double util::FuzzyIndices::AndersonIndices::wallaceIndex_4d1() const {
	return indexA/(indexA + indexB);
}

inline double util::FuzzyIndices::AndersonIndices::wallaceIndex_4d2() const {
	return indexA/(indexA + indexC);
}

inline double util::FuzzyIndices::AndersonIndices::fowlkesMallow_4e() const {
	return indexA / sqrt((indexA+indexB)*(indexA + indexC));
}

inline double util::FuzzyIndices::AndersonIndices::jacard_4f() const {
	return indexA / (indexA + indexB + indexC);
}

inline double util::FuzzyIndices::AndersonIndices::adjRandHA_4g() const {
	const double num = indexA - ((indexA + indexC)*(indexA + indexB))/(indexA + indexB + indexC + indexD);
	const double den = ((indexA + indexC)+(indexA + indexB))/2.0 - (indexA + indexC)*(indexA + indexB)/(indexA + indexB + indexC + indexD);
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::adjRandCampello_4h() const {
	const double num = indexA - ((indexA + indexC)*(indexA + indexB))/(indexD);
	const double den = ((indexA + indexC)+(indexA + indexB))/2.0 - (indexA + indexC)*(indexA + indexB)/(indexD);
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::adjRandBrouwer_4i() const {
	const double num = 2.0*(indexA*indexD - indexB*indexC);
	const double den = indexC*indexC + indexB*indexB + 2.0*indexA*indexD + (indexA + indexD)*(indexB + indexC);
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::minkowski_4j() const {
	return sqrt((indexB + indexC)/(indexB + indexA));
}

inline double util::FuzzyIndices::AndersonIndices::hGamma_4k() const {
	const int n = outer.elements()->size();
	const double combN2 = 0.5 * n * (n - 1.0);
	const double num = combN2 * indexA - (indexA + indexB)*(indexA + indexC);
	const double den = sqrt((indexA + indexB)*(indexA + indexC)*(combN2 - (indexA + indexB))*(combN2 - (indexA + indexC)));
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::yule_4l() const {
	const double num = (indexA * indexD) - (indexB * indexC);
	const double den = (indexA * indexB) + (indexC * indexD);
	return num/den;
}

inline double util::FuzzyIndices::AndersonIndices::chisq_4m() const {
	return chiSq_4m;
}

inline double util::FuzzyIndices::AndersonIndices::goodmnKrkl_4n() const {
	return gdmnKrkl_4n;
}

inline double util::FuzzyIndices::AndersonIndices::sFRHRIndex() const {
	return sFRHR;
}


#endif /* FUZZYINDICES_H_ */

/*
 * ConfusionMatrix.cpp
 *
 *  Created on: 24/08/2013
 *      Author: srmq
 */

#include "ConfusionMatrix.h"
#include <algorithm>
#include <memory>
#include <utility>
#include <cassert>
#include <cmath>
#include <set>

namespace util {

ConfusionMatrix::~ConfusionMatrix() {
	for (int i = 0; i < this->k; i++) {
		delete[] this->clusterCounts[i];
	}
	delete[] this->clusterCounts;

}

ConfusionMatrix::ConfusionMatrix(int k, int aprioriClassNumber,
		std::shared_ptr<std::vector<std::string> > classNames) :
		ConfusionMatrix(k, aprioriClassNumber) {
	this->classNames = classNames;
}

ConfusionMatrix::ConfusionMatrix(int k, int aprioriClassesNumber) :
		k(k), nClasses(aprioriClassesNumber), bestJ(-1), aprioriClasses(
				aprioriClassesNumber), clusters(k) {
	this->clusterCounts = new int*[k];
	for (int i = 0; i < k; i++) {
		this->clusterCounts[i] = new int[aprioriClassesNumber];
		std::fill(this->clusterCounts[i],
				this->clusterCounts[i] + aprioriClassesNumber, 0);
	}

}

void ConfusionMatrix::putObject(int indexObject, int indexCluster,
		int classIndexObject) {
	clusterCounts[indexCluster][classIndexObject]++;
	ObjectClassif classif;
	classif.indexClassApriori = classIndexObject;
	classif.indexCluster = indexCluster;
	objectClassifs[indexObject] = classif;
	aprioriClasses[classIndexObject].insert(indexObject);
	clusters[indexCluster].insert(indexObject);
}

int ConfusionMatrix::computeTotalOfClass(int classIndex) const {
	int total = 0;
	for (int i = 0; i < k; i++) {
		total += clusterCounts[i][classIndex];
	}
	return total;
}

int ConfusionMatrix::computeTotalOfCluster(int clusterIndex) const {
	int total = 0;
	for (int j = 0; j < nClasses; j++) {
		total += clusterCounts[clusterIndex][j];
	}
	return total;
}

int ConfusionMatrix::computeGrandTotal() const {
	int total = 0;
	for (int i = 0; i < k; i++) {
		total += computeTotalOfCluster(i);
	}
	assert((total == totalByLine()) && (total == (int) objectClassifs.size()));
	return total;
}

int ConfusionMatrix::totalByLine() const {
	int total = 0;
	for (int j = 0; j < nClasses; j++) {
		total += computeTotalOfClass(j);
	}
	return total;
}

double ConfusionMatrix::fMeasureGlobal() const {
	const double n = computeGrandTotal();
	double sum = 0;
	for (int i = 0; i < nClasses; i++) {
		sum += computeTotalOfClass(i) * maxFmeasure(i);
	}
	const double result = sum / n;
	return result;

}

double ConfusionMatrix::maxFmeasure(int partitionI) const {
	double result = -1;
	double currMeasure;
	for (int j = 0; j < k; j++) {
		currMeasure = fMeasure(partitionI, j);
		if (currMeasure > result)
			result = currMeasure;
	}
	return result;

}

int ConfusionMatrix::maxElemsClassInClusterJ(int j) const {
	int max = 0;
	int classIInClusterJ;
	for (int i = 0; i < nClasses; i++) {
		if ((classIInClusterJ = computeClassIClusterJ(i, j)) > max) {
			max = classIInClusterJ;
		}
	}
	return max;
}

double ConfusionMatrix::OERCIndex() const {
	double result = 0;
	for (int j = 0; j < k; j++) {
		result += maxElemsClassInClusterJ(j);
	}
	result /= totalByLine();
	result = 1.0 - result;
	return result;
}

double ConfusionMatrix::CRIndex() const {
	double a, b, c, d;
	a = b = c = d = 0.0;
	char agreeCluster, agreeClass;
	const int nElems = computeGrandTotal();
	for (int i = 0; i < nElems; i++) {
		const ObjectClassif classifI = objectClassifs.at(i);
		const int clusterI = classifI.indexCluster;
		const int classLabelI = classifI.indexClassApriori;
		for (int j = i + 1; j < nElems; j++) {
			const ObjectClassif classifJ = objectClassifs.at(j);
			const int clusterJ = classifJ.indexCluster;
			const int classLabelJ = classifJ.indexClassApriori;
			agreeCluster = (clusterI == clusterJ) ? (char) 1 : (char) 0;
			agreeClass = (classLabelI == classLabelJ) ? (char) 1 : (char) 0;
			a += agreeCluster * agreeClass;
			b += (1 - agreeCluster) * agreeClass;
			c += agreeCluster * (1 - agreeClass);
			d += (1 - agreeCluster) * (1 - agreeClass);
		}
	}
	const double p = a + b + c + d;
	double CR = (a + d) - ((a + b) * (a + c) + (c + d) * (b + d)) / p;
	CR /= p - ((a + b) * (a + c) + (c + d) * (b + d)) / p;
	return CR;
}

double ConfusionMatrix::nMIIndex() const {
	const double n = computeGrandTotal();
	double numerador = 0;
	double den1 = 0;
	double den2 = 0;
	for (int partitionH = 0; partitionH < nClasses; partitionH++) {
		const double nha = aprioriClasses.at(partitionH).size();
		den1 += nha / n * log(nha / n);
		for (int clusterL = 0; clusterL < k; clusterL++) {
			const double nlb = clusters.at(clusterL).size();
			std::set<int> nhlSet = aprioriClasses.at(partitionH);
			const std::set<int> &clust = clusters.at(clusterL);
			for (std::set<int>::iterator it = nhlSet.begin();
					it != nhlSet.end(); it++) {
				if (clust.find(*it) == clust.end()) {
					nhlSet.erase(it);
				}
			}
			const double nhl = nhlSet.size();
			if (nhl != 0) {
				if (nha * nlb != 0) {
					const double addToNumerador = nhl / n
							* log((n * nhl) / (nha * nlb));
					numerador += addToNumerador;
				}
			}
			if (partitionH == 0) {
				if (nlb != 0) {
					den2 += nlb / n * log(nlb / n);
				}
			}
		}
	}
	den1 *= -1.0;
	den2 *= -1.0;
	const double denominador = (den1 + den2) / 2.0;
	return numerador / denominador;

}

void ConfusionMatrix::printMatrix(std::ostream out) const {
	out << "Clusters\tClasses" << std::endl;
	if (classNames.use_count() != 0) {
		assert((int) (classNames.get()->size()) == (nClasses));
		for (std::vector<std::string>::const_iterator it =
				classNames.get()->begin(); it != classNames.get()->end();
				it++) {
			out << "\t" << (*it);
		}
	} else {
		for (int i = 0; i < nClasses; i++) {
			out << "\t" << i;
		}
	}
	out << std::endl;
	for (int i = 0; i < k; i++) {
		out << i;
		for (int c = 0; c < nClasses; c++) {
			out << "\t" << clusterCounts[i][c];
		}
		out << std::endl;
	}
}

} /* namespace util */

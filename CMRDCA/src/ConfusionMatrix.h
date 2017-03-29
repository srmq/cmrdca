/*
 * ConfusionMatrix.h
 *
 *  Created on: 24/08/2013
 *      Author: srmq
 */

#ifndef CONFUSIONMATRIX_H_
#define CONFUSIONMATRIX_H_

#include<string>
#include<vector>
#include<ostream>
#include<iostream>
#include<map>
#include <set>
#include <memory>

namespace util {

class ConfusionMatrix {
private:
	int **clusterCounts;
	int k;
	int nClasses;
	double bestJ;

	typedef struct object_classif {
		int indexCluster;
		int indexClassApriori;
	} ObjectClassif;

	std::map<int, ObjectClassif> objectClassifs;
	std::vector<std::set<int>> aprioriClasses;
	std::vector<std::set<int>> clusters;
	std::shared_ptr<std::vector<std::string> > classNames;




public:
	virtual ~ConfusionMatrix();
	ConfusionMatrix(int k, int aprioriClassesNumber);
	ConfusionMatrix(int k, int aprioriClassNumber, std::shared_ptr<std::vector<std::string> > classNames);

	void putObject(int indexObject, int indexCluster, int classIndexObject);
	int computeTotalOfClass(int classIndex) const;
	int computeTotalOfCluster(int clusterIndex) const;
	int computeClassIClusterJ(int i, int j) const {
		return clusterCounts[j][i];
	}
	int getK() const {
		return k;
	}
	int computeGrandTotal() const;
	int totalByLine() const;
	double getBestJ() const {
		return bestJ;
	}
	void setBestJ(double bestJ) {
		this->bestJ = bestJ;
	}
	double calcPrecision(int partitionI, int clusterJ) const {
		const double result = ((double)computeClassIClusterJ(partitionI, clusterJ))/computeTotalOfCluster(clusterJ);
		return result;
	}
	double calcRecall(int partitionI, int clusterJ) const {
		const double result = ((double)computeClassIClusterJ(partitionI, clusterJ))/computeTotalOfClass(partitionI);
		return result;

	}
	double fMeasure(int partitionI, int clusterJ) const {
		const double result = 2.0*(calcPrecision(partitionI, clusterJ)*calcRecall(partitionI, clusterJ))/(calcPrecision(partitionI, clusterJ)+calcRecall(partitionI, clusterJ));
		return result;
	}
	double fMeasureGlobal() const;
	double maxFmeasure(int partitionI) const;
	int maxElemsClassInClusterJ(int j) const;
	double OERCIndex() const;
	double CRIndex() const;
	double nMIIndex() const;
	void printMatrix(std::ostream &out) const;

};

} // namespace util

#endif /* CONFUSIONMATRIX_H_ */

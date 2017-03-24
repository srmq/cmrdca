/*
 * FWRDCAWine.cpp
 *
 *  Created on: 13/08/2013
 *      Author: srmq
 */
#include "DissimMatrix.h"
#include <string>
#include <memory>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <limits>
#include <iostream>
#include <sstream>
#include <thread>
#include <unistd.h>
#include "FuzzyCluster.h"
#include "FWRDCA.h"
#include "ConfusionMatrix.h"
#include "FuzzyIndices.h"

// Wine dataset has 178 elements
static	const int NELEM = 178;

static std::shared_ptr<util::DissimMatrix> parseFile(const std::string& fileName)  {
	std::shared_ptr<util::DissimMatrix> result = std::shared_ptr<util::DissimMatrix>(new util::DissimMatrix(NELEM));

	std::ifstream f(fileName);

	std::string line;
	const char* delimiters = ", ()";
	getline(f, line);
	while (line.find("DIST_MATRIX") == std::string::npos) {
		getline(f, line);
	}
	for (int i = 0; i < NELEM; i++) {
		getline(f, line);
		char *cstr = new char[line.length() + 1];
		strcpy(cstr, line.c_str());
		char *token = strtok(cstr, delimiters);
		for (int j = 0; j <= i; j++) {
			const double dissimValue = atof(token);

			result->putDissim(i, j, dissimValue);
			token = strtok(NULL, delimiters);
		}
		delete[] cstr;
	}
	f.close();
	return result;
}

void printIndices(int k, const int APRIORICLASSES,
		const std::shared_ptr<std::vector<util::FuzzyCluster> >& bestClusters) {
	std::vector<std::string>* classNamesVector = new std::vector<std::string>();
	classNamesVector->push_back("class 1");
	classNamesVector->push_back("class 2");
	classNamesVector->push_back("class 3");
	std::shared_ptr<std::vector<std::string> > classNames(classNamesVector);
	util::ConfusionMatrix confusionMatrix(k, APRIORICLASSES, classNames);
	util::FuzzyIndices fuzzyIndices(bestClusters);
	for (int i = 0; i < NELEM; i++) {
		int classlabel;
		if (i <= 58)
			classlabel = 0;
		else if (i <= 129)
			classlabel = 1;
		else
			classlabel = 2;

		confusionMatrix.putObject(i,
				clustering::FWRDCA::getBestClusterIndex(bestClusters, i),
				classlabel);
		fuzzyIndices.insertPertinenceDegree(std::make_pair(i, classlabel), 1.0);
	}
	std::cout << ">>>>>>>>>>>> The F-Measure is: ";
	std::cout << confusionMatrix.fMeasureGlobal() << std::endl;
	std::cout << ">>>>>>>>>>>> The CR-Index  is: ";
	std::cout << confusionMatrix.CRIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> OERC Index    is: ";
	std::cout << confusionMatrix.OERCIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> NMI  Index    is: ";
	std::cout << confusionMatrix.nMIIndex() << std::endl;
	util::FuzzyIndices::CampelloIndices campIndices(fuzzyIndices);
	std::cout << ">>>>>>>>>>>> Campello's Rand Index            is: ";
	std::cout << campIndices.randIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> Campello's Adjusted Rand Index   is: ";
	std::cout << campIndices.adjustedRandIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> Campello's Jaccard Coef          is: ";
	std::cout << campIndices.jaccardCoef() << std::endl;
	std::cout << ">>>>>>>>>>>> Campello's Fowlkes-Mallows Index is: ";
	std::cout << campIndices.fowlkesMallowsIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> Campello's Minkowski Measure     is: ";
	std::cout << campIndices.minkowskiMeasure() << std::endl;
	util::FuzzyIndices::AndersonIndices andIndices(fuzzyIndices);
	std::cout << ">>>>>>>>>>>> (Anderson)Fuz Rand Index Rifqi   is: ";
	std::cout << andIndices.sFRHRIndex() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Rand Index (4a)       is: ";
	std::cout << andIndices.randIndex_4a() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Johnson Index (4b)    is: ";
	std::cout << andIndices.jmcabIndex_4b() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Hubert  Index (4c)    is: ";
	std::cout << andIndices.hubertIndex_4c() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Wallace Index (4d1)   is: ";
	std::cout << andIndices.wallaceIndex_4d1() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Wallace Index (4d2)   is: ";
	std::cout << andIndices.wallaceIndex_4d2() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's FowlkesM Index (4e)   is: ";
	std::cout << andIndices.fowlkesMallow_4e() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Jacard Index (4f)     is: ";
	std::cout << andIndices.jacard_4f() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Wallace Index (4g)    is: ";
	std::cout << andIndices.adjRandHA_4g() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Adj Rand Campello(4h) is: ";
	std::cout << andIndices.adjRandCampello_4h() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Adj Rand Brower  (4i) is: ";
	std::cout << andIndices.adjRandBrouwer_4i() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Minkowski        (4j) is: ";
	std::cout << andIndices.minkowski_4j() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Hubert's Gamma   (4k) is: ";
	std::cout << andIndices.hGamma_4k() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Yule             (4l) is: ";
	std::cout << andIndices.yule_4l() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Chi-Squared      (4m) is: ";
	std::cout << andIndices.chisq_4m() << std::endl;
	std::cout << ">>>>>>>>>>>> Anderson's Goodman-Kruskal  (4n) is: ";
	std::cout << andIndices.goodmnKrkl_4n() << std::endl;
	std::cout << std::flush;
}

int main() {
	std::shared_ptr<util::DissimMatrix> tab2 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-2-N.sds");
	std::shared_ptr<util::DissimMatrix> tab3 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-3-N.sds");
	std::shared_ptr<util::DissimMatrix> tab4 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-4-N.sds");
	std::shared_ptr<util::DissimMatrix> tab5 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-5-N.sds");
	std::shared_ptr<util::DissimMatrix> tab6 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-6-N.sds");
	std::shared_ptr<util::DissimMatrix> tab7 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-7-N.sds");
	std::shared_ptr<util::DissimMatrix> tab8 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-8-N.sds");
	std::shared_ptr<util::DissimMatrix> tab9 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-9-N.sds");
	std::shared_ptr<util::DissimMatrix> tab10 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-10-N.sds");
	std::shared_ptr<util::DissimMatrix> tab11 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-11-N.sds");
	std::shared_ptr<util::DissimMatrix> tab12 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-12-N.sds");
	std::shared_ptr<util::DissimMatrix> tab13 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-13-N.sds");
	std::shared_ptr<util::DissimMatrix> tab14 = parseFile("/home/srmq/Dropbox/CIn/research/inria/dados/Wine/wine-L2-14-N.sds");
	std::vector<std::shared_ptr<util::DissimMatrix>> dissimMatrices = {tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10, tab11, tab12, tab13, tab14};
	const int NUMBER_OF_RUNS = 100;
	const int APRIORICLASSES = 3;
	int k = APRIORICLASSES;

	double bestJ = std::numeric_limits<double>::max();
	std::shared_ptr<std::vector<util::FuzzyCluster> > bestClusters;

	const unsigned int procCount = std::thread::hardware_concurrency();

	if (procCount > 1) {
		unsigned int procId = 0;
		bool quitwhile = false;
		int pipeParentChild[procCount][2]; // PARENT WRITES to CHILD, CHILD READS from PARENT
		int pipeChildParent[procCount][2]; // PARENT READS from CHILD, CHILD WRITES to PARENT
		pid_t cpid;
		do {
			if (pipe(pipeParentChild[procId]) == -1) {
				perror("pipe");
			    exit(EXIT_FAILURE);
			}
			if (pipe(pipeChildParent[procId]) == -1) {
				perror("pipe");
			    exit(EXIT_FAILURE);
			}

			cpid = fork();
			if (cpid == -1) {
			    perror("fork");
			    exit(EXIT_FAILURE);
			}

			if (cpid == 0) {    /* I am the child */
				// close WRITE in pipeParentChild
				close(pipeParentChild[procId][1]);
				// close READ in pipeChildParent
				close(pipeChildParent[procId][0]);
				quitwhile = true;
			} else {			/* I am the master */
				// close READ in pipeParentChild
				close(pipeParentChild[procId][0]);
				// close WRITE in pipeChildParent
				close(pipeChildParent[procId][1]);
				procId++;
			}
		} while (procId < procCount && !quitwhile);
		if (cpid == 0) {    /* I am the child do stuff */
			clustering::FWRDCA::seed_random_engine(2u*procId + 1u);
			for (int i = procId; i < NUMBER_OF_RUNS; i=i+procCount) {
				std::cout << "Run number ";
				std::cout << i;
				std::cout << std::endl;
				clustering::FWRDCA clust(dissimMatrices);
				clust.cluster(k);
				std::shared_ptr<std::vector<util::FuzzyCluster> > const myClusters = clust.getClusters();
				const double myJ = clust.calcJ(myClusters);
				std::cout << "J: ";
				std::cout << myJ;
				std::cout << std::endl;
				if (myJ < bestJ) {
					bestJ = myJ;
					bestClusters = clust.getClustersCopy();
				}
			}
			write(pipeChildParent[procId][1], &bestJ, sizeof(double));
			close(pipeChildParent[procId][1]);
			bool amITheBest;
			read(pipeParentChild[procId][0], &amITheBest, sizeof(bool)); // le resultado
			if (amITheBest) {
				printIndices(k, APRIORICLASSES, bestClusters);
			}
		} else {			/* I am the master get results */
			double overallBestJ;
			read(pipeChildParent[0][0], &overallBestJ, sizeof(double)); // le resultado
			unsigned int bestIndex = 0;
			for (unsigned int i = 1; i < procCount; i++) {
				double procBestJ;
				read(pipeChildParent[i][0], &procBestJ, sizeof(double)); // le resultado
				if (procBestJ > overallBestJ) {
					overallBestJ = procBestJ;
					bestIndex = i;
				}
			}
			const bool falseConst = false;
			const bool trueConst = true;
			for (unsigned int i = 0; i < procCount; i++) {
				if (i != bestIndex) {
					write(pipeParentChild[i][1], &falseConst, sizeof(bool));
					close(pipeParentChild[i][1]);
				} else {
					write(pipeParentChild[i][1], &trueConst, sizeof(bool));
					close(pipeParentChild[i][1]);
				}
			}
		}
	} else {
		for (int i = 0; i < NUMBER_OF_RUNS; i++) {
			std::cout << "Run number ";
			std::cout << i;
			std::cout << std::endl;
			clustering::FWRDCA clust(dissimMatrices);
			clust.cluster(k);
			std::shared_ptr<std::vector<util::FuzzyCluster> > const myClusters = clust.getClusters();
			const double myJ = clust.calcJ(myClusters);
			std::cout << "J: ";
			std::cout << myJ;
			std::cout << std::endl;
			if (myJ < bestJ) {
				bestJ = myJ;
				bestClusters = clust.getClustersCopy();
			}
		}
		printIndices(k, APRIORICLASSES, bestClusters);
	}
	return(0);

}

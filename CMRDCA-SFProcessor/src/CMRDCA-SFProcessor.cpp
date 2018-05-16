/*
 * FWRDCAIris.cpp
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
#include "CrispCluster.h"
#include "CMRDCA.h"
#include "CMRDCAGlobal.h"
#include "ConfusionMatrix.h"
#include <utility>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <unordered_map>
#include <cassert>
#include <omp.h>
#include "CMRDCA-SFProcessor.h"

static std::vector<std::string> inputFiles;
static int k;
static int numInicializacao;
static int numIteracoes = -1;
static int numMedoids = 1;
static std::string outputFile;
static int n;
static int numPrioriClusters;
static std::vector<std::shared_ptr<util::IDissimMatrix>> dissimMatrices;
static unsigned int procCount;
static bool useDissimFloats = false;
static bool useLocalMedoids = false;
static double blackListMedoidsPercentVar = 0.3;

static std::pair<std::shared_ptr<util::IDissimMatrix>, std::shared_ptr<std::vector<std::string> > >
	parseDissimMatrix(const std::string& fileName) {

	std::ifstream f(fileName);
	std::shared_ptr<std::vector<std::string> > objNamesPtr = std::shared_ptr<std::vector<std::string> >(new std::vector<std::string>());
	objNamesPtr->reserve(n);

	for (int i = 0; i < n; i++) {
		std::string line;
		getline(f, line);
		size_t pos = line.find_first_of(',');
		objNamesPtr->push_back(line.substr(pos+1, std::string::npos));
	}

	std::shared_ptr<util::IDissimMatrix> dissimMatrixPtr =
			std::shared_ptr<util::IDissimMatrix>(useDissimFloats ? static_cast<util::IDissimMatrix*>(new util::DissimMatrixFloat(n)) : static_cast<util::IDissimMatrix*>(new util::DissimMatrix(n)));

	const char* delimiters = ",\t\r\n\f";
	for (int i = 0; i < n; i++) {
		std::string line;
		getline(f, line);
		char cstr[line.length() + 1];
		strcpy(cstr, line.c_str());
		char *token = strtok(cstr, delimiters);
		for (int j = 0; j <= i; j++) {
			if (token == NULL)
				std::cerr << "token Should not be NULL";
			double dissimValue;
			std::string s(token);
			std::istringstream os(s);
			os >> dissimValue;

			dissimMatrixPtr->putDissim(i, j, dissimValue);
			token = strtok(NULL, delimiters);
		}
	}
	f.close();

	return std::make_pair(dissimMatrixPtr, objNamesPtr);

}

NameSort::NameSort(std::unordered_map<std::string*, int> &refOrder) : refOrder(refOrder) { }

bool NameSort::operator() (const std::pair<std::string*, int> &a, const std::pair<std::string*, int> &b) {
		return this->refOrder[a.first] < refOrder[b.first];
}


static std::shared_ptr<util::IDissimMatrix>
	reorderedDissim(std::unordered_map<std::string*, int> &refOrder,
			        std::vector<std::string> &currOrder,
			        const util::IDissimMatrix &dissim) {

	std::vector<std::pair<std::string*, int> > names;
	names.reserve(currOrder.size());
	for (int i = 0; i < (int)currOrder.size(); i++) {
		std::string *strP = &currOrder[i];
		names.push_back(std::make_pair(strP, i));
	}
	std::sort(names.begin(), names.end(), NameSort(refOrder));
	std::shared_ptr<util::IDissimMatrix> dissimMatrixPtr =
			std::shared_ptr<util::IDissimMatrix>(useDissimFloats ? static_cast<util::IDissimMatrix*>(new util::DissimMatrixFloat(n)) : static_cast<util::IDissimMatrix*>(new util::DissimMatrix(n)));
	for (int i =0; i < n; i++) {
		const int trI = names[i].second;
		for (int j = 0; j <= i; j++) {
			const int trJ = names[j].second;
			dissimMatrixPtr->putDissim(i, j, dissim.getDissim(trI, trJ));
		}
	}

	return dissimMatrixPtr;
}

static void parseDissimMatrices() {
	dissimMatrices.reserve(inputFiles.size());
	bool firstFile = true;
	std::unordered_map<std::string*, int> refOrder;
	std::shared_ptr<std::vector<std::string> > refNamesPtr;
	for (std::vector<std::string>::const_iterator fNameIter = inputFiles.begin(); fNameIter != inputFiles.end(); fNameIter++) {
		std::pair<std::shared_ptr<util::IDissimMatrix>, std::shared_ptr<std::vector<std::string> > >
			parseResult = parseDissimMatrix(*fNameIter);
		if (firstFile) {
			refNamesPtr = parseResult.second;
			refOrder.reserve(refNamesPtr->size());
			for (int i = 0; i < (int)refNamesPtr->size(); i++) {
				refOrder[&(refNamesPtr->at(i))] = i;
			}
			firstFile = false;
		} else {
			assert(refNamesPtr->size() == parseResult.second->size());
			if (!std::equal(refNamesPtr->begin(), refNamesPtr->end(), parseResult.second->begin())) {
				parseResult.first = reorderedDissim(refOrder, *(parseResult.second.get()), *(parseResult.first.get()));
			}
		}
		dissimMatrices.push_back(parseResult.first);
	}
}


void printIndices(int k,
		const std::shared_ptr<std::vector<util::CrispCluster> >& bestClusters,
		std::vector<int>& classLabels, std::vector<std::string>& objectNames, double bestJ, std::ostream &out) {

	util::ConfusionMatrix confusionMatrix(k, numPrioriClusters);

	out << "------CLUSTER CONTENTS-------" << std::endl;
	for (int i = 0; i < (int)bestClusters->size(); i++) {
			const util::CrispCluster& cluster = bestClusters->at(i);
			out << "Cluster " << i << std::endl;
			out << "Medoids: "  << std::endl;
			for (std::set<int>::const_iterator cIt = cluster.getMedoids()->begin(); cIt != cluster.getMedoids()->end(); cIt++) {
				out << objectNames[*cIt] << std::endl;
			}
			out << "-- BEGIN MEMBERS --" << std::endl;
			const std::set<int> &elements = *(cluster.getElements().get());
			for (std::set<int>::const_iterator elIt = elements.begin(); elIt != elements.end(); elIt++) {
				const int classlabel = classLabels[*elIt];
				confusionMatrix.putObject(*elIt, i, classlabel);
				out << objectNames[*elIt] << std::endl;
			}
			out << "-- END MEMBERS --" << std::endl;
	}

	out <<  "------CONFUSION MATRIX-------" << std::endl;
	confusionMatrix.printMatrix(out);
	out <<  "-----------------------------" << std::endl;
	out <<  "Cluster weights: " << std::endl;
	for (int i = 0; i < (int) bestClusters->size(); i++) {
		out <<  i << ": ";
		int p;
		double *weights = bestClusters->at(i).getWeights(&p).get();
		out << "[" << weights[0];
		for (int j = 1; j < p; j++) {
			out << ", " << weights[j];
		}
		out << "]" << std::endl;
	}

	out << ">>>>>>>>>>>> Best     J    is: ";
	out << bestJ << std::endl;
	out << ">>>>>>>>>>>> The F-Measure is: ";
	out << confusionMatrix.fMeasureGlobal() << std::endl;
	out << ">>>>>>>>>>>> The CR-Index  is: ";
	out << confusionMatrix.CRIndex() << std::endl;
	out << ">>>>>>>>>>>> OERC Index    is: ";
	out << confusionMatrix.OERCIndex() << std::endl;
	out << ">>>>>>>>>>>> NMI  Index    is: ";
	out << confusionMatrix.nMIIndex() << std::endl;
	out.flush();
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

static void readConfigFile(char configFileName[]) {
	std::ifstream f(configFileName);

	if (!f.good()) {
		std::cerr << "Error while trying to open file \"" << configFileName << "\"";
		exit(-2);
	}
	std::string line;
	while (getline(f, line)) {
		rtrim(line);
		if (line.length() == 0)
			continue;
		if (line.find("(numCluster)") != std::string::npos) {
			getline(f, line);
			k = atoi(line.c_str());
		} else if (line.find("(numInicializacao)") != std::string::npos) {
			getline(f, line);
			numInicializacao = atoi(line.c_str());
		} else if (line.find("(numIteracoes)") != std::string::npos) {
			getline(f, line);
			numIteracoes = atoi(line.c_str());
		} else if (line.find("(input)") != std::string::npos) {
			bool done = false;
			while (!done) {
				getline(f, line);
				rtrim(line);
				if (line.length() > 0) {
					inputFiles.push_back(line);
				} else {
					done = true;
				}
			}
		} else if (line.find("(output)") != std::string::npos) {
			getline(f, outputFile);
		} else if (line.find("(numIndividuos)") != std::string::npos) {
			getline(f, line);
			n = atoi(line.c_str());
		} else if (line.find("(numPrioriClusters)") != std::string::npos) {
			getline(f, line);
			numPrioriClusters = atoi(line.c_str());
		} else if (line.find("(numMedoids)") != std::string::npos) {
			getline(f, line);
			numMedoids = atoi(line.c_str());
		} else if (line.find("(numThreads)") != std::string::npos) {
			getline(f, line);
			procCount = atoi(line.c_str());
		} else if (line.find("(useDissimFloats)") != std::string::npos) {
			getline(f, line);
			useDissimFloats = (atoi(line.c_str()) != 0);
		} else if (line.find("(useLocalMedoids)") != std::string::npos) {
			getline(f, line);
			useLocalMedoids = (atoi(line.c_str()) != 0);
		} else if (line.find("(blackListPercentOfMeanVar)") != std::string::npos) {
			getline(f, line);
			blackListMedoidsPercentVar = atof(line.c_str());
		}
	}

}

static std::shared_ptr<std::pair<std::vector<int>, std::vector<std::string> > >
	classLabelsForObjects() {

	std::ifstream f(inputFiles[0]);

	if (!f.good()) {
		std::cerr << "Error while trying to open file \"" << inputFiles[0] << "\"";
	}

	std::shared_ptr<std::pair<typename std::vector<int>,
	                typename std::vector<std::string> > >
	                result(new std::pair<typename std::vector<int>,
	                	   typename std::vector<std::string> >(std::vector<int>(), std::vector<std::string>()));

	result->first.reserve(n);
	result->second.reserve(n);

	for (int i = 0; i < n; i++) {
		std::string line;
		getline(f, line);
		const size_t commaPos = line.find_first_of(',');
		const std::string token = line.substr(0, commaPos);
		result->first.push_back(std::stoi(token));
		const std::string name = line.substr(commaPos+1, line.size());
		result->second.push_back(name);
	}
	return result;
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Should give name of config file as argument" << std::endl;
		exit(-1);
	}
	procCount = std::thread::hardware_concurrency();
	readConfigFile(argv[1]);
	omp_set_num_threads(procCount);
	parseDissimMatrices();

	std::shared_ptr<std::pair<std::vector<int>, std::vector<std::string> > > classLabelsAndNames = classLabelsForObjects();
	//std::vector<int> &classLabels = classLabelsAndNames->first;
	//std::vector<std::string> &objNames = classLabelsAndNames->second;
	double bestJ = std::numeric_limits<double>::max();

	std::shared_ptr<std::vector<util::CrispCluster> > bestClusters;


	std::ofstream outFile(outputFile);
	if (!outFile.is_open()) {
		std::cerr << "Unable to open output file: " << outputFile << std::endl;
	}

	#pragma omp parallel for
	for (int i = 0; i < numInicializacao; i++) {
		#pragma omp critical
		{
			outFile << "Run number ";
			outFile << i;
			outFile << std::endl;
		}
		clustering::CMRDCA *clusteringAlgo;
		if (argc == 2) {
			clusteringAlgo = new clustering::CMRDCA(dissimMatrices);
		} else {
			clusteringAlgo = new clustering::CMRDCAGlobal(dissimMatrices);
			#pragma omp critical
			outFile << "RUNNING GLOBAL" << std::endl;
		}
		if (numIteracoes > 0) {
			clusteringAlgo->setIterationLimit(numIteracoes);
		}
		clusteringAlgo->setNumbOfMedoids(numMedoids);
		clusteringAlgo->setUseLocalMedoids(useLocalMedoids);
		clusteringAlgo->setBlackListPercentOfMeanVariance(blackListMedoidsPercentVar);
		clusteringAlgo->cluster(k);
		std::shared_ptr<std::vector<util::CrispCluster> > const myClusters = clusteringAlgo->getClusters();
		const double myJ = clusteringAlgo->calcJ(myClusters);
		#pragma omp critical
		{
			outFile << "J: ";
			outFile << myJ;
			outFile << std::endl;
			if (myJ < bestJ) {
				bestJ = myJ;
				bestClusters = clusteringAlgo->getClustersCopy();
			}
		}
		delete(clusteringAlgo);
	}
	printIndices(k, bestClusters, classLabelsAndNames->first, classLabelsAndNames->second, bestJ, outFile);

	outFile.flush();
	outFile.close();
	return(0);
}

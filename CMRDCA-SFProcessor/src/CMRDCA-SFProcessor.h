/*
 * CMRDCA-SFProcessor.h
 *
 *  Created on: Mar 28, 2017
 *      Author: srmq
 */

#ifndef CMRDCA_SFPROCESSOR_H_
#define CMRDCA_SFPROCESSOR_H_

#include<unordered_map>
#include <utility>

class NameSort {
public:
	NameSort(std::unordered_map<std::string*, int> &refOrder);
	bool operator() (std::pair<std::string*, int> &a, std::pair<std::string*, int> &b);

private:
	std::unordered_map<std::string*, int> & refOrder;
};



#endif /* CMRDCA_SFPROCESSOR_H_ */

/*
 * FWRDCAGlobal.h
 *
 *  Created on: Sep 2, 2013
 *      Author: srmq
 */

#ifndef FWRDCAGLOBAL_H_
#define FWRDCAGLOBAL_H_

#include "FWRDCA.h"

namespace clustering {

class FWRDCAGlobal : public FWRDCA {
public:
	FWRDCAGlobal(const std::vector<std::shared_ptr<util::DissimMatrix>>& dissimMatrices);
	virtual ~FWRDCAGlobal();
	virtual void cluster(int Kclusters);

private:
	double updateWeights(std::shared_ptr<std::vector<util::FuzzyCluster> > &clusters, double maxValue);
};

}

#endif /* FWRDCAGLOBAL_H_ */

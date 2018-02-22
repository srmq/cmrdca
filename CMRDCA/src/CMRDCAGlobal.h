/*
 * CMRDCAGlobal.h
 *
 *  Created on: Sep 2, 2013
 *      Author: srmq
 */

#ifndef CMRDCAGLOBAL_H_
#define CMRDCAGLOBAL_H_

#include "CMRDCA.h"

namespace clustering {

class CMRDCAGlobal : public CMRDCA {
public:
	CMRDCAGlobal(const std::vector<std::shared_ptr<util::IDissimMatrix>>& dissimMatrices);
	virtual ~CMRDCAGlobal();
	virtual void cluster(int Kclusters);

private:
	double updateWeights(std::shared_ptr<std::vector<util::CrispCluster> > &clusters, double maxValue);
};

}

#endif /* CMRDCAGLOBAL_H_ */

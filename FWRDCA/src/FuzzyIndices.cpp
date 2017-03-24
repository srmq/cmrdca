/*
 * FuzzyIndices.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: srmq
 */

#include "FuzzyIndices.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <cstdio>

namespace util {

FuzzyIndices::FuzzyIndices(
		const std::shared_ptr<std::vector<util::FuzzyCluster> > &fuzzyClusters) :
		fuzzyClusters(fuzzyClusters) {

}

FuzzyIndices::CampelloIndices::CampelloIndices(FuzzyIndices& outer) :
		outer(outer) {
	std::shared_ptr<std::vector<int>> elements = outer.elements();
	a = compute_a(*elements.get());
	b = compute_b(*elements.get());
	c = compute_c(*elements.get());
	d = compute_d(*elements.get());
}

double FuzzyIndices::CampelloIndices::computeV(const int j1, const int j2) const {
	double result = 0.0;

	for (std::set<int>::const_iterator it = outer.partitionIndices.begin();
			it != outer.partitionIndices.end(); ++it) {
		const int classNum = *it;
		result = std::max(result,
				std::min(
						outer.getPertinenceDegree(std::make_pair(j1, classNum)),
						outer.getPertinenceDegree(
								std::make_pair(j2, classNum))));
	}
	//printf(" V = %lf ", result);
	return result;
}

double FuzzyIndices::CampelloIndices::computeX(const int j1, const int j2) const {
	double result = 0.0;

	for (std::set<int>::const_iterator it = outer.partitionIndices.begin();
			it != outer.partitionIndices.end(); ++it) {
		const int i1 = *it;
		for (std::set<int>::const_iterator it2 = outer.partitionIndices.begin();
				it2 != outer.partitionIndices.end(); ++it2) {
			const int i2 = *it2;
			if (i1 != i2) {
				result = std::max(result,
						std::min(
								outer.getPertinenceDegree(
										std::make_pair(j1, i1)),
								outer.getPertinenceDegree(
										std::make_pair(j2, i2))));
			}

		}

	}
	//printf(" X = %lf ", result);
	return result;
}

double FuzzyIndices::CampelloIndices::computeY(const int j1, const int j2) const {
	double result = 0.0;
	for (unsigned int i = 0; i < outer.fuzzyClusters->size(); i++) {
		result = std::max(result,
				std::min(outer.fuzzyClusters->at(i).getMembershipDegree(j1),
						outer.fuzzyClusters->at(i).getMembershipDegree(j2)));
	}
	//printf(" Y = %lf ", result);
	return result;
}

double FuzzyIndices::CampelloIndices::computeZ(const int j1, const int j2) const {
	double result = 0.0;

	for (unsigned int i = 0; i < outer.fuzzyClusters->size(); i++) {
		for (unsigned int j = 0; j < outer.fuzzyClusters->size(); j++) {
			if (i != j) {
				result = std::max(result,
						std::min(
								outer.fuzzyClusters->at(i).getMembershipDegree(
										j1),
								outer.fuzzyClusters->at(j).getMembershipDegree(
										j2)));

			}
		}

	}
	//printf(" Z = %lf ", result);

	return result;
}

std::shared_ptr<std::vector<int>> FuzzyIndices::elements() const {
	std::shared_ptr<std::set<int>> elements = fuzzyClusters->at(0).getElements();
	std::shared_ptr<std::vector<int>> result(new std::vector<int>(elements->begin(), elements->end()));
	return result;
}

double FuzzyIndices::CampelloIndices::compute_a(const std::vector<int> &elements) const {
	double result = 0.0;
	for (size_t j1 = 0; j1 < elements.size(); j1++) {
		for(size_t j2 = j1+1; j2 < elements.size(); j2++) {
			result += std::min(computeV(elements[j1], elements[j2]), computeY(elements[j1], elements[j2]));
		}
	}

	return result;
}

double FuzzyIndices::CampelloIndices::compute_b(const std::vector<int> &elements) const {
	double result = 0.0;
	for (size_t j1 = 0; j1 < elements.size(); j1++) {
		for(size_t j2 = j1+1; j2 < elements.size(); j2++) {
			result += std::min(computeV(elements[j1], elements[j2]), computeZ(elements[j1], elements[j2]));
		}
	}

	return result;
}

double FuzzyIndices::CampelloIndices::compute_c(const std::vector<int> &elements) const {
	double result = 0.0;
	for (size_t j1 = 0; j1 < elements.size(); j1++) {
		for(size_t j2 = j1+1; j2 < elements.size(); j2++) {
			result += std::min(computeX(elements[j1], elements[j2]), computeY(elements[j1], elements[j2]));
		}
	}

	return result;
}

double FuzzyIndices::CampelloIndices::randIndex() const {
	const double num = a + d;
	const double den = a + b + c + d;
	return num/den;
}

double FuzzyIndices::CampelloIndices::adjustedRandIndex() const {
	printf(" (a = %lf, b = %lf, c = %lf, d = %lf) ", a, b, c, d);
	const double num = a - ((a+c)*(a+b))/d;
	const double den = ((a+c)+(a+b))/2.0 - ((a+c)*(a+b))/d;
	return num/den;
}

double FuzzyIndices::CampelloIndices::jaccardCoef() const {
	return a / (a + b + c);
}

double FuzzyIndices::CampelloIndices::fowlkesMallowsIndex() const {
	const double den = sqrt((a + b)*(a + c));
	return a / den;
}

double FuzzyIndices::CampelloIndices::minkowskiMeasure() const {
	const double num = b + c;
	const double den = b + a;
	return sqrt(num/den);
}

double FuzzyIndices::CampelloIndices::compute_d(const std::vector<int> &elements) const {
	double result = 0.0;
	for (size_t j1 = 0; j1 < elements.size(); j1++) {
		for(size_t j2 = j1+1; j2 < elements.size(); j2++) {
			result += std::min(computeX(elements[j1], elements[j2]), computeZ(elements[j1], elements[j2]));
		}
	}

	return result;
}

FuzzyIndices::AndersonIndices::AndersonIndices(FuzzyIndices& outer) : outer(outer) {
	std::shared_ptr<std::vector<int>> elements = outer.elements();
	const size_t n = elements->size();
	arma::mat U(outer.numberOfClusters(), n);
	arma::mat V(outer.numberOfPartitionsReference(), n);
	const int r = outer.numberOfClusters();
	for (int i = 0; i < r; i++) {
		const util::FuzzyCluster &fc = outer.fuzzyClusters->at(i);
		for (size_t k = 0; k < elements->size(); k++) {
			U(i, k) = fc.getMembershipDegree(elements->at(k));
		}
	}
	int c = 0;
	{

		for (std::set<int>::const_iterator rIterator =
				outer.partitionIndices.begin();
				rIterator != outer.partitionIndices.end(); rIterator++) {
			for (size_t k = 0; k < n; k++) {
				V(c, k) = outer.getPertinenceDegree(std::make_pair(elements->at(k), (*rIterator)));
			}
			c++;
		}
	}
	arma::mat N = U * V.st();

	double n_i_dot_sum = 0.0;
	for (int i = 0; i < r; i++) {
		const double val = accu(N.row(i));
		n_i_dot_sum += val;
	}

	double phi = n/n_i_dot_sum;
	std::cout << "phi is " << phi << std::endl;
	N *= phi;
	indexA = 0.0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			double val = N(i,j);
			indexA += val * (val - 1.0);
		}
	}
	indexA *= 0.5;
	double Nsquare = accu(square(N));

	double n_i_dot_squared_sum = 0.0;
	for (int i = 0; i < r; i++) {
		const double val = accu(N.row(i));
		n_i_dot_squared_sum += val*val;
	}

	double n_dot_j_squared_sum = 0.0;
	for (int j = 0; j < c; j++) {
		const double val = accu(N.col(j));
		n_dot_j_squared_sum += val*val;
	}
	const double n_squared = n*n;
	indexD = 0.5*(n_squared + Nsquare - (n_i_dot_squared_sum + n_dot_j_squared_sum));
	indexB = (n_dot_j_squared_sum - Nsquare)*0.5;
	indexC = (n_i_dot_squared_sum - Nsquare)*0.5;

	chiSq_4m = 0.0;
	for (int j = 0; j < c; ++j) {
		for (int i = 0; i < r; ++i) {
			const double val = N(i, j);
			chiSq_4m += ((val*val)/(accu(N.row(i))*accu(N.col(j))));
		}
	}
	chiSq_4m -= 1.0;
	chiSq_4m *= n;

	gdmnKrkl_4n = 0.0;
	for (int j = 0; j < c; ++j) {
			for (int i = 0; i < r; ++i) {
				const double val = N(i, j);
				gdmnKrkl_4n += (val*val)/accu(N.row(i));
			}
	}
	gdmnKrkl_4n *= n;
	double minusSum = 0.0;
	for (int j = 0; j < c; ++j) {
		const double val = accu(N.col(j));
		minusSum += val*val;
	}
	gdmnKrkl_4n -= minusSum;
	double den = n_squared - n_dot_j_squared_sum;
	gdmnKrkl_4n /= den;

	sFRHR = 0.0;
	for (size_t i = 0; i < n-1; ++i) {
		for (size_t j = i+1; j < n; ++j) {
			sFRHR += fabs(accu(V.col(i) - V.col(j)) - accu(U.col(i) - U.col(j)));
		}
	}
	sFRHR /= n*(n-1.0)*0.5;
	sFRHR = 1.0 - sFRHR;
}
}

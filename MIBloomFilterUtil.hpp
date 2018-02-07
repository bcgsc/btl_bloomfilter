#ifndef MIBLOOMFILTERUTIL
#define MIBLOOMFILTERUTIL 1

#include "MIBloomFilter.hpp"
//#include "ntHashIterator.hpp"
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include <boost/math/distributions/binomial.hpp>
#include <vector>
#include <limits>

using namespace std;
using boost::math::binomial;

namespace MIBloomFilterUtil {

static const unsigned s_emptyID = 0;

//from https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
static unsigned nChoosek(unsigned n, unsigned k) {
	if (k > n)
		return 0;
	if (k * 2 > n)
		k = n - k;
	if (k == 0)
		return 1;

	int result = n;
	for (unsigned i = 2; i <= k; ++i) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}

//helper methods
//calculates the per frame probability of a random match for single value
static double calcProbSingleFrame(double occupancy, unsigned hashNum,
		double freq, unsigned allowedMisses = 0) {
	double probTotal = 0.0;
	for (unsigned i = hashNum - allowedMisses; i <= hashNum; i++) {
		double prob = nChoosek(hashNum, i);
		prob *= pow(occupancy, i);
		prob *= pow(1.0 - occupancy, hashNum - i);
		prob *= (1.0 - pow(1.0 - freq, i));
		probTotal += prob;
	}
	return probTotal;
}

//calculates the per frame probability of a random multi match given a significant result
static double calcProbMultiMatchSingleFrame(double occupancy, unsigned hashNum,
		double freq) {
	double prob = 1.0
			- pow(1.0 - freq, hashNum * (1 + occupancy / log(1 - occupancy)));
	return prob;
}

/*
 * Max value is the largest value seen in your set of possible values
 */
template<typename T>
static vector<double> calcPerMultiMatchFrameProb(MIBloomFilter<T> &miBF,
		T maxValue) {
	double occupancy = double(miBF.getPop()) / double(miBF.size());
	unsigned hashNum = miBF.getHashNum();
	vector<size_t> countTable = vector<size_t>(maxValue + 1, 0);
	miBF.getIDCounts(countTable);
	size_t sum = 0;
	for (vector<size_t>::const_iterator itr = countTable.begin();
			itr != countTable.end(); ++itr) {
		sum += *itr;
	}
	vector<double> perFrameProb = vector<double>(maxValue + 1, 0.0);
	for (size_t i = 0; i < countTable.size(); ++i) {
		perFrameProb[i] = MIBloomFilterUtil::calcProbMultiMatchSingleFrame(
				occupancy, hashNum, double(countTable[i]) / double(sum));
//		cerr << double(countTable[i]) / double(sum) << " " << perFrameProb[i] << endl;
	}

	return perFrameProb;
}

/*
 * Max value is the largest value seen in your set of possible values
 */
template<typename T>
static vector<double> calcPerFrameProb(MIBloomFilter<T> &miBF, T maxValue) {
	double occupancy = double(miBF.getPop()) / double(miBF.size());
	unsigned hashNum = miBF.getHashNum();
	vector<size_t> countTable = vector<size_t>(maxValue + 1, 0);
	miBF.getIDCounts(countTable);
	size_t sum = 0;
	for (vector<size_t>::const_iterator itr = countTable.begin();
			itr != countTable.end(); ++itr) {
		sum += *itr;
	}
	vector<double> perFrameProb = vector<double>(maxValue + 1, 0.0);
	for (size_t i = 0; i < countTable.size(); ++i) {
		perFrameProb[i] = MIBloomFilterUtil::calcProbSingleFrame(occupancy,
				hashNum, double(countTable[i]) / double(sum));
//		cerr << double(countTable[i]) / double(sum) << " " << perFrameProb[i] << endl;
	}

	return perFrameProb;
}

/*
 * Utility function for querying MiBF for a sequence in a hash itr object
 * Computes probability of a match for each possible value queried.
 * In this scheme saturated regions are ignored.
 *
 * Parameters:
 * perFrameProb - per frame probability of each possible value
 * alpha - significance threshold
 * itr - hash iterator of type H (ntHash)
 * maxPos - max number of positions to move on the sequence
 *
 * Returns significant values (smaller than alpha threshold)
 */
//TODO return saturated frame counts?
//TODO return pVals?
template<typename T, typename H>
static vector<T> query(MIBloomFilter<T> &miBF, H &itr,
		const vector<double> &perFrameProb, const vector<double> &perMultiMatchFrameProb,
		size_t maxPos = numeric_limits<size_t>::max(), double alpha = 0.0001,
		double multimapAlpha = 0.001) {
	unsigned evaluatedSeeds = 0;

	google::dense_hash_map<T, unsigned> counts;
	counts.set_empty_key(s_emptyID);
	while (itr != itr.end() && itr.pos() < maxPos) {
		bool saturated = true;
		vector<T> results = miBF.at(*itr, saturated);
		//to determine if already added for this frame
		google::dense_hash_set<T> tempIDs;
		tempIDs.set_empty_key(s_emptyID);

		if (!saturated) {
			for (typename vector<T>::const_iterator j = results.begin();
					j != results.end(); j++) {
				if (*j != s_emptyID) {
					if (tempIDs.find(*j) == tempIDs.end()) {
						typename google::dense_hash_map<T, unsigned>::iterator tempItr =
								counts.find(*j);
						assert(*j > 0);
						if (tempItr == counts.end()) {
							counts[*j] = 1;
						} else {
							++(tempItr->second);
						}
						tempIDs.insert(*j);
					}
				}
			}
			++evaluatedSeeds;
		}
		++itr;
	}

	//potential signifResults
	vector<T> potSignifResults;
	vector<T> signifResults;

	double adjustedPValThreshold = 1.0
			- pow(1.0 - alpha, 1.0 / double(perFrameProb.size() - 1));
	T bestSignifVal = counts.begin()->first;
	for (typename google::dense_hash_map<T, unsigned>::const_iterator itr =
			counts.begin(); itr != counts.end(); itr++) {
		//TODO use complement cdf? so I don't have to subtract?
		binomial bin(evaluatedSeeds, 1.0 - perFrameProb.at(itr->first));
		double cumProb = cdf(bin, evaluatedSeeds - itr->second);
		if (adjustedPValThreshold > cumProb) {
			if (counts[bestSignifVal] < counts[itr->first]) {
				bestSignifVal = itr->first;
			}
			potSignifResults.push_back(itr->first);
		}
//		cerr << unsigned(itr->first) << " " << adjustedPValThreshold << " "
//				<< cumProb << " " << itr->second << " " << miBF.getPop()
//				<< endl;
	}

	adjustedPValThreshold = 1.0
			- pow(1.0 - multimapAlpha, 1.0 / double(perFrameProb.size() - 1));
	//TODO: generalized because this assumes a = 0, fix me?
	for (typename vector<T>::const_iterator itr = potSignifResults.begin();
			itr != potSignifResults.end(); ++itr) {
		//compute single frame prob
		binomial bin(counts[bestSignifVal], 1.0 - perMultiMatchFrameProb.at(*itr));
		double cumProb = cdf(bin, counts[bestSignifVal] - counts[*itr]);
		if (adjustedPValThreshold > cumProb) {
			signifResults.push_back(*itr);
		}
//		cerr << unsigned(*itr) << " " << counts[bestSignifVal] << " "
//				<< counts[*itr] << " " << adjustedPValThreshold << " "
//				<< multimapAlpha << " " << cumProb << endl;
	}

	//Best hit considered the class with the most hits
	return signifResults;
}

}

#endif

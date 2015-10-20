/*
 * BloomFilterInfo.hpp
 * Intended to be used to collect/store and compute information about bloom filters
 * Can output information in IEE format text file
 *
 *  Created on: Aug 20, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTERINFO_H_
#define BLOOMFILTERINFO_H_
#include <string>
#include <vector>
//#include </property_tree/ptree.hpp>
//#include </property_tree/ini_parser.hpp>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "SimpleIni.h"

using namespace std;

// For Calculating aspects of bloom filter
// todo: Tweak calculations as they are approximations and may not be 100% optimal
// see http://en.wikipedia.org/wiki/Bloom_filter

/*
 * Calculation assumes optimal ratio of bytes per entry given a fpr
 */
//Note: Rounded down because in practice you want to calculate as few hash values as possible
//static uint8_t calcOptiHashNum(float fpr) {
//	return uint8_t(-log(fpr) / log(2));
//}
class BloomFilterInfo {
public:

	BloomFilterInfo(string const &filterID, unsigned kmerSize, unsigned hashNum,
			double desiredFPR, size_t expectedNumEntries,
			const vector<string> &seqSrcs) :
			m_filterID(filterID), m_kmerSize(kmerSize), m_desiredFPR(desiredFPR), m_seqSrcs(
					seqSrcs), m_hashNum(hashNum), m_expectedNumEntries(
					expectedNumEntries)
	{
		m_runInfo.size = calcOptimalSize(expectedNumEntries, desiredFPR, hashNum);
		m_runInfo.redundantSequences = 0;
	}

	/*
	 * loads bloom filter information from a file
	 */
	//Todo: convert to having variables stored in property tree for more modularity
	BloomFilterInfo(string const &fileName)
	{
		CSimpleIniA reader;
		reader.SetUnicode();
		reader.LoadFile(fileName.c_str());

		m_filterID = reader.GetValue("user_input_options", "filter_id", "");
		m_kmerSize = reader.GetLongValue("user_input_options", "kmer_size", 0);
		m_desiredFPR = reader.GetDoubleValue("user_input_options", "desired_false_positve_rate", 0);
		string tempSeqSrcs = reader.GetValue("user_input_options", "sequence_sources", "");
		m_seqSrcs = convertSeqSrcString(tempSeqSrcs);
		m_hashNum = reader.GetLongValue("user_input_options", "number_of_hash_functions", 0);

		//runtime params
		m_runInfo.size = reader.GetLongValue("runtime_options", "size", 0);
		m_runInfo.numEntries = reader.GetLongValue("runtime_options", "num_entries", 0);

		m_runInfo.redundantSequences = reader.GetLongValue("runtime_options", "redundant_sequences", 0);
		m_runInfo.redundantFPR = reader.GetDoubleValue("runtime_options", "redundant_fpr", 0);
		m_expectedNumEntries = reader.GetLongValue("user_input_options", "expected_num_entries", 0);
		m_runInfo.FPR = reader.GetDoubleValue("runtime_options", "approximate_false_positive_rate", 0);
	}

	/**
	 * Sets number of redundant sequences found in file. Also calculate the approximate
	 * error that may be happening to this value.
	 */
	void setRedundancy(size_t redunSeq)
	{
		assert(m_runInfo.numEntries > 0);
		m_runInfo.redundantSequences = redunSeq;
		m_runInfo.redundantFPR = calcRedunancyFPR(m_runInfo.size, m_runInfo.numEntries,
				m_hashNum);

		m_runInfo.FPR = calcApproxFPR(m_runInfo.size,
				m_runInfo.numEntries, m_hashNum);
	}

	/**
	 * Sets number of elements inserted into filter
	 */
	void setTotalNum(size_t totalNum)
	{
		m_runInfo.numEntries = totalNum;
	}

	/*
	 * Prints out INI format file
	 */
	//Todo: research better method of outputting INI format file
	void printInfoFile(const string &fileName) const
	{

		//cannot output unless runtime values set
		assert(m_runInfo.size !=0);
		assert(m_runInfo.numEntries !=0);
		assert(m_expectedNumEntries !=0);
		assert(m_runInfo.FPR !=0);
		assert(m_hashNum > 0);

		ofstream output(fileName.c_str(), ios::out);
		//user specified
		output << "[user_input_options]\nfilter_id=" << m_filterID << "\nkmer_size="
				<< m_kmerSize << "\ndesired_false_positve_rate=" << m_desiredFPR
				<< "\nnumber_of_hash_functions=" << m_hashNum
				<< "\nexpected_num_entries=" << m_expectedNumEntries
				<< "\nsequence_sources=";

		//print out sources as a list
		for (vector<string>::const_iterator it = m_seqSrcs.begin();
				it != m_seqSrcs.end(); ++it)
		{
			output << *it;
			output << " ";
		}

		//runtime determined options
		output << "\n\n[runtime_options]\nsize=" << m_runInfo.size << "\nnum_entries="
				<< m_runInfo.numEntries << "\napproximate_false_positive_rate="
				<< m_runInfo.FPR << "\nredundant_sequences="
				<< m_runInfo.redundantSequences << "\nredundant_fpr="
				<< m_runInfo.redundantFPR << "\n";
		//print out hash functions as a list

		output.close();
	}

	//getters

	unsigned getKmerSize() const
	{
		return m_kmerSize;
	}

	unsigned getHashNum() const
	{
		return m_hashNum;
	}

	size_t getCalcuatedFilterSize() const
	{
		return m_runInfo.size;
	}

	const string &getFilterID() const
	{
		return m_filterID;
	}

	double getRedundancyFPR() const
	{
		return m_runInfo.redundantFPR;
	}

	double getFPR() const
	{
		return m_runInfo.FPR;
	}

	~BloomFilterInfo()
	{
	}

private:
	//user specified input
	string m_filterID;
	unsigned m_kmerSize;
	double m_desiredFPR;
	vector<string> m_seqSrcs;
	unsigned m_hashNum;
	size_t m_expectedNumEntries;

	//determined at run time
	struct runtime {
		size_t size;
		size_t numEntries;
		double FPR;
		size_t redundantSequences;
		double redundantFPR;
	};

	runtime m_runInfo;

	const vector<string> convertSeqSrcString(
			string const &seqSrcStr) const
	{
		vector<string> inputs;
		string currentFileName = "";
		string temp;
		stringstream converter(seqSrcStr);
		while (converter >> temp) {
			inputs.push_back(temp);
		}
		return inputs;
	}

	// functions for calculations regarding bloomfilter
	// todo: Tweak calculations as they are approximations and may not be 100% optimal
	// see http://en.wikipedia.org/wiki/Bloom_filter

	//Private functions
	/*
	 * Calculate FPR based on hash functions, size and number of entries
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	double calcApproxFPR(size_t size, size_t numEntr,
			unsigned hashFunctNum) const
	{
		return pow(
				1.0 - pow(1.0 - 1.0 / double(size), double(numEntr) * hashFunctNum),
				double(hashFunctNum));
	}

	/*
	 * Calculates redundancy FPR
	 */
	double calcRedunancyFPR(size_t size, size_t numEntr,
			unsigned hashFunctNum) const
	{
		double total = log(calcApproxFPR(size, 1, hashFunctNum));
		for (size_t i = 2; i < numEntr; ++i) {
			total = log(exp(total) + calcApproxFPR(size, i, hashFunctNum));
		}
		return exp(total) / numEntr;
	}

	/*
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * assuming optimal # of hash functions used
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	//NOTE: Not currently used.
	size_t calcOptimalSize(size_t entries, float fpr) const
	{
		size_t non64ApproxVal = size_t(entries * -log(fpr) / pow(log(2), 2));
		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}

	/*
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * given the number of hash functions
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	size_t calcOptimalSize(size_t entries, float fpr,
			unsigned hashNum) const
	{
		size_t non64ApproxVal = size_t(
				-double(entries) * double(hashNum)
						/ log(1.0 - pow(fpr, float(1 / (float(hashNum))))));

		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}
};

#endif /* BLOOMFILTERINFO_H_ */

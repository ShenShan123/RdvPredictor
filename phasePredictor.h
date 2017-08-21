#ifndef __MEM_TRACE_SIMPLE_H__
#define __MEM_TRACE_SIMPLE_H__

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <map>
#include <random>
#include <unordered_map>
#include <assert.h>
#include <stdlib.h> 
#include <iomanip>
#include <deque>
#include <list>
#include <float.h>
#include "pin.H"

static uint64_t Truncation = 0;
static uint64_t InterCount = 0;
static uint64_t NumMemAccs = 0;
static uint64_t NumIntervals = 0;
static uint64_t IntervalSize = 0;
/* block size is 64 byte, so the block offset bits is lower 6 bits */
static uint64_t BlkBits = 6;

/* parse the command line arguments */
KNOB<string> KnobOutputFile(KNOB_MODE_WRITEONCE, "pintool", "o", "RDV.txt", "specify output file name");
KNOB<UINT64> KnobTruncDist(KNOB_MODE_WRITEONCE, "pintool", "m", "16384", "the truncation distance of SD");
KNOB<UINT64> KnobIntervalSize(KNOB_MODE_WRITEONCE, "pintool", "i", "10000000", "the interval size");
//KNOB<UINT64> KnobSampleRate(KNOB_MODE_WRITEONCE, "pintool", "s", "10000", "the sample rate");
KNOB<UINT32> KnobRdvThreshold(KNOB_MODE_WRITEONCE, "pintool", "t", "5", "the maximum normalized manhattan distance of two RD vector");
KNOB<UINT32> KnobPhaseTableSize(KNOB_MODE_WRITEONCE, "pintool", "p", "10000", "phase table size");
KNOB<UINT32> KnobHistoryLength(KNOB_MODE_WRITEONCE, "pintool", "l", "10", "history register length");
KNOB<UINT32> KnobHistoryTableSize(KNOB_MODE_WRITEONCE, "pintool", "s", "65536", "history signature table size");


#define LOG2
//#define SAMPLE

/* calculate log2(s) + 1 */
template<class T>
inline T log2p1(T s)
{
	T result = 0;
	while (s) {
		s >>= 1;
		++result;
	}

	return result;
}

/* fast to calculate log2(x)+1 */
#ifdef LOG2
#define DOLOG(x) log2p1(x)
#else
#define DOLOG(x) x
#endif

/* for recording distribution into a Histogram, 
   Accur is the accuracy of transforming calculation */
template <class B = int64_t>
class Histogram
{
protected:
	B * bins;
	B samples;
	int _size;

public:
	Histogram() : bins(nullptr), samples(0), _size(0) {};

	Histogram(int s);

	Histogram(const Histogram<B> & rhs);

	~Histogram();

	void setSize(int s);

	const int size() const;

	void clear();

	void normalize();

	B & operator[] (const int idx) const;

	Histogram<B> & operator=(const Histogram<B> & rhs);

	Histogram<B> & operator+=(const Histogram<B> & rhs);

	void sample(int x);

	const B getSamples() const; 

	void print(std::ofstream & file);
};

/* do reuse distance statistics */
class ReuseDist
{
	std::map<uint64_t, long> addrMap;
	long index;

public:
	ReuseDist() : index(0) {};

	~ReuseDist() {};

	void calReuseDist(uint64_t addr, Histogram<> & rdv);
};

class PhaseTable
{
private:
	template<class B>
	class Entry : public Histogram<B>
	{
	private:
		friend class PhaseTable;
		uint32_t id;
		uint32_t occur;
	
	public:
		Entry(const Histogram<B> & rdv, const uint32_t _id);
	
		~Entry() {};

		double manhattanDist(const Histogram<B> & rhs);
	};

	/* list for LRU replacement policy */
	std::list<Entry<int64_t> *> pt;
	double threshold;
	uint32_t index;
	uint32_t ptSize;

public:
	PhaseTable() : threshold(0.0), index(0), ptSize(0) {};

	~PhaseTable();

	void init(double t, uint32_t s);

	/* do LRU replacement */
	void lruRepl(Entry <int64_t> * ent, std::list<Entry<int64_t> *>::iterator & p);

	uint32_t find(const Histogram<> & rdv);

	const int size() const { return pt.size(); }
};

class Predictor
{
private:
	/* history registers */
	std::deque<uint32_t> histRegs;
	uint32_t histLength;
	/* history signature table */
	uint32_t * hst;
	uint32_t hstSize;
	uint64_t index;
	uint32_t predOutCome;

public:
	uint32_t misPredicts;

	Predictor() : histLength(0), hstSize(0), index(0), predOutCome(0), misPredicts(0) {};

	~Predictor() { delete [] hst; }

	void init(uint32_t l, uint32_t s);

	uint32_t predict(uint32_t id);
};

VOID PIN_FAST_ANALYSIS_CALL
RecordMemRefs(ADDRINT ea);

/*
 * Insert code to write data to a thread-specific buffer for instructions
 * that access memory.
 */
VOID Trace(TRACE trace, VOID *v);

/* output the results, and free the poiters */
VOID Fini(INT32 code, VOID *v);

/* ===================================================================== */
/* Print Help Message                                                    */
/* ===================================================================== */
INT32 Usage();


/* global variates */
std::ofstream fout;
Histogram<> currRDD;
ReuseDist reuseDist;
PhaseTable phaseTable;
Predictor predictor;

#endif
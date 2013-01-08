// ==========================================================================
//                                    dmx
// ==========================================================================
// Copyright (c) 2012, Jay DePasse, University of Pittsburgh
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jay DePasse <jvd10@pitt.edu>
// ==========================================================================

#ifndef SANDBOX_JVD_APPS_DMX_DMXCORE_H_
#define SANDBOX_JVD_APPS_DMX_DMXCORE_H_


#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/index.h>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>
#include <seqan/refinement.h>
#include <seqan/score.h>

#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_do.h>
#include <tbb/parallel_sort.h>
#include <tbb/concurrent_priority_queue.h>

#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "dmxBarcode.h"
#include "dmxIO.h"
#include "dmxRead.h"

#include <ltiClustering.h>
#include <ltiL2Distance.h>
#include <ltiKMeansClustering.h>
#include <ltiDBScan.h>
#include <ltiAdaptiveKMeans.h>

using namespace seqan;
using namespace tbb;

//enum barcodeAssignmentType { BOTH, FWD, REV, NO_MATCH, MISMATCH } ;


struct fastqPair {
  std::string id1, id2, sq1, sq2, ql1, ql2;
  unsigned num;
  void reset() {
    id1.erase();
    id2.erase();
    sq1.erase();
    sq2.erase();
    ql1.erase();
    ql2.erase();
  }
};

struct dmxMatch {
  unsigned min;
  int index;
};

typedef concurrent_vector< dmxRead * > dmxReadVector; 
typedef concurrent_priority_queue< dmxRead *, dmxReadCompare > dmxReadPriQ; 
typedef std::vector< dmxRead * > dmxReadSerialVector; 


typedef StringSet< String< char > > barcodeStringSetType;
typedef Index< StringSet< String< char > > > barcodeStringSetIndexType;
typedef Finder< Index< StringSet< String< char > > > > barcodeStringSetIndexFinderType;


class dmx {

  public:
    static const int seqTagLength = 20;

    dmx( char* barcodeFile );
    void initFastq( unsigned _maxDistance, unsigned _chunkSize, unsigned _trimSize );
    void runFastq( char* pair1FileName, char* pair2FileName );
    
    unsigned readCount; 
    unsigned maxDistance;
    void digest();
    void digest(barcodeStringSetIndexFinderType * _barcodeFinder, std::vector< fastqPair > * fastqFeedChunk);

    void parallelDigest2();

    void printBarcodeResults( dmxReadVector resultVector );
    void printBarcodeResults( dmxReadPriQ &resultVector );

    void printConcordantBarcodeResults();
    void printFwdOnlyBarcodeResults();
    void printRevOnlyBarcodeResults();
    void printDiscordantBarcodeResults();
    void printUnidentifiableBarcodeResults();

    void printFasta( dmxReadSerialVector, std::ofstream & fh );
    void printFastq( dmxReadSerialVector, std::ofstream & fh );

    void printFastq( dmxReadSerialVector, std::ofstream & fh, int barcode_index );

    void printFasta( dmxReadSerialVector, std::ofstream & fh, int barcode_index );

    void printFasta( dmxReadPriQ&, std::ofstream & fh );
    void printFasta( dmxReadPriQ&, std::ofstream & fh, int barcode_index );

    void printGoodFasta( std::string filename );
    void printGoodFastq( std::string filename );

    void printAllFasta( std::string filename );

    void printPerBarcodeFasta( std::string outfilePrefix );

    unsigned chunkSize, trimSize; 
    unsigned int distance(const std::string s1, const std::string s2);
    void read2FilePairedFastq( char * pair1FileName, char * pair2FileName );

    int readBarcodeFile(char* barcodeFile);
    dmxMatch getMatch(std::string seq);
    dmxMatch getExactMatch(std::string seq);
    dmxMatch getMatchMyersInfix(std::string seq);
    dmxMatch getMatchIndexFinder(std::string seq,  barcodeStringSetIndexFinderType * _barcodeFinder);

    bool pairedEnd;
    atomic< bool > finishedReading;
    atomic< bool > spoon;

    MultiSeqFile pair1File;
    MultiSeqFile pair2File;

    std::map< std::string, barcode > barcodes;
    std::vector< std::string > barcodeNames;

    dmxReadPriQ fwdBarcode;
    dmxReadPriQ revBarcode;
    dmxReadPriQ conBarcode;
    dmxReadPriQ disBarcode;
    dmxReadPriQ nonBarcode;

    dmxReadSerialVector fwdBarcodeSerVec;
    dmxReadSerialVector revBarcodeSerVec;
    dmxReadSerialVector conBarcodeSerVec;
    dmxReadSerialVector disBarcodeSerVec;
    dmxReadSerialVector nonBarcodeSerVec;

    void convertPriorityQueueToVector( dmxReadPriQ & q, dmxReadSerialVector & v );
    void convertPriorityQueuesToVectors();

    concurrent_queue< std::vector< fastqPair > * > fastqChunks;
    concurrent_vector< std::vector< fastqPair > * > fastqFeed;

    barcodeStringSetType barcodeStringSet;
    barcodeStringSetIndexType barcodeStringSetIndex;
    barcodeStringSetIndexFinderType barcodeFinder;

    dmxIO * dmxio;

    void cluster_test();
    
    void groupReduce();
    void groupReduce( dmxReadSerialVector * v, dmxReadPriQ * q );

    struct groupReduceFunctor {
      dmx * d;
      dmxReadSerialVector * drsv;
      dmxReadPriQ * drpq;
      tbb::atomic< bool > * group_reduce_spoon;

      groupReduceFunctor( dmx * _d, dmxReadSerialVector * _drsv, dmxReadPriQ * _drpq, tbb::atomic< bool > * _spoon ) {
        d = _d;
        drsv = _drsv;
        drpq = _drpq;
        group_reduce_spoon = _spoon;
      }

      void operator()( dmxReadSerialVector * r, parallel_do_feeder< dmxReadSerialVector * >& feeder ) const;
    };

    void getClusters( std::map< int, dmxReadSerialVector > & clusterMap, dmxReadSerialVector * rv );
    dmxRead * condenseGroup( std::vector< dmxRead * > & rv );
    void test_consensus();
    std::string computeConsensus( std::string & matrix, size_t nrow );
  };

  struct digestFunctor2 {

    typedef std::vector< fastqPair > * argument_type;

    dmx * d;
    void operator()( argument_type at, parallel_do_feeder< argument_type >& feeder ) const {

    barcodeStringSetIndexFinderType _barcodeFinder( d->barcodeFinder );
   
    if ( d->spoon.compare_and_swap( false, true ) ) {
      // I have the spoon, so I must be the feeder; following loop only entered by the feeder:
      while ( !(d->finishedReading) ) {
        int numChunks = d->fastqChunks.unsafe_size();
        for ( int i = 0; i < numChunks; ++i ) {
          std::vector< fastqPair > * chunk;
          if ( d->fastqChunks.try_pop( chunk ) ) {
            feeder.add( chunk );
          }
        }
        // check safe for concurrent use; compact the internal representation of the concurrent vector
        d->fastqFeed.shrink_to_fit();
      }
    }
    // executed by all tasks, including, eventually, the feeder:
    d->digest( & _barcodeFinder, at ); 
    at->clear();
    delete at;
  }
};


#endif  // #ifndef SANDBOX_JVD_APPS_DMX_DMXCORE_H_









































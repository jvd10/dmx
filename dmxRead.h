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

#ifndef SANDBOX_JVD_APPS_DMX_DMXREAD_H_
#define SANDBOX_JVD_APPS_DMX_DMXREAD_H_

#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <fstream>

enum barcodeAssignmentType { BOTH, FWD, REV, NO_MATCH, MISMATCH } ;

class dmxRead {

public:

  /*
   * Constructors, Destructor, clone methods
   */
  dmxRead();
  dmxRead * newClone();
  dmxRead( barcodeAssignmentType _descriptionCode, std::string _tag, unsigned _fIdx, unsigned _rIdx );
 
  /*
   * Initialization/Set methods
   */
  void fwd( int _fBCidx, std::string _fSeq );
  void rev( int _rBCidx, std::string _rSeq );

  void fwd( int _fBCidx, std::string _fSeq, std::string _fQual );
  void rev( int _rBCidx, std::string _rSeq, std::string _rQual );
   
  /*
   * Access/Get methods
   */
  std::string getDescription();
  std::string getShortDescription();
  void getDinucleotideFreqs( std::vector< double > & kmer );
  void getDinucleotideFreqs( std::string s, std::vector< double > & kmer );

  int getFwdBCidx() { return fBCidx; }
  int getRevBCidx() { return rBCidx; }

  /*
   * Operators
   */
  bool operator== ( dmxRead & other );

  /*
   * Printing and debugging
   */
  void print();
  void printFFasta( unsigned i, std::ofstream & fh );
  void printRFasta( unsigned i, std::ofstream & fh );
  void printFFastq( unsigned i, std::ofstream & fh );
  void printRFastq( unsigned i, std::ofstream & fh );
  void printFasta( unsigned i, std::ofstream & fh );
  void printFastq( unsigned i, std::ofstream & fh );

  //TODO Eventually all data members below should be private TODO//

  unsigned fIdx, rIdx;
  std::string tag, fSeq, rSeq;
  std::string fQual, rQual;
  barcodeAssignmentType descriptionCode;
  /* 
   * groupSize is the number of reads that share barcode, random counter and random primer 
   * cluster size is the size of each of the clusters (based on remaining sequence) identified
   * in the alignment of the group and then clustered based on kmer distribution
   */
  unsigned groupSize, clusterSize;

  //std::string fId, rId;

private:

  short int fBCidx;
  short int rBCidx;

};

struct dmxReadCompare {
  
  enum cmp { LT, EQ, GT };

  bool operator() ( dmxRead * x, dmxRead * y ) const {
    if (x->descriptionCode != y->descriptionCode) {
      std::cout << "Trying to compare " << x->descriptionCode << " to " << y->descriptionCode << " in read sort." << std::endl;
    }
    cmp fCmp = EQ;
    cmp rCmp = EQ;
    if ( x->descriptionCode != REV ) {
      if ( x->getFwdBCidx() < y->getFwdBCidx() ) { fCmp = LT; }
      else if ( x->getFwdBCidx() > y->getFwdBCidx() ) { fCmp = GT; }
    }
    if ( x->descriptionCode != FWD ) {
      if ( x->getRevBCidx() < y->getRevBCidx() ) { rCmp = LT; }
      else if ( x->getRevBCidx() > y->getRevBCidx() ) { rCmp = GT; }
    }
    switch (fCmp) {
      case LT:
        return true;
      case GT:
        return false;
      default:
        switch (rCmp) {
          case LT:
            return true;
          case EQ:
            return x->tag < y->tag;
          default:
            return false;
        }
        return false;
    }
  }
};


#endif  // #ifndef SANDBOX_JVD_APPS_DMX_DMXREAD_H_









































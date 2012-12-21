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


#include "dmxRead.h"


dmxRead::dmxRead() {
  groupSize = 0;
  clusterSize = 0;
}

dmxRead * dmxRead::newClone() {
  dmxRead * clone = new dmxRead();
  clone->fId = fId;
  clone->rId = rId;
  clone->tag = tag;
  clone->fSeq = fSeq;
  clone->rSeq = rSeq;
  clone->fQual = fQual;
  clone->rQual = rQual;
  clone->fIdx = fIdx;
  clone->descriptionCode = descriptionCode;
  clone->groupSize = groupSize;
  clone->clusterSize = clusterSize;
  return clone;
}

dmxRead::dmxRead( barcodeAssignmentType _descriptionCode, std::string _tag, unsigned _fIdx, unsigned _rIdx ) {
  descriptionCode = _descriptionCode;
  tag = _tag;
  fIdx = _fIdx;
  rIdx = _rIdx;
  groupSize = 0;
  clusterSize = 0;
}

void dmxRead::fwd( std::string _fId, std::string _fSeq ) {
  fId = _fId;
  fSeq = _fSeq;
  fQual.assign( _fSeq.length(), 'J' );
}

void dmxRead::rev( std::string _rId, std::string _rSeq ) {
  rId = _rId;
  rSeq = _rSeq;
  rQual.assign( _rSeq.length(), 'J' );
}

void dmxRead::fwd( std::string _fId, std::string _fSeq, std::string _fQual ) {
  fId = _fId;
  fSeq = _fSeq;
  fQual = _fQual;
}

void dmxRead::rev( std::string _rId, std::string _rSeq, std::string _rQual ) {
  rId = _rId;
  rSeq = _rSeq;
  rQual = _rQual;
}

std::string dmxRead::getDescription() {

  switch (descriptionCode) {
    case BOTH:
      return "Concordant"; break;
    case FWD:
      return "Forward"; break;
    case REV:
      return "Reverse"; break;
    case NO_MATCH:
      return "Unidentifiable"; break;
    case MISMATCH:
      return "Discordant"; break;
    default:
      return "ERROR!"; break;
  }
}

std::string dmxRead::getShortDescription() {

  switch (descriptionCode) {
    case BOTH:
      return "CON"; break;
    case FWD:
      return "FWD"; break;
    case REV:
      return "REV"; break;
    case NO_MATCH:
      return "NON"; break;
    case MISMATCH:
      return "DIS"; break;
    default:
      return "ERR"; break;
  }
}

void dmxRead::print() {

  std::string description = getDescription();

  std::cout 
    << description << " " 
    << tag << " " 
    << fIdx << " " 
    << rIdx << " " 
    << fId << " " 
    << fSeq << " " 
    << rId << " " 
    << " groupSize " << groupSize
    << " clusterSize " << clusterSize
    << rSeq 
    << std::endl;
}

void dmxRead::printFFasta( unsigned i, std::ofstream & fh ) {

  std::string description = getShortDescription();

  fh 
    << ">"
    << description << "_" << i << "_1 " 
    << tag << " " 
    << fIdx << " " 
    << fId
    << " groupSize " << groupSize
    << " clusterSize " << clusterSize
    << std::endl 
    << fSeq 
    << std::endl;
}

void dmxRead::printRFasta( unsigned i, std::ofstream & fh ) {

  std::string description = getShortDescription();

  fh
    << ">"
    << description << "_" << i << "_2 " 
    << tag << " " 
    << rIdx << " " 
    << rId
    << " groupSize " << groupSize
    << " clusterSize " << clusterSize
    << std::endl
    << rSeq 
    << std::endl;
}

void dmxRead::printFFastq( unsigned i, std::ofstream & fh ) {

  std::string description = getShortDescription();

  fh 
    << "@"
    << description << "_" << i << "_1 " 
    << tag << " " 
    << fIdx << " " 
    << fId
    << " groupSize " << groupSize
    << " clusterSize " << clusterSize
    << std::endl 
    << fSeq
    << std::endl << "+" << std::endl
    << fQual
    << std::endl;
}

void dmxRead::printRFastq( unsigned i, std::ofstream & fh ) {

  std::string description = getShortDescription();

  fh
    << "@"
    << description << "_" << i << "_2 " 
    << tag << " " 
    << rIdx << " " 
    << rId
    << " groupSize " << groupSize
    << " clusterSize " << clusterSize
    << std::endl
    << rSeq 
    << std::endl << "+" << std::endl
    << fQual
    << std::endl;
}

void dmxRead::printFasta( unsigned i, std::ofstream & fh ) {
  printFFasta( i, fh );
  printRFasta( i, fh );
}

void dmxRead::printFastq( unsigned i, std::ofstream & fh ) {
  printFFastq( i, fh );
  printRFastq( i, fh );
}

void dmxRead::getDinucleotideFreqs( std::vector< double > & kmer ) {
  std::vector< double > fkmer, rkmer;
  getDinucleotideFreqs( fSeq, fkmer );
  kmer.insert( kmer.begin(), fkmer.begin(), fkmer.end() );
  getDinucleotideFreqs( rSeq, rkmer );
  kmer.insert( kmer.end(), rkmer.begin(), rkmer.end() );
}

void dmxRead::getDinucleotideFreqs( std::string s, std::vector< double > & kmer ) {

  std::map< std::string, int > counts;
  counts["AA"]=0;counts["AC"]=0;counts["AG"]=0;counts["AT"]=0;
  counts["CA"]=0;counts["CC"]=0;counts["CG"]=0;counts["CT"]=0;
  counts["GA"]=0;counts["GC"]=0;counts["GG"]=0;counts["GT"]=0;
  counts["TA"]=0;counts["TC"]=0;counts["TG"]=0;counts["TT"]=0;

  for ( size_t i = 1; i < s.size(); ++i ) {
    counts[ s.substr( i - 1, 1 ) ]++;
  }

  for ( std::map< std::string, int >::iterator it = counts.begin(); it != counts.end(); ++it ) {
    kmer.push_back( (double) (*it).second );
  }
}

bool dmxRead::operator== ( dmxRead & other ) {

  enum cmp { LT, EQ, GT };

  if (descriptionCode != other.descriptionCode) {
    std::cout << "Trying to compare " << descriptionCode << " to " << other.descriptionCode << " in read sort." << std::endl;
  }

  cmp fCmp = EQ;
  cmp rCmp = EQ;

  if ( descriptionCode != REV ) {
    if ( fId < other.fId ) {
      fCmp = LT;
    }
    else if ( fId > other.fId ) {
      fCmp = GT;
    }
  }
  if ( descriptionCode != FWD ) {
    if ( rId < other.rId ) {
      rCmp = LT;
    }
    else if ( rId > other.rId ) {
      rCmp = GT;
    }
  }

  switch (fCmp) {
    case LT:
      return false;
    case GT:
      return false;
    default:
      switch (rCmp) {
        case LT:
          return false;
        case GT:
          return false;
        default:
          return tag == other.tag;
      }
  }
}

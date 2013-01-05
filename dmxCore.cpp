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

#include "dmxCore.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

dmx::dmx( char* barcodeFile ) { 
  finishedReading = false;
  spoon = true;
  readBarcodeFile(barcodeFile);
}

void dmx::initFastq( unsigned _maxDistance, unsigned _chunkSize, unsigned _trimSize ) {
  maxDistance = _maxDistance;
  chunkSize = _chunkSize;
  trimSize = _trimSize;
}

void dmx::runFastq ( char* pair1FileName, char* pair2FileName ) {
  pairedEnd = true;
  fastqChunks.clear();
  fastqFeed.clear();

  dmxio = new dmxIO( pair1FileName, pair2FileName, chunkSize, 4 );

  #pragma omp parallel sections
  {
    #pragma omp section
    { parallelDigest2(); }
    #pragma omp section
    { dmxio->buffer(); }
    #pragma omp section
    { read2FilePairedFastq(pair1FileName, pair2FileName); }
  }
  std::cout << "finished reading and digesting..." << std::endl;
}

void dmx::read2FilePairedFastq( char * pair1FileName, char * pair2FileName ) {
  using namespace std;

  printf( "Begin reading Paired Fastq Files...\n" );

  int l = 0;
  string line1, line2;
  fastqPair fqp;
  bool breaker = false;

  vector< fastqPair > * chunk = new vector< fastqPair >();
  chunk->reserve( chunkSize );
  size_t itemsPerChunk = chunkSize;

  unsigned n = 0;
  unsigned n_c = 0;

  while ( !dmxio->ready() ) {
    usleep(10);
  }

  while ( !dmxio->isEmpty() ) {
    if ( !dmxio->getline( pair1FileName, line1 ) ) breaker = true;
    if ( !dmxio->getline( pair2FileName, line2 ) ) breaker = true;
    if (breaker) { break; }
    l++;
    //printf( "l = %d, n = %d, n_c = %d,  %s , %s\n", l, n, n_c, line1.c_str(), line2.c_str() );
    switch ( l ) {
      case 1:
        fqp.id1 = line1;
        fqp.id2 = line2;
        break;
      case 2:
        fqp.sq1 = line1.substr( trimSize, line1.size() );
        fqp.sq2 = line2.substr( trimSize, line2.size() );
        break;
      case 3:
        break;
      case 4:
        fqp.ql1 = line1.substr( trimSize, line1.size() );
        fqp.ql2 = line2.substr( trimSize, line2.size() );
        fqp.num = n;
        n++;
        chunk->push_back( fqp ); 
        fqp = fastqPair();
        l = 0;
        break;
      default:
        break;
    }
    if ( l == 0 && chunk->size() >= itemsPerChunk ) {
      fastqChunks.push( chunk );
      chunk = new vector< fastqPair >;
      chunk->reserve( chunkSize );
      ++n_c;
    }
  }

  if ( chunk->size() > 0 ) {
    fastqChunks.push( chunk );
  }
  finishedReading = true;
  printf( "Finished reading Paired Fastq Files. Contained %d reads, %d chunks...\n", n, n_c );
}


int dmx::readBarcodeFile(char* barcodeFileName) {
  printf( "Begin reading barcode file...\n");
  std::ifstream barcodeFile (barcodeFileName);

  while (barcodeFile) {
    std::string line;
    if (!getline( barcodeFile, line )) break;
    std::istringstream ss(line);
    std::vector <std::string> elements; 
    while (ss) {
      std::string s;
      if (!getline( ss, s, ' ' )) break;
      elements.push_back(s);
    }  
    barcodes[elements[0]] = barcode();
    barcodeNames.push_back( elements[0] );
    barcodes[elements[0]].loadBarcode(elements[1],elements[2]);
  } 
  if (!barcodeFile.eof()) {
    std::cerr << "Inconceivable!\n";
  }

  resize( barcodeStringSet, barcodes.size() );
  for ( size_t i = 0; i < barcodeNames.size(); ++i ) {
    barcodeStringSet[ i ] = barcodes[ barcodeNames[ i ] ].barcodeString;
  }
  barcodeStringSetIndex = barcodeStringSetType( barcodeStringSet );
  barcodeFinder = barcodeStringSetIndexType( barcodeStringSetIndex );
  printf( "Finished reading barcode file...\n");
  return 0;
}

dmxMatch dmx::getMatchIndexFinder( std::string seq, barcodeStringSetIndexFinderType * _barcodeFinder ) {
  dmxMatch m;
  m.min = 6;
  m.index = 0;

  clear( * _barcodeFinder );
  while ( find( * _barcodeFinder, seq.substr(27,6) ) ) {
    m.index = position( * _barcodeFinder ).i1;
    m.min = 0;
    return m;
  }
  return m;
}

void dmx::digest() {
  for ( concurrent_vector< std::vector< fastqPair > * >::iterator chunk = fastqFeed.begin();
      chunk != fastqFeed.end(); ++chunk ) {
    digest( &barcodeFinder, *chunk );
  }
}

void dmx::digest( barcodeStringSetIndexFinderType *_barcodeFinder, std::vector< fastqPair > * fastqFeedChunk ) {

  using namespace std;

  for ( std::vector< fastqPair >::iterator pairIt = (*fastqFeedChunk).begin();
      pairIt != (*fastqFeedChunk).end(); ++pairIt ) {

    string & fwdMate = (*pairIt).sq1;
    string & revMate = (*pairIt).sq2;
    string & fwdMateQual = (*pairIt).ql1;
    string & revMateQual = (*pairIt).ql2;
    int r = (*pairIt).num;

    dmxMatch fwdMatch = getMatchIndexFinder( fwdMate, _barcodeFinder );
    //dmxMatch fwdMatch = getExactMatch( fwdMate );

    unsigned fwdMin = fwdMatch.min;
    int fwdMinIndex = fwdMatch.index;

    dmxMatch revMatch = getMatchIndexFinder( revMate, _barcodeFinder );
    //dmxMatch revMatch = getExactMatch( revMate );

    unsigned revMin = revMatch.min;
    int revMinIndex = revMatch.index;

    barcodeAssignmentType BCA = NO_MATCH;
    barcode & fBC = barcodes[ barcodeNames[ fwdMinIndex ] ];
    barcode & rBC = barcodes[ barcodeNames[ revMinIndex ] ];

    if ( fwdMin > fBC.maxBarcodeDistance && revMin > rBC.maxBarcodeDistance ) {
      BCA = NO_MATCH;
      dmxRead * read = new dmxRead( NO_MATCH, "", r, r );
      read->fwd( -1, fwdMate, fwdMateQual );
      read->rev( -1, revMate, revMateQual );
      nonBarcode.push( read );
    }
    else if (fwdMinIndex == revMinIndex) { 
      if (fwdMin <= fBC.maxBarcodeDistance || 
          revMin <= rBC.maxBarcodeDistance ) {
        BCA = BOTH;
        dmxRead * read = new dmxRead(
            BOTH,
            fwdMate.substr( fBC.randTagStart, fBC.randTagLength ) +
            fwdMate.substr( fBC.randPrimerStart, fBC.randPrimerLength ) +
            fwdMate.substr( fBC.seqStart, seqTagLength ) +
            revMate.substr( rBC.randTagStart, rBC.randTagLength ) +
            revMate.substr( rBC.randPrimerStart, rBC.randPrimerLength ) +
            revMate.substr( rBC.seqStart, seqTagLength ),
            r, r );
        read->fwd( fwdMinIndex, fwdMate.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ), fwdMateQual.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ) );
        read->rev( revMinIndex, revMate.substr( rBC.seqStart, revMate.length() - rBC.seqStart ), revMateQual.substr( rBC.seqStart, revMate.length() - rBC.seqStart ) );
        conBarcode.push( read );
      }
    }
    else {
      if (fwdMin <= fBC.maxBarcodeDistance && 
          revMin > rBC.maxBarcodeDistance ) {
        BCA = FWD;
        dmxRead * read = new dmxRead(
            FWD,
            fwdMate.substr( fBC.randTagStart, fBC.randTagLength ) +
            fwdMate.substr( fBC.randPrimerStart, fBC.randPrimerLength ) +
            fwdMate.substr( fBC.seqStart, seqTagLength ), 
            r, r );
        read->fwd( fwdMinIndex, fwdMate.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ), fwdMateQual.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ) );
        read->rev( -1, revMate, revMateQual );
        fwdBarcode.push( read );
      }
      else if (fwdMin > fBC.maxBarcodeDistance && 
          revMin <= rBC.maxBarcodeDistance ) {
        BCA = REV;
        dmxRead * read = new dmxRead(
            REV,
            revMate.substr( rBC.randTagStart, rBC.randTagLength ) +
            revMate.substr( rBC.randPrimerStart, rBC.randPrimerLength ) +
            revMate.substr( rBC.seqStart, seqTagLength ),
            r, r );
        read->fwd( -1, fwdMate, fwdMateQual );
        read->rev( revMinIndex, revMate.substr( rBC.seqStart, revMate.length() - rBC.seqStart ), revMateQual.substr( rBC.seqStart, revMate.length() - rBC.seqStart ) );
        revBarcode.push( read );
      }
      else if (fwdMin <= fBC.maxBarcodeDistance && 
          revMin <= rBC.maxBarcodeDistance ) {
        BCA = MISMATCH;
        dmxRead * read = new dmxRead(
            MISMATCH,
            fwdMate.substr( fBC.randTagStart, fBC.randTagLength ) +
            fwdMate.substr( fBC.randPrimerStart, fBC.randPrimerLength ) +
            fwdMate.substr( fBC.seqStart, seqTagLength ) +
            revMate.substr( rBC.randTagStart, rBC.randTagLength ) +
            revMate.substr( rBC.randPrimerStart, rBC.randPrimerLength ) +
            revMate.substr( rBC.seqStart, seqTagLength ),
            r, r );       
        read->fwd( fwdMinIndex, fwdMate.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ), fwdMateQual.substr( fBC.seqStart, fwdMate.length() - fBC.seqStart ) );
        read->rev( revMinIndex, revMate.substr( fBC.seqStart, revMate.length() - rBC.seqStart ), revMateQual.substr( fBC.seqStart, revMate.length() - rBC.seqStart ) );
        disBarcode.push( read );
      }
    }
  }
}

void dmx::parallelDigest() {
  digestFunctor df;
  df.d = this;

  std::vector< fastqPair > * chunk;

  while ( !(fastqChunks.empty()) ) { 
    if ( fastqChunks.try_pop( chunk ) ) {
      fastqFeed.push_back( chunk );
      chunk = new std::vector< fastqPair >();
    }
  }

  if ( fastqFeed.size() > 4 ) {
    parallel_for( blocked_range< int >( 0, fastqFeed.size() ), df );
  }
  else {
    for ( size_t i = 0; i < fastqFeed.size(); ++i ) {
      digest( & barcodeFinder, fastqFeed[ i ] );
    }
  }
}

void dmx::parallelDigest2() {
  digestFunctor2 df;
  df.d = this;

  // wait for some work to build up...
  while ( fastqChunks.empty() ) {
    usleep( 10 ); // TODO: CLEANUP
  }

  printf( "Digesting...\n" );

  fastqFeed.reserve( fastqChunks.unsafe_size() * 2 );

  int numChunks = fastqChunks.unsafe_size();
  for ( int i = 0; i < numChunks; ++i ) {
    std::vector< fastqPair > * chunk;
    fastqChunks.try_pop( chunk );
    fastqFeed.push_back( chunk );
  }

  parallel_do( fastqFeed.begin(), fastqFeed.end(), df );
  while ( !finishedReading ) { 
    parallel_do( fastqFeed.begin(), fastqFeed.end(), df );
  } // TODO: clean this up ...
  fastqFeed.clear();
  fastqFeed.shrink_to_fit();
  fastqChunks.clear();
  spoon = true;

  // TODO these should be elsewhere...
  convertPriorityQueuesToVectors();
  groupReduce();
  convertPriorityQueuesToVectors();

}

void dmx::convertPriorityQueueToVector( dmxReadPriQ & q, dmxReadSerialVector & v ) {

  v.clear();
  dmxRead * r;

  while ( q.try_pop(r) ) {
    v.push_back( r );
  }

}

void dmx::convertPriorityQueuesToVectors() {
  // popping from priority queues is slow and inherently serial
  // copy everything into vectors in preparation for subsequent 
  // parallel processing.  Seems to be significantly faster if
  // the copying is done in parallel
#pragma omp parallel sections
  {
#pragma omp section
    {
      convertPriorityQueueToVector( fwdBarcode, fwdBarcodeSerVec );  
      printf( "FWD %lu\n", fwdBarcodeSerVec.size() );
    }
#pragma omp section
    {
      convertPriorityQueueToVector( revBarcode, revBarcodeSerVec );  
      printf( "REV %lu\n", fwdBarcodeSerVec.size() );

    }
#pragma omp section
    {
      convertPriorityQueueToVector( conBarcode, conBarcodeSerVec );  
      printf( "CON %lu\n", fwdBarcodeSerVec.size() );

    }
#pragma omp section
    {
      convertPriorityQueueToVector( disBarcode, disBarcodeSerVec );  
      printf( "DIS %lu\n", fwdBarcodeSerVec.size() );

    }
#pragma omp section
    {
      convertPriorityQueueToVector( nonBarcode, nonBarcodeSerVec );  
      printf( "NON %lu\n", fwdBarcodeSerVec.size() );

    }
  }
}

void dmx::groupReduce() {
  std::cout << "group reduce" << std::endl;
#pragma omp parallel sections
  {
#pragma omp section
    { groupReduce( &fwdBarcodeSerVec, &fwdBarcode ); }
#pragma omp section
    { groupReduce( &revBarcodeSerVec, &revBarcode ); }
#pragma omp section
    { groupReduce( &conBarcodeSerVec, &conBarcode ); }
  }
}

void dmx::groupReduce( dmxReadSerialVector * drsv, dmxReadPriQ * drpq ) {
  // runs through barcode vectors and groups, then clusters, then reduces
  // each stack to a consensus, inserts into priority queues...
  // uses the 'spoon feeding' parallel_do approach...
  tbb::atomic< bool > * my_spoon = new tbb::atomic< bool >();
  (*my_spoon) = true;
  tbb::concurrent_vector< dmxReadSerialVector * > groupReduceFeed; 

  size_t bufferCount = 0;
  size_t bufferSize = 10;

  dmxReadSerialVector * stack = new dmxReadSerialVector();
  
  for ( dmxReadSerialVector::iterator it = drsv->begin(); it != drsv->end(); ++it ) {
    if ( stack->empty() ) {
      stack->push_back( *it );
    }
    else {
      if ( *(*it) == *( stack->back() ) ) {
        stack->push_back( *it );
      }
      else {
        groupReduceFeed.push_back( stack );
        stack = new dmxReadSerialVector();
        bufferCount++;
      }
    }
    if ( bufferCount == bufferSize || bufferCount == drsv->size() ) {
      break;
    }
  }
  parallel_do( groupReduceFeed.begin(), groupReduceFeed.end(), dmx::groupReduceFunctor( this, drsv, drpq, my_spoon ) );
  drsv->clear();
}

inline void dmx::groupReduceFunctor::operator() ( dmxReadSerialVector * r, parallel_do_feeder< dmxReadSerialVector * >& feeder ) const {
  
  if ( group_reduce_spoon->compare_and_swap( false, true ) ) {
    // I have the group_reduce_spoon, so I must be the feeder...
    // following loop only entered by the feeder
    dmxReadSerialVector * stack = new dmxReadSerialVector();

    for ( dmxReadSerialVector::iterator it = drsv->begin(); it != drsv->end(); ++it ) {
      if ( stack->empty() ) {
        stack->push_back( *it );
      }
      else {
        if ( *(*it) == *( stack->back() ) ) {
          stack->push_back( *it );
        }
        else {
          feeder.add( stack );
          stack = new dmxReadSerialVector();
        }
      }
    }


  }

  // executed by all tasks, including, eventually, the feeder:
  if ( r->size() > 0 ) {

    std::map< int, dmxReadSerialVector > clusterMap;
    // TODO parameterize the minimumum size to condense; randomly select from vector when below size

    dmxRead * processedRead;

    if ( r->size() > 100 ) {
      d->getClusters( clusterMap, r );

      for ( std::map< int, dmxReadSerialVector >::iterator it = clusterMap.begin(); it != clusterMap.end(); ++it ) {
        if ( (*it).second.size() > 10 ) {
          processedRead = d->condenseGroup( (*it).second );
        }
        else {
          processedRead = (*it).second.front()->newClone() ;
        }
      }
    }
    else {
      processedRead = r->front()->newClone();
    }

    processedRead->groupSize = r->size();
    //std::cout << &it << " numclusters: " << clusterMap.size() << " groupSize: " << r->size() << " clusterSize: " << (*it).second.size() << std::endl;
    drpq->push( processedRead );
    for ( dmxReadSerialVector::iterator i = r->begin(); i != r->end(); ++i ) {
      (*i) = NULL;
    }
  }
  r->clear();
}




void dmx::getClusters( std::map< int, dmxReadSerialVector > & clusterMap, dmxReadSerialVector * rv ) {
  // clusters reads then loads cluster membership into map passed as reference
  std::vector< double > groupKmers, readKmers;

  for ( dmxReadSerialVector::iterator it = rv->begin(); it != rv->end(); ++it ) {
    (*it)->getDinucleotideFreqs( readKmers );
    groupKmers.insert( groupKmers.end(), readKmers.begin(), readKmers.end() );
    readKmers.clear();
  }
  double * groupKmerArray = &groupKmers[0];
  lti::matrix<double> groupKmerMatrix( rv->size(), 32, groupKmerArray );

  typedef lti::l2Distantor< lti::vector< double > > distanceType;
  lti::DBScan< distanceType >::parameters clusteringParameters;
  clusteringParameters.eps = 20;
  clusteringParameters.minPts = 1; // if greater than '1', points can be classified as noise and placed in cluster '0'
  lti::DBScan< distanceType > clustering;
  clustering.setParameters( clusteringParameters );
  lti::ivector clusteringResult;
/*
  lti::adaptiveKMeans::parameters clusteringParameters;
  clusteringParameters.detectNeighborhood = true;
  lti::adaptiveKMeans clustering;
  clustering.setParameters( clusteringParameters );
  lti::ivector clusteringResult; 
*/  
  clustering.train( groupKmerMatrix, clusteringResult );

  int i = 0;
  for ( lti::ivector::iterator it = clusteringResult.begin(); it != clusteringResult.end(); ++it ) {
    // exclude noise (cluster==0)
    if ( (*it) > 0 ) {
      clusterMap[ (*it) ].push_back( rv->at( i ) ); 
    }
    i++;
  }
}

dmxRead * dmx::condenseGroup( std::vector< dmxRead * > & rv ) {
  typedef String< Dna5 > TSequence;
  StringSet< TSequence > fSeq;
  StringSet< TSequence > rSeq;

  for ( std::vector< dmxRead * >::iterator i = rv.begin(); i != rv.end(); ++i ) {
    appendValue( fSeq, (*i)->fSeq );
    appendValue( rSeq, (*i)->rSeq );
  }

  Graph< Alignment< StringSet< TSequence, Dependent<> > > > fAliG( fSeq );
  Graph< Alignment< StringSet< TSequence, Dependent<> > > > rAliG( rSeq );

  globalMsaAlignment( fAliG, Score<int>(0, -1, -1, -2) );
  globalMsaAlignment( rAliG, Score<int>(0, -1, -1, -2) );

  std::string fMatrix;
  std::string rMatrix;

  convertAlignment( fAliG, fMatrix );
  convertAlignment( rAliG, rMatrix );

  std::string fCon = computeConsensus( fMatrix, rv.size() );
  std::string rCon = computeConsensus( rMatrix, rv.size() );

  // TODO modify dmxRead struct to record consensus info (reads that go into consensus, etc.)... halfway done...
  dmxRead * r = new dmxRead( rv.front()->descriptionCode, rv.front()->tag, rv.front()->fIdx, rv.front()->rIdx ); 
  r->fwd( rv.front()->getFwdBCidx(), fCon );
  r->rev( rv.front()->getRevBCidx(), rCon );
  r->clusterSize = rv.size(); 
  return r;
}

std::string dmx::computeConsensus( std::string & matrix, size_t nrow ) {
  size_t skip = matrix.size() / nrow;
  std::string consensus = "";

  std::map< char, int > baseCount;
  std::map< int, char > countBase;
  
  for ( size_t i = 0; i < skip; ++i ) {

    baseCount.clear();
    countBase.clear();
    
    for ( size_t j = nrow - 1; j != 0; --j ) {
      size_t index = ( j * skip ) + i;
      baseCount[ matrix.at( index ) ]++;
      countBase.insert( std::pair< int, char >( baseCount[ matrix.at( index ) ], matrix.at( index ) ) );
    }
 
    if ( ( * ( countBase.rbegin() ) ).second != '-' ) {
      consensus.push_back( ( * ( countBase.rbegin() ) ).second );
    }

  }
  return consensus;
}





/////////////// PRINTING /////////////////////

void dmx::printBarcodeResults( concurrent_vector< dmxRead * > resultVector ) {
  for ( concurrent_vector< dmxRead * >::iterator i = resultVector.begin(); i != resultVector.end(); ++i ) {
    (*i)->print();
  }
}

void dmx::printBarcodeResults( dmxReadPriQ &resultVector ) {
  dmxRead * read;
  while (resultVector.try_pop(read)) {
    read->print();
  }
}

void dmx::printConcordantBarcodeResults() {
  printBarcodeResults( conBarcode );
}

void dmx::printFwdOnlyBarcodeResults() {
  printBarcodeResults( fwdBarcode );
}

void dmx::printRevOnlyBarcodeResults() {
  printBarcodeResults( revBarcode );
}

void dmx::printDiscordantBarcodeResults() {
  printBarcodeResults( disBarcode );
}

void dmx::printUnidentifiableBarcodeResults() {
  printBarcodeResults( nonBarcode );
}

void dmx::printFasta( dmxReadSerialVector drv, std::ofstream & fh ) {
  for ( unsigned i = 0; i < drv.size(); ++i ) {
    drv[ i ]->printFasta( i, fh );
  }
}

void dmx::printFastq( dmxReadSerialVector drv, std::ofstream & fh ) {
  for ( unsigned i = 0; i < drv.size(); ++i ) {
    drv[ i ]->printFastq( i, fh );
  }
}

void dmx::printFasta( dmxReadSerialVector drv, std::ofstream & fh, int barcode_index ) {
  for ( unsigned i = 0; i < drv.size(); ++i ) {
    if ( drv[ i ]->getFwdBCidx() == barcode_index ) {
      drv[ i ]->printFFasta( i, fh );
    }
    if ( drv[ i ]->getRevBCidx() == barcode_index ) {
      drv[ i ]->printRFasta( i, fh );
    }
  }
}

void dmx::printFastq( dmxReadSerialVector drv, std::ofstream & fh, int barcode_index ) {
  for ( unsigned i = 0; i < drv.size(); ++i ) {
    if ( drv[ i ]->getFwdBCidx() == barcode_index ) {
      drv[ i ]->printFFastq( i, fh );
    }
    if ( drv[ i ]->getRevBCidx() == barcode_index ) {
      drv[ i ]->printRFastq( i, fh );
    }
  }
}

void dmx::printFasta( dmxReadPriQ & drq, std::ofstream & fh ) {
  dmxRead * read;
  unsigned i = 0;
  while ( drq.try_pop(read) ) {
    read->printFasta( i, fh );
    ++i;
  }
}

void dmx::printFasta( dmxReadPriQ & drq, std::ofstream & fh, int barcode_index ) {
  dmxRead * read;
  unsigned i = 0;
  while ( drq.try_pop(read) ) {
    if ( read->getFwdBCidx() == barcode_index ) {
      read->printFFasta( i, fh );
    }
    if ( read->getRevBCidx() == barcode_index ) {
      read->printRFasta( i, fh );
    }
    ++i;
  }
}

void dmx::printGoodFasta( std::string goodFastaOutfile ) {
  std::ofstream fh ( goodFastaOutfile.c_str() );
  printFasta( conBarcodeSerVec, fh );
  printFasta( fwdBarcodeSerVec, fh );
  printFasta( revBarcodeSerVec, fh );
  fh.close();
}

void dmx::printGoodFastq( std::string goodFastqOutfile ) {
  std::ofstream fh ( goodFastqOutfile.c_str() );
  printFastq( conBarcodeSerVec, fh );
  printFastq( fwdBarcodeSerVec, fh );
  printFastq( revBarcodeSerVec, fh );
  fh.close();
}

void dmx::printAllFasta( std::string goodFastaOutfile ) {
  std::ofstream fh ( goodFastaOutfile.c_str() );
  printFasta( conBarcodeSerVec, fh );
  printFasta( fwdBarcodeSerVec, fh );
  printFasta( revBarcodeSerVec, fh );
  printFasta( disBarcode, fh );
  printFasta( nonBarcode, fh );
  fh.close();
}












////////////// JUNKYARD /////////////////////////////////////////////

unsigned int dmx::distance(const std::string s1, const std::string s2) {
  using namespace std;
  const size_t len1 = s1.size(), len2 = s2.size();
  vector<unsigned int> col(len2+1), prevCol(len2+1);

  for (unsigned int i = 0; i < prevCol.size(); i++)
    prevCol[i] = i;
  for (unsigned int i = 0; i < len1; i++) {
    col[0] = i+1;
    for (unsigned int j = 0; j < len2; j++) {
      col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]),
          prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
    }
    col.swap(prevCol);
    if ( col[i] > maxDistance ) {
      return col[i];
    }
  }
  return prevCol[len2];
}

dmxMatch dmx::getMatch( std::string seq ) {

  dmxMatch m;
  m.min = seq.size();
  m.index = -1;
  unsigned dist;

  for (unsigned i = 0; i < barcodeNames.size(); ++i) {

    barcode & b = barcodes[ barcodeNames[ i ] ];

    dist = distance( seq.substr(b.barcodeStart,b.barcodeLength), b.barcodeString );
    if ( dist < m.min ) {
      m.min = dist;
      m.index = i;
      if ( dist < b.maxBarcodeDistance ) {
        break;
      }
    }
  }
  return m;
}

dmxMatch dmx::getExactMatch( std::string seq ) {

  dmxMatch m;
  m.min = seq.size();
  m.index = -1;
  unsigned dist;

  for (unsigned i = 0; i < barcodeNames.size(); ++i) {

    barcode & b = barcodes[ barcodeNames[ i ] ];

    //dist = distance( seq.substr(b.barcodeStart,b.barcodeLength), b.barcodeString );

    if ( seq.substr(b.barcodeStart,b.barcodeLength) == b.barcodeString ) {
      dist = 0;
    }

    //printf( "%s %s\n",  seq.substr(b.barcodeStart,b.barcodeLength).c_str(), b.barcodeString.c_str() );
    if ( dist < m.min ) {
      m.min = dist;
      m.index = i;
      if ( dist < b.maxBarcodeDistance ) {
        break;
      }
    }
  }
  return m;
}

dmxMatch dmx::getMatchMyersInfix( std::string seq ) {
  dmxMatch m;
  m.min = seq.size();
  m.index = -1;

  for (unsigned i = 0; i < barcodeNames.size(); ++i) {

    barcode b = barcodes[ barcodeNames[ i ] ];

    String< char > haystack = seq.substr(b.barcodeStart,b.barcodeLength);
    Finder< String< char > > finder( haystack );
    String< char > needle = b.barcodeString;
    Pattern< String< char >, Myers< FindInfix > > pattern( needle );

    int min = b.barcodeLength;
    while (find(finder, pattern, -b.maxBarcodeDistance)) {
      //std::cout << "end: " << endPosition(finder) << std::endl;
      while (findBegin(finder, pattern, getScore(pattern))) {
        //std::cout << "begin: " << beginPosition(finder) << std::endl;
        //std::cout << b.barcodeString << std::endl;
        //std::cout << infix(finder) << " matches with score ";
        if ( -getBeginScore(pattern) < min) {
          min = -getBeginScore(pattern);
        }
        //std::cout << (unsigned) min << std::endl;
      }
    }
    if ( (unsigned) min < m.min ) {
      m.min = (unsigned) min;
      m.index = i;
      if (min == 0) {
        return m;
      }
    }
  }
  return m;
}

void dmx::test_consensus() {
    typedef String< Dna5 > TSequence;
    StringSet< TSequence > seq;
    appendValue(seq,"AAACTCCGAT");
    appendValue(seq,"AAAGTCCGAT");
    appendValue(seq,"AAACTCCGA");
    appendValue(seq,"AAAGTCCG");
    appendValue(seq,"AAAGTCCG");
    Graph< Alignment< StringSet< TSequence, Dependent<> > > > aliG( seq );
    //globalMsaAlignment(aliG, Blosum62(-1, -11) );
    globalMsaAlignment( aliG, Score<int>(0, -1, -1, -2) );
    //::std::cout << aliG << ::std::endl;
    std::string matrix;
    convertAlignment( aliG, matrix );
    std::cout << aliG << std::endl;
    std::cout << matrix << std::endl;
    std::cout << computeConsensus( matrix, 5 ) << std::endl;
}

void dmx::cluster_test() {

  using namespace lti;

  double rawData[] = {
    1,1,0,0,0,0,0,0,1,1,
    1,1,0,0,0,0,0,0,1,1,
    1,1,0,0,0,0,0,0,1,1,
    1,1,0,0,0,0,0,0,0,0,
    1,1,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,1,1,1,0,0,
    0,0,1,1,1,0,0,0,0,0
  };

  lti::matrix<double> inData( 10, 10, rawData );

  typedef l2Distantor< lti::vector< double > > distanceType;

  lti::DBScan< distanceType >::parameters clusteringParm;
  //adaptiveKMeans::parameters clusteringParm;

  clusteringParm.eps = 1;
  clusteringParm.minPts = 1;

  lti::DBScan< distanceType > clustering;

  //adaptiveKMeans clustering;

  clustering.setParameters(clusteringParm);

  lti::ivector clusteringResult;

  clustering.train( inData, clusteringResult );

  std::cerr << clusteringResult << std::endl;

}









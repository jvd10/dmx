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

#ifndef SANDBOX_JVD_APPS_DMX_DMXIO_H_
#define SANDBOX_JVD_APPS_DMX_DMXIO_H_


#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/index.h>

#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_do.h>
#include <tbb/parallel_sort.h>

#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <deque>
#include <list>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


typedef tbb::spin_mutex dmxIOBufferMutexT;


struct dmxIOBuffer {

  dmxIOBufferMutexT dmxIOBufferMutex;
  dmxIOBufferMutexT dmxIOGetLineMutex;

  size_t bufferSize, refreshSize;
  std::deque < std::string > buffer;
  std::ifstream file;
  tbb::atomic< bool > fileEmpty;

  dmxIOBuffer( size_t chunkSize, size_t bufferFactor, char * filename ); 
  ~dmxIOBuffer();
  void refresh();
  bool fileIsEmpty();
  bool isEmpty();
  bool getline( std::string & line );

  void pop_buffer();

  size_t size() {
    return buffer.size();
  }

  std::ifstream fileStream;
  boost::iostreams::filtering_istream gzStream;  
};

typedef std::map< char * , dmxIOBuffer * > bufferMap;

class dmxIO {

  public:

    dmxIO( char * fileName1, char * fileName2, size_t chunkSize, size_t bufferFactor );

    bufferMap buffers;

    bool getline( char * fileName, std::string & line );
    
    void buffer();

    bool isEmpty();

    bool ready() { return ready_flag; }

  private:

    tbb::atomic< bool > ready_flag;
};


#endif  // #ifndef SANDBOX_JVD_APPS_DMX_DMXIO_H_









































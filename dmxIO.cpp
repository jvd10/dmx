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

#include "dmxIO.h"
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <tbb/tbb.h>

////////// dmxIOBuffer //////////////



dmxIOBuffer::dmxIOBuffer( size_t chunkSize, size_t bufferFactor, char * filename ) {
  refreshSize = chunkSize * 4 * bufferFactor;
  bufferSize = refreshSize; 
  buffer.clear();
  
  //file.open( filename );
  
  fileStream.open(filename, std::ios_base::in | std::ios_base::binary);
  gzStream.push( boost::iostreams::gzip_decompressor() );
  gzStream.push( fileStream );

  fileEmpty = false;
}

dmxIOBuffer::~dmxIOBuffer() {
  file.close();
}

void dmxIOBuffer::refresh() {

  using namespace std;

  dmxIOBufferMutexT::scoped_lock lock(dmxIOBufferMutex);
  std::string line;

  //printf( "current buffer size: %lu requested buffer size %lu\n", buffer.size(), bufferSize );

  if ( !fileEmpty ) {
    //printf( "file is empty?\n");
    for (size_t i = buffer.size(); i < bufferSize; ++i) {
      //if ( std::getline( file, line ) ) {
      if ( std::getline( gzStream, line ) ) {
        //std::cout << line << std::endl;
        buffer.push_back( line );
        //printf( "wtf %lu\n", buffer.size() );
      }
      //else if ( !file.eof() ) {
      //  cout << "problem reading file!" << endl;
      //  fileEmpty = true;
      //  break;
      //}
      else {
        fileEmpty = true;
        break;
      }
    }
  }
}

bool dmxIOBuffer::fileIsEmpty() {
  dmxIOBufferMutexT::scoped_lock lock(dmxIOBufferMutex);
  return fileEmpty;
}

bool dmxIOBuffer::isEmpty() {
  dmxIOBufferMutexT::scoped_lock lock(dmxIOBufferMutex);
  if ( fileEmpty ) {
    if ( buffer.size() == 0 ) {
      return true;
    }
  }
  return false;
}

void dmxIOBuffer::pop_buffer() {
  dmxIOBufferMutexT::scoped_lock lock(dmxIOBufferMutex);
  buffer.pop_front();
}

bool dmxIOBuffer::getline( std::string & line ) {
  dmxIOBufferMutexT::scoped_lock lock(dmxIOGetLineMutex);
  if ( buffer.size() == 0 ) {
    refresh();
  }
  
  if ( buffer.size() > 0 ) {
    line = buffer.front();
    pop_buffer();
    return true;
  }
  return false;
}

//////////// dmxIO ////////////////////


dmxIO::dmxIO( char * fileName1, char * fileName2, size_t chunkSize, size_t bufferFactor ) {

  buffers.clear();
  buffers[ fileName1 ] = new dmxIOBuffer( chunkSize, bufferFactor, fileName1 );
  buffers[ fileName2 ] = new dmxIOBuffer( chunkSize, bufferFactor, fileName2 );

  ready_flag = true;
}

void dmxIO::buffer() {
  while ( !( isEmpty() ) ) {
    for ( std::map< char * , dmxIOBuffer * >::iterator bit = buffers.begin(); bit != buffers.end(); ++bit ) {
      (*bit).second->refresh();
    }
  }
}

bool dmxIO::getline( char * fileName, std::string & line ) {
  if ( buffers.find( fileName ) != buffers.end() ) {
    return buffers[ fileName ]->getline( line );
  }
  else {
    printf( "Buffer not found!\n" );
    return false;
  }
}

bool dmxIO::isEmpty() {
 
  for ( bufferMap::iterator itr = buffers.begin(); itr != buffers.end(); ++itr ) {
    if ( (*itr).second->isEmpty() == false ) {
      return false;
    }
  } 
  return true;
}


//////////// PRINTING /////////////////





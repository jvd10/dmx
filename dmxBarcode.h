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

#ifndef SANDBOX_JVD_APPS_DMX_DMXBARCODE_H_
#define SANDBOX_JVD_APPS_DMX_DMXBARCODE_H_

#include <string>

struct barcode {

  std::string layout;
  std::string barcodeString, ampPrimerString, randPrimerString, randTagString;
  std::string barcodeStringRC, ampPrimerStringRC;

  unsigned digestionCutSite;

  unsigned barcodeStart, ampPrimerStart, randTagStart, randPrimerStart, seqStart;

  unsigned barcodeLength, ampPrimerLength, randTagLength, randPrimerLength;

  unsigned maxBarcodeDistance;

  void loadBarcode(std::string _layout, std::string sequence);

  inline std::string reverseComplement(const std::string s); 

  void parseBarcodeLayout();

  bool valid_type( char c );

  bool valid_value( char c );

  void set_value( char type, std::string & value, int & current_index );

  void print();
};



#endif  // #ifndef SANDBOX_JVD_APPS_DMX_DMXBARCODE_H_









































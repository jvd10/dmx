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


#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include <seqan/file.h>
#include <iostream>

#include "dmx.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
  // Setup command line parser.
  CommandLineParser parser;
  Options options;
  setupCommandLineParser(parser, options);

  // parse command line and handle the cases where help requested
  // or bad parameters were given.
  int ret = parseCommandLineAndCheck(options, parser, argc, argv);

  if (ret != 0)
  {
    std::cerr << "\n!!!!!!Invalid or missing arguments!!!!!!\n\nThese are the arguments that were given:" << std::endl;
    ret = mainWithOptions(options);
    return ret;
  }

  if (options.showHelp || options.showVersion)
    return 0;

  // Finally, launch the program.
  ret = mainWithOptions(options);

  dmx * d;

  try {
    d = new dmx( toCString(options.barcodeFile) ); 
  } catch (int e) {
    std::cout << "Unable to open one or more input files" << std::endl;
    return 1;
  }

  //d->cluster_test();
  //d->test_consensus();

  d->initFastq( 2, options.chunkSize, options.trimSize );
  d->runFastq( toCString(options.inputFiles[0]), toCString(options.inputFiles[1]) );
  
  //d->digest(0, d->readCount);
  //d->parallelDigest();
  //d->printConcordantBarcodeResults();
  //d->printFwdOnlyBarcodeResults();
  //d->printRevOnlyBarcodeResults();
  //d->printDiscordantBarcodeResults();
  //d->printUnidentifiableBarcodeResults();
  
  
  std::string outputPrefix(toCString(options.outputPrefix));
  std::string goodFastaOutfile = outputPrefix+".good.fasta";
  std::string goodFastqOutfile = outputPrefix+".good.interleaved.fastq";

  // everything in one file
  d->printGoodFastq(goodFastqOutfile);
  
  // separate files for each barcode
  //d->printPerBarcodeFasta( outputPrefix );

  return ret;
}




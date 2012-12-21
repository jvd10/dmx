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

#ifndef SANDBOX_JVD_APPS_DMX_DMX_H_
#define SANDBOX_JVD_APPS_DMX_DMX_H_

#include <stdio.h>
#include <time.h>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include "dmxCore.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Options
{
  bool showHelp;
  bool showVersion;
 
  bool pairedEnd, combinedPairs, sortedPairs;
  CharString barcodeFile;
  CharString outputPrefix;
  int chunkSize, trimSize;

  String<CharString> inputFiles;

  Options()
  {
    // Set defaults.
    showHelp = false;
    showVersion = false;
    pairedEnd = false;
    combinedPairs = false;
    sortedPairs = false;
    std::ostringstream oss;
    oss << "DMX_OUTPUT_" << time(NULL);
    outputPrefix = oss.str();
    chunkSize = 10000;
    trimSize = 0;
  }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

  void setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
  addVersionLine(parser, "0.1");

  addTitleLine(parser, "**********************");
  addTitleLine(parser, "********* dmx ********");
  addTitleLine(parser, "**********************");
  addTitleLine(parser, "");
  addTitleLine(parser, "(c) 2012 by Jay DePasse <jvd10@pitt.edu>");

  addUsageLine(parser, "[pcst] -b <barcode file> <fastq input files>");

  addSection(parser, "Barcode File Format:");
  addHelpLine(parser, "<barcode id> <barcode layout string> <barcode sequence>");

  addSection(parser, "Barcode Layout String (similar to CIGAR format; type character, followed by length):");
  addHelpLine(parser, "P = final amplification primer sequence");
  addHelpLine(parser, "R = random primer ID sequence");
  addHelpLine(parser, "B = barcode sequence");
  addHelpLine(parser, "N = random SISPA primer");
  addHelpLine(parser,"");
  addHelpLine(parser, "EXAMPLE: \"P16R4B6N6\"  <-- 16 base amplification primer, 4 base random primer ID, 6 base barcode, 6 base random SISPA primer");

  addSection(parser, "Options:");
  addOption(parser, CommandLineOption("p",  "paired", "Files contain (some) paired-end reads.", OptionType::Boolean));
  addOption(parser, CommandLineOption("c",  "combined", "Paired-end reads contained in a single file.", OptionType::Boolean));
  addOption(parser, CommandLineOption("s",  "sorted", "Paired-end reads are in sorted order.", OptionType::Boolean));
  addOption(parser, CommandLineOption("k",  "chunk", "Paired-end reads are in sorted order.", OptionType::Integer));
  addOption(parser, CommandLineOption("t",  "trim", "Number of bases to trim from beginngin of all reads before barcode search.", OptionType::Integer));
  addOption(parser, CommandLineOption("o",  "outputPrefix", "Prefix for all output files.", OptionType::String, options.outputPrefix));
  addOption(parser, CommandLineOption("b",  "barcodeFile", "Mandatory barcode file.", OptionType::String | OptionType::Mandatory));

  requiredArguments(parser, 1);
}

int parseCommandLineAndCheck(Options & options,
    CommandLineParser & parser,
    int argc,
    char const ** argv)
{

  int ret = !(parse(parser, argc, argv));

  if (ret) 
  {
    if (isSetLong(parser, "help"))
    {
      options.showHelp = true;
      ret = 0;
    }
    if (isSetLong(parser, "version"))
    {
      options.showVersion = true;
      ret = 0;
    }
  }

  getOptionValueLong(parser, "paired", options.pairedEnd);
  getOptionValueLong(parser, "combined", options.combinedPairs);
  getOptionValueLong(parser, "sorted", options.sortedPairs);
  getOptionValueLong(parser, "outputPrefix", options.outputPrefix);
  getOptionValueLong(parser, "barcodeFile", options.barcodeFile);
  getOptionValueLong(parser, "chunk", options.chunkSize);
  getOptionValueLong(parser, "trim", options.trimSize);


  options.inputFiles = getArgumentValues(parser);

  return ret;
}

int mainWithOptions(Options & options)
{
  typedef Iterator<String<CharString> >::Type TIterator;
  std::cout << "\nOptional Arguments:" << std::endl;
  std::cout << "  paired:          \"" << options.pairedEnd << "\"" << std::endl;
  std::cout << "  combined:        \"" << options.combinedPairs << "\"" << std::endl;
  std::cout << "  sorted:          \"" << options.sortedPairs << "\"" << std::endl;
  std::cout << "  output prefix:   \"" << options.outputPrefix << "\"" << std::endl;
  std::cout << "  chunk size:      \"" << options.chunkSize << "\"" << std::endl;
  std::cout << "  trim size:       \"" << options.trimSize << "\"" << std::endl;

  std::cout << "\nRequired Arguments:" << std::endl;

  std::cout << "  barcode file:     \"" << options.barcodeFile << "\"" << std::endl;
  std::cout << "  input files:    ";
  for (TIterator it = begin(options.inputFiles); it != end(options.inputFiles); ++it)
  {
    std::cout << "  \"" << *it << "\"";
  }
  std::cout << std::endl << std::endl;;

  return 0;
}



#endif  // #ifndef SANDBOX_JVD_APPS_DMX_DMX_H_

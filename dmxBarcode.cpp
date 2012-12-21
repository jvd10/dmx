#include "dmxBarcode.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sstream>

void barcode::loadBarcode(std::string _layout, std::string sequence) {
  layout = _layout;

/*
  barcodeStart = 29; // 2
  barcodeLength = 6; // 20
  ampPrimerStart = 1; // 1
  ampPrimerLength = 21; // 1
  randTagStart = 22; // 2
  randTagLength = 6; // 1
  randPrimerStart = 36; // 23
  randPrimerLength = 6; // 6
  seqStart = 43; // 29
*/

/*
  ampPrimerStart = 0; // 1
  ampPrimerLength = 21; // 1
  
  randTagStart = 21; // 2
  randTagLength = 6; // 1
 
  barcodeStart = 27; // 2
  barcodeLength = 6; // 20

  randPrimerStart = 33; // 23
  randPrimerLength = 6; // 6
  
  seqStart = 39; // 29
*/
  parseBarcodeLayout();

  barcodeString = sequence.substr(barcodeStart, barcodeLength);
  barcodeStringRC = reverseComplement(barcodeString);
  ampPrimerString = sequence.substr(ampPrimerStart, ampPrimerLength);
  randPrimerString = sequence.substr( randPrimerStart, randPrimerLength );
  randTagString = sequence.substr( randTagStart, randTagLength );
  maxBarcodeDistance = 1;

  print();
}

void barcode::print() {
  printf( "ampPrimerStart %d ampPrimerLength %d ampPrimerString %s\n", ampPrimerStart, ampPrimerLength, ampPrimerString.c_str() );
  
  printf( "randTagStart %d randTagLength %d randTagString %s\n", randTagStart, randTagLength, randTagString.c_str() );

  printf( "barcodeStart %d barcodeLength %d barcodeString %s\n", barcodeStart, barcodeLength, barcodeString.c_str() );

  printf( "randPrimerStart %d randPrimerLength %d randPrimerString %s\n\n", randPrimerStart, randPrimerLength, randPrimerString.c_str() );

}

std::string barcode::reverseComplement(const std::string s) {
  std::string r = "";
  for (int i = s.size() - 1; i >= 0; i--) {
    switch (s[i]) {
      case 'A':
        r += 'T'; break;
      case 'C':
        r += 'G'; break;
      case 'G':
        r += 'C'; break;
      case 'T':
        r += 'A'; break;
      default:
        r += s[i]; break;
    }
  }
  return r;
}

void barcode::parseBarcodeLayout() {
  int current_index = 0;
  char current_type = '\0';
  std::string current_value;

  std::string::iterator itr = layout.begin();
  for ( ; ; ++itr ) {

    if ( valid_type( *itr ) || itr == layout.end() ) {
      if ( current_type != '\0' ) {
        set_value( current_type, current_value, current_index );
      }
      else if ( itr != layout.begin() ) {
        std::cerr << "help!!!" << std::endl;
      }
      if ( itr == layout.end() ) {
        break;
      }
      current_type = *itr;
      current_value.clear();
    }
    else if ( itr != layout.end() && valid_value( *itr ) ) {
      current_value.append( 1, *itr );
    }
    else {
      std::cerr << "Invalid Barcode Layout: " << layout << " at character: " << *itr << std::endl;
      std::exit( 1 );
    }
  }

}

void barcode::set_value( char type, std::string & value, int & current_index ) {
  if ( type == '\0' || value.empty() ) {
    std::cerr << "Invalid Barcode Layout Value: " << layout << " type " << type << " value " << value << std::endl;
    std::exit( 1 );
  }
  std::istringstream ss( value );
  switch (type) {
    case 'P':
      ampPrimerStart = current_index;
      ss >> ampPrimerLength;
      current_index += ampPrimerLength;
      break;
    case 'R':
      randTagStart = current_index;
      ss >> randTagLength;
      current_index += randTagLength;
      break;
    case 'B':
      barcodeStart = current_index;
      ss >> barcodeLength;
      current_index += barcodeLength;
      break;
    case 'N':
      randPrimerStart = current_index;
      ss >> randPrimerLength;
      current_index += randPrimerLength;
      seqStart = current_index;
      break;
    default:
      std::cerr << "Unknown Barcode Feature Type! " << type << std::endl;
      std::exit( 1 );
  }
}

bool barcode::valid_type( char c ) {
  switch (c) {
    case 'P':
      break;
    case 'R':
      break;
    case 'B':
      break;
    case 'N':
      break;
    default:
      return false;
  }
  return true;
}

bool barcode::valid_value( char c ) {
  switch (c) {
    case '0':
      break;
    case '1':
      break;
    case '2':
      break;
    case '3':
      break;
    case '4':
      break;
    case '5':
      break;
    case '6':
      break;
    case '7':
      break;
    case '8':
      break;
    case '9':
      break;
    default:
      return false;
  }
  return true;
}



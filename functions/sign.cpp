//sign.cpp: this function returns the sign of some input value 
//Marco, 29 June 2011

#include "wet.h"
using namespace std;

double wet::sign(double value)
{
  if (value < 0.0) return -1.0;
  return 1.0;
}

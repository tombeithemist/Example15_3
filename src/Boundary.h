/*
 * Boundary.h
 *
 *  Created on: 22 Oct 2014
 *      Author: David
 */

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <string>
using namespace std;


class Boundary {
public:
  Boundary();
  virtual ~Boundary();
  string name_;
  string type_;
  int N_;
  int* indices_;
  double value_;

};

#endif /* BOUNDARY_H_ */

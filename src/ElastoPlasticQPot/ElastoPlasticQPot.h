/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_H
#define ELASTOPLASTICQPOT_H

// --------------------------------------- include libraries ---------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <cppmat/cppmat.h>

// ---------------------------------------- dummy operation ----------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------- version information --------------------------------------

#define ELASTOPLASTICQPOT_WORLD_VERSION 0
#define ELASTOPLASTICQPOT_MAJOR_VERSION 0
#define ELASTOPLASTICQPOT_MINOR_VERSION 1

#define ELASTOPLASTICQPOT_VERSION_AT_LEAST(x,y,z) \
  (ELASTOPLASTICQPOT_WORLD_VERSION>x || (ELASTOPLASTICQPOT_WORLD_VERSION>=x && \
  (ELASTOPLASTICQPOT_MAJOR_VERSION>y || (ELASTOPLASTICQPOT_MAJOR_VERSION>=y && \
                              ELASTOPLASTICQPOT_MINOR_VERSION>=z))))

#define ELASTOPLASTICQPOT_VERSION(x,y,z) \
  (ELASTOPLASTICQPOT_WORLD_VERSION==x && \
   ELASTOPLASTICQPOT_MAJOR_VERSION==y && \
   ELASTOPLASTICQPOT_MINOR_VERSION==z)

// ------------------------------------------ alias types ------------------------------------------

namespace ElastoPlasticQPot {

typedef cppmat::matrix <size_t> ArrS;
typedef cppmat::matrix <double> ArrD;
typedef cppmat::matrix2<size_t> MatS;
typedef cppmat::matrix2<double> MatD;
typedef cppmat::vector <size_t> ColS;
typedef cppmat::vector <double> ColD;
}

// ---------------------------------------- include headers ----------------------------------------

#include "Cartesian2d.h"

// ---------------------------------------- include scripts ----------------------------------------

#include "Cartesian2d.cpp"
#include "Cartesian2d_Elastic.cpp"
#include "Cartesian2d_Cusp.cpp"
#include "Cartesian2d_Smooth.cpp"
#include "Cartesian2d_Matrix.cpp"

// -------------------------------------------------------------------------------------------------

#endif

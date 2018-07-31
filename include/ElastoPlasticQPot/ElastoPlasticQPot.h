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
#define ELASTOPLASTICQPOT_MINOR_VERSION 4

#define ELASTOPLASTICQPOT_VERSION_AT_LEAST(x,y,z) \
  (ELASTOPLASTICQPOT_WORLD_VERSION>x || (ELASTOPLASTICQPOT_WORLD_VERSION>=x && \
  (ELASTOPLASTICQPOT_MAJOR_VERSION>y || (ELASTOPLASTICQPOT_MAJOR_VERSION>=y && \
                              ELASTOPLASTICQPOT_MINOR_VERSION>=z))))

#define ELASTOPLASTICQPOT_VERSION(x,y,z) \
  (ELASTOPLASTICQPOT_WORLD_VERSION==x && \
   ELASTOPLASTICQPOT_MAJOR_VERSION==y && \
   ELASTOPLASTICQPOT_MINOR_VERSION==z)

// ---------------------------------------- include headers ----------------------------------------

#include "Cartesian2d.h"

// ---------------------------------------- include scripts ----------------------------------------

#include "Cartesian2d.hpp"
#include "Cartesian2d_Elastic.hpp"
#include "Cartesian2d_Cusp.hpp"
#include "Cartesian2d_Smooth.hpp"
#include "Cartesian2d_Matrix.hpp"

// -------------------------------------------------------------------------------------------------

#endif

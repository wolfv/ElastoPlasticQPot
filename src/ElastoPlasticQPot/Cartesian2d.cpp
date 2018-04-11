/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_CPP
#define ELASTOPLASTICQPOT_CARTESIAN2D_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace ElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ mean strain ------------------------------------------

inline double epsm(const T2s &Eps)
{
  return Eps.trace()/2.;
}

// ---------------------------------- equivalent strain deviator -----------------------------------

inline double epsd(const T2s &Eps)
{
  T2d    I    = cm::identity2<double>();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;

  return std::sqrt(.5*Epsd.ddot(Epsd));
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

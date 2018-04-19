/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_ELASTIC_CPP
#define ELASTOPLASTICQPOT_CARTESIAN2D_ELASTIC_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace ElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

// ----------------------------------------- yield stress ------------------------------------------

inline double Elastic::epsy(size_t i) const
{
  UNUSED(i);

  return std::numeric_limits<double>::infinity();
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Elastic::find(const T2s &Eps) const
{
  UNUSED(Eps);

  return 0;
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Elastic::find(double epsd) const
{
  UNUSED(epsd);

  return 0;
}

// -------------------------------------------- stress ---------------------------------------------

inline T2s Elastic::stress(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2<double>();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;

  // return stress tensor
  return ( m_K * epsm ) * I + m_G * Epsd;
}

// -------------------------------------------- energy ---------------------------------------------

inline double Elastic::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2<double>();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // hydrostatic and deviatoric part of the energy
  double U = m_K * std::pow(epsm,2.);
  double V = m_G * std::pow(epsd,2.);

  // return total strain energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

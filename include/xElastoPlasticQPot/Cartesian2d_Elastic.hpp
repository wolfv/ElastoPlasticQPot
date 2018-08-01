/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef XELASTOPLASTICQPOT_CARTESIAN2D_ELASTIC_CPP
#define XELASTOPLASTICQPOT_CARTESIAN2D_ELASTIC_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace xElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

// ------------------------------------------ parameters -------------------------------------------

inline double Elastic::K() const
{
  return m_K;
}

// ------------------------------------------ parameters -------------------------------------------

inline double Elastic::G() const
{
  return m_G;
}

// ---------------------------------- equivalent deviator strain -----------------------------------

inline double Elastic::epsd(const T2s &Eps) const
{
  T2s  I    = eye();
  auto Epsd = Eps - trace(Eps)/2. * I;

  return std::sqrt(.5*ddot(Epsd,Epsd));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Elastic::epsp(const T2s &Eps) const
{
  UNUSED(Eps);

  return 0.0;
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Elastic::epsp(double epsd) const
{
  UNUSED(epsd);

  return 0.0;
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

inline T2s Elastic::Sig(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I    = eye();
  auto epsm = trace(Eps)/2.;
  auto Epsd = Eps - epsm * I;

  // return stress tensor
  return m_K * epsm * I + m_G * Epsd;
}

// -------------------------------------------- energy ---------------------------------------------

inline double Elastic::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I    = eye();
  auto epsm = trace(Eps)/2.;
  auto Epsd = Eps - epsm * I;
  auto epsd = std::sqrt(.5*ddot(Epsd,Epsd));

  // hydrostatic part of the energy
  auto U = m_K * std::pow(epsm,2.);
  // deviatoric part of the energy
  auto V = m_G * std::pow(epsd,2.);

  // return total energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

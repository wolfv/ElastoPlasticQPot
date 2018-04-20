/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_CUSP_CPP
#define ELASTOPLASTICQPOT_CARTESIAN2D_CUSP_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace ElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Cusp::Cusp(double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;

  // copy input - yield strains
  std::vector<double> vec = epsy;
  // sort input
  std::sort(vec.begin(), vec.end());

  // check to add item to have an initial elastic response
  if ( init_elastic )
    if ( vec[0] == -vec[1] )
      init_elastic = false;

  // copy input
  // - counters
  size_t N = vec.size();
  size_t i = 0;
  // - add yield strain to have an initial elastic response
  if ( init_elastic ) ++N;
  // - allocate
  m_epsy.resize(N);
  // - add yield strain to have an initial elastic response
  if ( init_elastic ) { m_epsy[i] = -vec[0]; ++i; }
  // - copy the rest
  for ( auto &j : vec ) { m_epsy[i] = j; ++i; }

  // check the number of yield strains
  if ( m_epsy.size() < 2 )
    throw std::runtime_error("Specify at least two yield strains 'epsy'");
}

// ---------------------------------- equivalent deviator strain -----------------------------------

inline double Cusp::epsd(const T2s &Eps) const
{
  T2s Epsd = Eps - Eps.trace()/2. * cppmat::cartesian2d::identity2<double>();

  return std::sqrt(.5*Epsd.ddot(Epsd));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Cusp::epsp(const T2s &Eps) const
{
  return epsp(epsd(Eps));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Cusp::epsp(double epsd) const
{
  size_t i = find(epsd);

  return ( m_epsy[i+1] + m_epsy[i] ) / 2.;
}

// ----------------------------------------- yield stress ------------------------------------------

inline double Cusp::epsy(size_t i) const
{
  return m_epsy[i];
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Cusp::find(const T2s &Eps) const
{
  return find(epsd(Eps));
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Cusp::find(double epsd) const
{
  // check extremes
  if ( epsd < m_epsy.front() or epsd >= m_epsy.back() )
    throw std::runtime_error("Insufficient 'epsy'");

  // set initial search bounds and index
  size_t n = m_epsy.size()-1;  // upper-bound
  size_t z = 1;                // lower-bound
  size_t l = 0;                // left-bound
  size_t r = n;                // right-bound
  size_t i = r / 2;            // average

  // loop until found
  while ( true )
  {
    // check if found, unroll once to speed-up
    if ( epsd >= m_epsy[i-1] and epsd < m_epsy[i  ] ) return i-1;
    if ( epsd >= m_epsy[i  ] and epsd < m_epsy[i+1] ) return i;
    if ( epsd >= m_epsy[i+1] and epsd < m_epsy[i+2] ) return i+1;

    // correct the left- and right-bound
    if ( epsd >= m_epsy[i] ) l = i;
    else                     r = i;

    // set new search index
    i = ( r + l ) / 2;
    i = std::max(i,z);
  }
}

// -------------------------------------------- stress ---------------------------------------------

inline T2s Cusp::stress(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2<double>();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // no deviatoric strain -> only hydrostatic stress
  if ( epsd <= 0. ) return (m_K*epsm) * I;

  // read current yield strains
  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;

  // return stress tensor
  return (m_K*epsm)*I + ( m_G * (1.-eps_min/epsd) ) * Epsd;
}

// -------------------------------------------- energy ---------------------------------------------

inline double Cusp::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2<double>();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // hydrostatic part of the energy
  double U = m_K * std::pow(epsm,2.);

  // read current yield strain
  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // deviatoric part of the energy
  double V = m_G * ( std::pow(epsd-eps_min,2.) - std::pow(deps_y,2.) );

  // return total energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

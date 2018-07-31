/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef XELASTOPLASTICQPOT_CARTESIAN2D_SMOOTH_CPP
#define XELASTOPLASTICQPOT_CARTESIAN2D_SMOOTH_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace xElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Smooth::Smooth(double K, double G, const std::vector<double> &epsy, bool init_elastic)
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

// ------------------------------------------ parameters -------------------------------------------

inline double Smooth::K() const
{
  return m_K;
}

// ------------------------------------------ parameters -------------------------------------------

inline double Smooth::G() const
{
  return m_G;
}

// ---------------------------------- equivalent deviator strain -----------------------------------

inline double Smooth::epsd(const T2s &Eps) const
{
  auto Epsd = Eps - trace(Eps)/ND * xt::eye(ndim);

  return std::sqrt(.5*ddot(Epsd,Epsd));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Smooth::epsp(const T2s &Eps) const
{
  return epsp(epsd(Eps));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Smooth::epsp(double epsd) const
{
  size_t i = find(epsd);

  return ( m_epsy[i+1] + m_epsy[i] ) / 2.;
}

// ----------------------------------------- yield stress ------------------------------------------

inline double Smooth::epsy(size_t i) const
{
  return m_epsy[i];
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Smooth::find(const T2s &Eps) const
{
  return find(epsd(Eps));
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Smooth::find(double epsd) const
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

inline T2s Smooth::Sig(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I    = xt::eye(ndim);
  auto epsm = trace(Eps)/ND;
  auto Epsd = Eps - epsm * I;
  auto epsd = std::sqrt(.5*ddot(Epsd,Epsd));

  // no deviatoric strain -> only hydrostatic stress
  if ( epsd <= 0. ) return m_K * epsm * I;

  // read current yield strains
  auto i       = find(epsd);
  auto eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  auto deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // return stress tensor
  return m_K*epsm*I + (m_G/epsd)*(deps_y/M_PI)*sin(M_PI/deps_y*(epsd-eps_min))*Epsd;
}

// -------------------------------------------- energy ---------------------------------------------

inline double Smooth::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I    = xt::eye(ndim);
  auto epsm = trace(Eps)/ND;
  auto Epsd = Eps - epsm * I;
  auto epsd = std::sqrt(.5*ddot(Epsd,Epsd));

  // hydrostatic part of the energy
  auto U = m_K * std::pow(epsm,2.);

  // read current yield strain
  auto i       = find(epsd);
  auto eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  auto deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // deviatoric part of the energy
  auto V = -2.*m_G * std::pow(deps_y/M_PI,2.) * ( 1. + cos( M_PI/deps_y * (epsd-eps_min) ) );

  // return total energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

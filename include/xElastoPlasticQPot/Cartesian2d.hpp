/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef XELASTOPLASTICQPOT_CARTESIAN2D_HPP
#define XELASTOPLASTICQPOT_CARTESIAN2D_HPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace xElastoPlasticQPot {
namespace Cartesian2d {

// ======================================== TENSOR PRODUCTS ========================================

// --------------------------------------------- trace ---------------------------------------------

template<class T>
inline double trace(const T &A)
{
  return A(0,0) + A(1,1);
}

// -------------------------------------- double dot product ---------------------------------------

template<class T>
inline double ddot(const T &A, const T &B)
{
  return A(0,0) * B(0,0) + 2.0 * A(0,1) * B(0,1) + A(1,1) * B(1,1);
}

// ------------------------------------------ unit tensor ------------------------------------------

inline T2s eye()
{
  T2s I;

  I[0] = I[3] = 1.;
  I[1] = I[2] = 0.;

  return I;
}

// ===================================== TENSOR DECOMPOSITION ======================================

// -------------------------------------- hydrostatic strain ---------------------------------------

inline double epsm(const T2s &Eps)
{
  return trace(Eps)/2.;
}

// --------------------------------- equivalent deviatoric strain ----------------------------------

inline double epsd(const T2s &Eps)
{
  T2s I    = eye();
  T2s Epsd = Eps - trace(Eps)/2. * I;

  return std::sqrt(.5*ddot(Epsd,Epsd));
}

// ---------------------------------------- strain deviator ----------------------------------------

inline T2s Epsd(const T2s &Eps)
{
  T2s I = eye();

  return Eps - trace(Eps)/2. * I;
}

// -------------------------------------- hydrostatic stress ---------------------------------------

inline double sigm(const T2s &Sig)
{
  return trace(Sig)/2.;
}

// --------------------------------- equivalent deviatoric stress ----------------------------------

inline double sigd(const T2s &Sig)
{
  T2s I    = eye();
  T2s Sigd = Sig - trace(Sig)/2. * I;

  return std::sqrt(2.*ddot(Sigd,Sigd));
}

// ---------------------------------------- stress deviator ----------------------------------------

inline T2s Sigd(const T2s &Sig)
{
  T2s I = eye();

  return Sig - trace(Sig)/2. * I;
}

// ================================= TENSOR DECOMPOSITION - MATRIX =================================

// -------------------------------------- hydrostatic strain ---------------------------------------

inline void epsm(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsm)
{
  // check input
  assert( a_Eps.shape()[2] == 2 );
  assert( a_Eps.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Eps.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Eps.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
        // - equivalent value
        a_epsm(e,k) = trace(Eps)/2.;
      }
    }
  }
}

// -------------------------------- hydrostatic strain - interface ---------------------------------

inline xt::xtensor<double,2> epsm(const xt::xtensor<double,4> &a_Eps)
{
  xt::xtensor<double,2> out = xt::empty<double>({a_Eps.shape()[0], a_Eps.shape()[1]});

  epsm(a_Eps, out);

  return out;
}

// --------------------------------- equivalent deviatoric strain ----------------------------------

inline void epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsd)
{
  // check input
  assert( a_Eps.shape()[2] == 2 );
  assert( a_Eps.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2s I = eye();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Eps.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Eps.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
        // - strain deviator
        auto Epsd = Eps - trace(Eps)/2. * I;
        // - equivalent value
        a_epsd(e,k) = std::sqrt(.5*ddot(Epsd,Epsd));
      }
    }
  }
}

// --------------------------- equivalent deviatoric strain - interface ----------------------------

inline xt::xtensor<double,2> epsd(const xt::xtensor<double,4> &a_Eps)
{
  xt::xtensor<double,2> out = xt::empty<double>({a_Eps.shape()[0], a_Eps.shape()[1]});

  epsd(a_Eps, out);

  return out;
}

// ---------------------------------------- strain deviator ----------------------------------------

inline void Epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Epsd)
{
  // check input
  assert( a_Eps .shape()[2] == 2 );
  assert( a_Eps .shape()[3] == 2 );
  assert( a_Epsd.shape()[2] == 2 );
  assert( a_Epsd.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2s I = eye();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Eps.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Eps.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Eps  = xt::view(a_Eps , e, k, xt::all(), xt::all());
        auto Epsd = xt::view(a_Epsd, e, k, xt::all(), xt::all());
        // - strain deviator
        Epsd = Eps - trace(Eps)/2. * I;
      }
    }
  }
}

// ---------------------------------- strain deviator - interface ----------------------------------

inline xt::xtensor<double,4> Epsd(const xt::xtensor<double,4> &a_Eps)
{
  xt::xtensor<double,4> out = xt::empty<double>(a_Eps.shape());

  Epsd(a_Eps, out);

  return out;
}

// -------------------------------------- hydrostatic strain ---------------------------------------

inline void sigm(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigm)
{
  // check input
  assert( a_Sig.shape()[2] == 2 );
  assert( a_Sig.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Sig.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Sig.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Sig = xt::view(a_Sig, e, k, xt::all(), xt::all());
        // - equivalent value
        a_sigm(e,k) = trace(Sig)/2.;
      }
    }
  }
}

// -------------------------------- hydrostatic strain - interface ---------------------------------

inline xt::xtensor<double,2> sigm(const xt::xtensor<double,4> &a_Sig)
{
  xt::xtensor<double,2> out = xt::empty<double>({a_Sig.shape()[0], a_Sig.shape()[1]});

  sigm(a_Sig, out);

  return out;
}

// --------------------------------- equivalent deviatoric strain ----------------------------------

inline void sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigd)
{
  // check input
  assert( a_Sig.shape()[2] == 2 );
  assert( a_Sig.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2s I = eye();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Sig.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Sig.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Sig = xt::view(a_Sig, e, k, xt::all(), xt::all());
        // - strain deviator
        auto Sigd = Sig - trace(Sig)/2. * I;
        // - equivalent value
        a_sigd(e,k) = std::sqrt(2.*ddot(Sigd,Sigd));
      }
    }
  }
}

// --------------------------- equivalent deviatoric strain - interface ----------------------------

inline xt::xtensor<double,2> sigd(const xt::xtensor<double,4> &a_Sig)
{
  xt::xtensor<double,2> out = xt::empty<double>({a_Sig.shape()[0], a_Sig.shape()[1]});

  sigd(a_Sig, out);

  return out;
}

// ---------------------------------------- strain deviator ----------------------------------------

inline void Sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,4> &a_Sigd)
{
  // check input
  assert( a_Sig .shape()[2] == 2 );
  assert( a_Sig .shape()[3] == 2 );
  assert( a_Sigd.shape()[2] == 2 );
  assert( a_Sigd.shape()[3] == 2 );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2s I = eye();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < a_Sig.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < a_Sig.shape()[1] ; ++k )
      {
        // - strain tensor
        auto Sig  = xt::view(a_Sig , e, k, xt::all(), xt::all());
        auto Sigd = xt::view(a_Sigd, e, k, xt::all(), xt::all());
        // - strain deviator
        Sigd = Sig - trace(Sig)/2. * I;
      }
    }
  }
}

// ---------------------------------- strain deviator - interface ----------------------------------

inline xt::xtensor<double,4> Sigd(const xt::xtensor<double,4> &a_Sig)
{
  xt::xtensor<double,4> out = xt::empty<double>(a_Sig.shape());

  Sigd(a_Sig, out);

  return out;
}

// ============================================ MAXIMUM ============================================

// -------------------------------------- hydrostatic strain ---------------------------------------

inline double epsm_max(const xt::xtensor<double,4> &a_Eps)
{
  // check input
  assert( a_Eps.shape()[2] == 2 );
  assert( a_Eps.shape()[3] == 2 );

  // allocate maximum
  double out;

  // compute one point
  {
    // - strain tensor
    auto Eps = xt::view(a_Eps, 0, 0, xt::all(), xt::all());
    // - equivalent value
    out = trace(Eps)/2.;
  }

  // loop over all points
  for ( size_t e = 0 ; e < a_Eps.shape()[0] ; ++e )
  {
    for ( size_t k = 0 ; k < a_Eps.shape()[1] ; ++k )
    {
      // - strain tensor
      auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
      // - equivalent value
      out = std::max(out, trace(Eps)/2.);
    }
  }

  return out;
}

// -------------------------------------- hydrostatic stress ---------------------------------------

inline double sigm_max(const xt::xtensor<double,4> &a_Sig)
{
  // check input
  assert( a_Sig.shape()[2] == 2 );
  assert( a_Sig.shape()[3] == 2 );

  // allocate maximum
  double out;

  // compute one point
  {
    // - stress tensor
    auto Sig = xt::view(a_Sig, 0, 0, xt::all(), xt::all());
    // - equivalent value
    out = trace(Sig)/2.;
  }

  // loop over all points
  for ( size_t e = 0 ; e < a_Sig.shape()[0] ; ++e )
  {
    for ( size_t k = 0 ; k < a_Sig.shape()[1] ; ++k )
    {
      // - stress tensor
      auto Sig = xt::view(a_Sig, e, k, xt::all(), xt::all());
      // - equivalent value
      out = std::max(out, trace(Sig)/2.);
    }
  }

  return out;
}

// --------------------------------- equivalent deviatoric strain ----------------------------------

inline double epsd_max(const xt::xtensor<double,4> &a_Eps)
{
  // check input
  assert( a_Eps.shape()[2] == 2 );
  assert( a_Eps.shape()[3] == 2 );

  // identity tensor
  T2s I = eye();

  // allocate maximum
  double out;

  // compute one point
  {
    // - strain tensor
    auto Eps = xt::view(a_Eps, 0, 0, xt::all(), xt::all());
    // - strain deviator
    auto Epsd = Eps - trace(Eps)/2. * I;
    // - equivalent value
    out = std::sqrt(.5*ddot(Epsd,Epsd));
  }

  // loop over all points
  for ( size_t e = 0 ; e < a_Eps.shape()[0] ; ++e )
  {
    for ( size_t k = 0 ; k < a_Eps.shape()[1] ; ++k )
    {
      // - strain tensor
      auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
      // - strain deviator
      auto Epsd = Eps - trace(Eps)/2. * I;
      // - equivalent value
      out = std::max(out, std::sqrt(.5*ddot(Epsd,Epsd)));
    }
  }

  return out;
}

// --------------------------------- equivalent deviatoric stress ----------------------------------

inline double sigd_max(const xt::xtensor<double,4> &a_Sig)
{
  // check input
  assert( a_Sig.shape()[2] == 2 );
  assert( a_Sig.shape()[3] == 2 );

  // identity tensor
  T2s I = eye();

  // allocate maximum
  double out;

  // compute one point
  {
    // - stress tensor
    auto Sig = xt::view(a_Sig, 0, 0, xt::all(), xt::all());
    // - stress deviator
    auto Sigd = Sig - trace(Sig)/2. * I;
    // - equivalent value
    out = std::sqrt(.5*ddot(Sigd,Sigd));
  }

  // loop over all points
  for ( size_t e = 0 ; e < a_Sig.shape()[0] ; ++e )
  {
    for ( size_t k = 0 ; k < a_Sig.shape()[1] ; ++k )
    {
      // - stress tensor
      auto Sig = xt::view(a_Sig, e, k, xt::all(), xt::all());
      // - stress deviator
      auto Sigd = Sig - trace(Sig)/2. * I;
      // - equivalent value
      out = std::max(out, std::sqrt(2.*ddot(Sigd,Sigd)));
    }
  }

  return out;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

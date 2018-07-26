/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_HPP
#define ELASTOPLASTICQPOT_CARTESIAN2D_HPP

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
  T2s Epsd = Eps - Eps.trace()/2. * T2d::I();

  return std::sqrt(.5*Epsd.ddot(Epsd));
}

// ---------------------------------------- strain deviator ----------------------------------------

inline T2s Epsd(const T2s &Eps)
{
  return Eps - Eps.trace()/2. * T2d::I();
}

// ------------------------------------------ mean stress ------------------------------------------

inline double sigm(const T2s &Sig)
{
  return Sig.trace()/2.;
}

// ---------------------------------- equivalent stress deviator -----------------------------------

inline double sigd(const T2s &Sig)
{
  T2s Sigd = Sig - Sig.trace()/2. * T2d::I();

  return std::sqrt(2.*Sigd.ddot(Sigd));
}

// ---------------------------------------- stress deviator ----------------------------------------

inline T2s Sigd(const T2s &Sig)
{
  return Sig - Sig.trace()/2. * T2d::I();
}

// ------------------------------------- mean strain - matrix --------------------------------------

inline ArrD epsm(const ArrD &a_Eps)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Eps.size() / ncomp;

  // check input
  assert( a_Eps.rank()    >= 2     );
  assert( a_Eps.shape(-1) == ncomp );

  // allocate output: matrix of scalars (shape of the input matrix, without tensor-components)
  ArrD a_epsm(cppmat::del(a_Eps.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Eps = T2s::Copy(a_Eps.index(i*ncomp));
      // compute/store the mean strain
      a_epsm[i] = Eps.trace()/2.;
    }
  }

  return a_epsm;
}

// ------------------------------ equivalent strain deviator - matrix ------------------------------

inline ArrD epsd(const ArrD &a_Eps)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Eps.size() / ncomp;

  // check input
  assert( a_Eps.rank()    >= 2     );
  assert( a_Eps.shape(-1) == ncomp );

  // allocate output: matrix of scalars (shape of the input matrix, without tensor-components)
  ArrD a_epsd(cppmat::del(a_Eps.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Eps = T2s::Copy(a_Eps.index(i*ncomp));
      // compute the strain deviator
      T2s Epsd = Eps - Eps.trace()/2. * T2d::I();
      // compute/store the equivalent strain deviator
      a_epsd[i] = std::sqrt(.5*Epsd.ddot(Epsd));
    }
  }

  return a_epsd;
}

// ----------------------------------- strain deviator - matrix ------------------------------------

inline ArrD Epsd(const ArrD &a_Eps)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Eps.size() / ncomp;

  // check input
  assert( a_Eps.rank()    >= 2     );
  assert( a_Eps.shape(-1) == ncomp );

  // allocate output: matrix of tensors
  ArrD a_Epsd(a_Eps.shape());

  // iterators
  auto i_Epsd = a_Epsd.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Eps = T2s::Copy(a_Eps.index(i*ncomp));
      // compute the strain deviator
      T2s Epsd = Eps - Eps.trace()/2. * T2d::I();
      // store the strain deviator
      std::copy(Epsd.begin(), Epsd.end(), i_Epsd+i*ncomp);
    }
  }

  return a_Epsd;
}

// ------------------------------------- mean stress - matrix --------------------------------------

inline ArrD sigm(const ArrD &a_Sig)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Sig.size() / ncomp;

  // check input
  assert( a_Sig.rank()    >= 2     );
  assert( a_Sig.shape(-1) == ncomp );

  // allocate output: matrix of scalars (shape of the input matrix, without tensor-components)
  ArrD a_sigm(cppmat::del(a_Sig.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Sig = T2s::Copy(a_Sig.index(i*ncomp));
      // compute/store the mean stress
      a_sigm[i] = Sig.trace()/2.;
    }
  }

  return a_sigm;
}

// ------------------------------ equivalent stress deviator - matrix ------------------------------

inline ArrD sigd(const ArrD &a_Sig)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Sig.size() / ncomp;

  // check input
  assert( a_Sig.rank()    >= 2     );
  assert( a_Sig.shape(-1) == ncomp );

  // allocate output: matrix of scalars (shape of the input matrix, without tensor-components)
  ArrD a_sigd(cppmat::del(a_Sig.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Sig = T2s::Copy(a_Sig.index(i*ncomp));
      // compute the stress deviator
      T2s Sigd = Sig - Sig.trace()/2. * T2d::I();
      // compute/store the equivalent stress deviator
      a_sigd[i] = std::sqrt(2.*Sigd.ddot(Sigd));
    }
  }

  return a_sigd;
}

// ----------------------------------- stress deviator - matrix ------------------------------------

inline ArrD Sigd(const ArrD &a_Sig)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Sig.size() / ncomp;

  // check input
  assert( a_Sig.rank()    >= 2     );
  assert( a_Sig.shape(-1) == ncomp );

  // allocate output: matrix of tensors
  ArrD a_Sigd(a_Sig.shape());

  // iterators
  auto i_Sigd = a_Sigd.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < N ; ++i )
    {
      // map from matrix of strains
      T2s Sig = T2s::Copy(a_Sig.index(i*ncomp));
      // compute the stress deviator
      T2s Sigd = Sig - Sig.trace()/2. * T2d::I();
      // store the stress deviator
      std::copy(Sigd.begin(), Sigd.end(), i_Sigd+i*ncomp);
    }
  }

  return a_Sigd;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

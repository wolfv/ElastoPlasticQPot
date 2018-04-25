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
  T2s Epsd = Eps - Eps.trace()/2. * T2d::I();

  return std::sqrt(.5*Epsd.ddot(Epsd));
}

// ---------------------------------------- strain deviator ----------------------------------------

inline T2s Epsd(const T2s &Eps)
{
  return Eps - Eps.trace()/2. * T2d::I();
}

// ----------------------------------- mean & equivalent stress ------------------------------------

inline double sigm(const T2s &Sig) { return epsm(Sig); }
inline double sigd(const T2s &Sig) { return epsd(Sig); }
inline T2s    Sigd(const T2s &Sig) { return Epsd(Sig); }

// ------------------------------------- mean strain - matrix --------------------------------------

inline ArrD epsm(const ArrD &a_Eps)
{
  // number of tensor-components
  static const size_t ncomp = 3;
  // number of entries
  size_t N = a_Eps.size() / ncomp;

  // check input
  assert( a_Eps.ndim()    >= 2     );
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
      vT2s Eps = vT2s::Map(&a_Eps[i*ncomp]);
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
  assert( a_Eps.ndim()    >= 2     );
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
      vT2s Eps  = vT2s::Map(&a_Eps[i*ncomp]);
      // compute the strain deviator
      T2s  Epsd = Eps - Eps.trace()/2. * T2d::I();
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
  assert( a_Eps.ndim()    >= 2     );
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
      vT2s Eps  = vT2s::Map(&a_Eps[i*ncomp]);
      // compute the strain deviator
      T2s  Epsd = Eps - Eps.trace()/2. * T2d::I();
      // store the strain deviator
      std::copy(Epsd.begin(), Epsd.end(), i_Epsd+i*ncomp);
    }
  }

  return a_Epsd;
}

// ----------------------------------- mean & equivalent stress ------------------------------------

inline ArrD sigm(const ArrD &a_Sig) { return epsm(a_Sig); }
inline ArrD sigd(const ArrD &a_Sig) { return epsd(a_Sig); }
inline ArrD Sigd(const ArrD &a_Sig) { return epsd(a_Sig); }

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

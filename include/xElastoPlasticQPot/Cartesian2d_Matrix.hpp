/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef XELASTOPLASTICQPOT_CARTESIAN2D_MATRIX_CPP
#define XELASTOPLASTICQPOT_CARTESIAN2D_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace xElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Matrix::Matrix(const std::vector<size_t> &shape)
{
  // resize type and index look-up
  m_type .resize(shape);
  m_index.resize(shape);

  // initialize everything to unassigned
  m_type = Type::Unset * xt::ones<size_t>(shape);
}

// --------------------------------------------- shape ---------------------------------------------

inline std::vector<size_t> Matrix::shape() const
{
  return std::vector<size_t>(m_type.shape().begin(), m_type.shape().end());
}

// --------------------------------------------- shape ---------------------------------------------

inline size_t Matrix::shape(size_t i) const
{
  return m_type.shape()[i];
}

// --------------------------------------------- type ----------------------------------------------

inline xt::xtensor<size_t,2> Matrix::type() const
{
  return m_type;
}

// ------------------------------------------ parameters -------------------------------------------

inline xt::xtensor<double,2> Matrix::K() const
{
  xt::xtensor<double,2> a_K = xt::empty<double>(m_type.shape());

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      switch ( m_type[i] )
      {
        case Type::Elastic: a_K[i] = m_Elastic[m_index[i]].K(); break;
        case Type::Cusp   : a_K[i] = m_Cusp   [m_index[i]].K(); break;
        case Type::Smooth : a_K[i] = m_Smooth [m_index[i]].K(); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return a_K;
}

// ------------------------------------------ parameters -------------------------------------------

inline xt::xtensor<double,2> Matrix::G() const
{
  xt::xtensor<double,2> a_G = xt::empty<double>(m_type.shape());

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      switch ( m_type[i] )
      {
        case Type::Elastic: a_G[i] = m_Elastic[m_index[i]].G(); break;
        case Type::Cusp   : a_G[i] = m_Cusp   [m_index[i]].G(); break;
        case Type::Smooth : a_G[i] = m_Smooth [m_index[i]].G(); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return a_G;
}

// ---------------------------------------- type is plastic ----------------------------------------

inline xt::xtensor<size_t,2> Matrix::isPlastic() const
{
  xt::xtensor<size_t,2> out = xt::zeros<size_t>(m_type.shape());

  for ( size_t i = 0 ; i < m_type.size() ; ++i )
    if ( m_type[i] != Type::Unset and m_type[i] != Type::Elastic )
      out[i] = 1;

  return out;
}

// --------------------------- check that a type has been set everywhere ---------------------------

inline void Matrix::check() const
{
  for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
        if ( m_type(e,k) == Type::Unset )
          throw std::runtime_error("No type set for: "+std::to_string(e)+", "+std::to_string(k));
}

// ---------------------------------- set Elastic material points ----------------------------------

inline void Matrix::setElastic(const xt::xtensor<size_t,2> &I, double K, double G)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
  #endif

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Elastic;
      m_index[i] = m_Elastic.size();
    }
  }
  // store material definition
  m_Elastic.push_back(Elastic(K, G));
}

// ----------------------------------- set Cusp material points ------------------------------------

inline void Matrix::setCusp(const xt::xtensor<size_t,2> &I,
  double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
  #endif

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Cusp;
      m_index[i] = m_Cusp.size();
    }
  }
  // store material definition
  m_Cusp.push_back(Cusp(K, G, epsy, init_elastic));
}

// ---------------------------------- set Smooth material points -----------------------------------

inline void Matrix::setSmooth(const xt::xtensor<size_t,2> &I,
  double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
  #endif

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Smooth;
      m_index[i] = m_Smooth.size();
    }
  }
  // store material definition
  m_Smooth.push_back(Smooth(K, G, epsy, init_elastic));
}

// ---------------------------------- set Elastic material points ----------------------------------

inline void Matrix::setElastic(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
  const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    assert( I.shape() == idx   .shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
    // - extent
    assert( xt::amax(idx)[0] == K.size()-1 );
    // - consistency
    assert( K.size() == G.size() );
  #endif

  // start index
  size_t N = m_Elastic.size();

  // store material definition
  for ( size_t i = 0 ; i < K.size() ; ++i )
    m_Elastic.push_back(Elastic(K(i), G(i)));

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Elastic;
      m_index[i] = idx[i] + N;
    }
  }
}

// ----------------------------------- set Cusp material points ------------------------------------

inline void Matrix::setCusp(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
  const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G,
  const xt::xtensor<double,2> &epsy, bool init_elastic)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    assert( I.shape() == idx   .shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
    // - extent
    assert( xt::amax(idx)[0] == K.size()-1 );
    // - consistency
    assert( K.size() == G.size()        );
    assert( K.size() == epsy.shape()[0] );
  #endif

  // start index
  size_t N = m_Cusp.size();

  // store material definition
  for ( size_t i = 0 ; i < K.size() ; ++i ) {
    std::vector<double> y(epsy.begin()+i*epsy.shape()[1], epsy.begin()+(i+1)*epsy.shape()[1]);
    m_Cusp.push_back(Cusp(K(i), G(i), y, init_elastic));
  }

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Cusp;
      m_index[i] = idx[i] + N;
    }
  }
}

// ----------------------------------- set Smooth material points ------------------------------------

inline void Matrix::setSmooth(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
  const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G,
  const xt::xtensor<double,2> &epsy, bool init_elastic)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_type.shape() );
    assert( I.shape() == idx   .shape() );
    // - empty material definitions
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
      if ( I[i] )
        assert( m_type[i] == Type::Unset );
    // - extent
    assert( xt::amax(idx)[0] == K.size()-1 );
    // - consistency
    assert( K.size() == G.size()        );
    assert( K.size() == epsy.shape()[0] );
  #endif

  // start index
  size_t N = m_Smooth.size();

  // store material definition
  for ( size_t i = 0 ; i < K.size() ; ++i ) {
    std::vector<double> y(epsy.begin()+i*epsy.shape()[1], epsy.begin()+(i+1)*epsy.shape()[1]);
    m_Smooth.push_back(Smooth(K(i), G(i), y, init_elastic));
  }

  // set type and position in material vector
  for ( size_t i = 0 ; i < m_type.size() ; ++i ) {
    if ( I[i] ) {
      m_type [i] = Type::Smooth;
      m_index[i] = idx[i] + N;
    }
  }
}

// -------------------------------- compute stress for all entries ---------------------------------

inline void Matrix::Sig(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Sig) const
{
  // check input
  assert( a_Eps.shape()[0] == m_type.shape()[0] );
  assert( a_Eps.shape()[1] == m_type.shape()[1] );
  assert( a_Eps.shape()[2] == 2                 );
  assert( a_Eps.shape()[3] == 2                 );
  assert( a_Eps.shape()    == a_Sig.shape()     );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
      {
        // - alias
        auto Eps = xt::view(a_Eps, e, k);
        auto Sig = xt::view(a_Sig, e, k);
        // - compute
        switch ( m_type(e,k) )
        {
          case Type::Elastic: Sig = m_Elastic[m_index(e,k)].Sig(Eps); break;
          case Type::Cusp   : Sig = m_Cusp   [m_index(e,k)].Sig(Eps); break;
          case Type::Smooth : Sig = m_Smooth [m_index(e,k)].Sig(Eps); break;
          default: std::runtime_error("Unknown material");
        }
      }
    }
  }
}


// -------------------------------- compute energy for all entries ---------------------------------

inline void Matrix::energy(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_energy) const
{
  // check input
  assert( a_Eps.shape()[0] == m_type.shape()[0] );
  assert( a_Eps.shape()[1] == m_type.shape()[1] );
  assert( a_Eps.shape()[2] == 2                 );
  assert( a_Eps.shape()[3] == 2                 );
  assert( a_energy.shape() == m_type.shape()    );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
      {
        // - alias
        auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
        // - compute
        switch ( m_type(e,k) )
        {
          case Type::Elastic: a_energy(e,k) = m_Elastic[m_index(e,k)].energy(Eps); break;
          case Type::Cusp   : a_energy(e,k) = m_Cusp   [m_index(e,k)].energy(Eps); break;
          case Type::Smooth : a_energy(e,k) = m_Smooth [m_index(e,k)].energy(Eps); break;
          default: std::runtime_error("Unknown material");
        }
      }
    }
  }
}

// ------------------------- find the current yield strain for all entries -------------------------

inline void Matrix::find(const xt::xtensor<double,4> &a_Eps, xt::xtensor<size_t,2> &a_idx) const
{
  // check input
  assert( a_Eps.shape()[0] == m_type.shape()[0] );
  assert( a_Eps.shape()[1] == m_type.shape()[1] );
  assert( a_Eps.shape()[2] == 2                 );
  assert( a_Eps.shape()[3] == 2                 );
  assert( a_idx.shape()    == m_type.shape()    );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
      {
        // - alias
        auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
        // - compute
        switch ( m_type(e,k) )
        {
          case Type::Elastic: a_idx(e,k) = m_Elastic[m_index(e,k)].find(Eps); break;
          case Type::Cusp   : a_idx(e,k) = m_Cusp   [m_index(e,k)].find(Eps); break;
          case Type::Smooth : a_idx(e,k) = m_Smooth [m_index(e,k)].find(Eps); break;
          default: std::runtime_error("Unknown material");
        }
      }
    }
  }
}

// ----------------------------- get the yield strain for all entries ------------------------------

inline void Matrix::epsy(const xt::xtensor<size_t,2> &a_idx, xt::xtensor<double,2> &a_epsy) const
{
  // check input
  assert( a_idx.shape()  == m_type.shape() );
  assert( a_epsy.shape() == m_type.shape() );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
      {
        switch ( m_type(e,k) )
        {
          case Type::Elastic: a_epsy(e,k) = m_Elastic[m_index(e,k)].epsy(a_idx(e,k)); break;
          case Type::Cusp   : a_epsy(e,k) = m_Cusp   [m_index(e,k)].epsy(a_idx(e,k)); break;
          case Type::Smooth : a_epsy(e,k) = m_Smooth [m_index(e,k)].epsy(a_idx(e,k)); break;
          default: std::runtime_error("Unknown material");
        }
      }
    }
  }
}

// ----------------------- get the equivalent plastic strain for all entries -----------------------

inline void Matrix::epsp(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsp) const
{
  // check input
  assert( a_Eps.shape()[0] == m_type.shape()[0] );
  assert( a_Eps.shape()[1] == m_type.shape()[1] );
  assert( a_Eps.shape()[2] == 2                 );
  assert( a_Eps.shape()[3] == 2                 );
  assert( a_epsp.shape()   == m_type.shape()    );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_type.shape()[0] ; ++e )
    {
      for ( size_t k = 0 ; k < m_type.shape()[1] ; ++k )
      {
        // - alias
        auto Eps = xt::view(a_Eps, e, k, xt::all(), xt::all());
        // - compute
        switch ( m_type(e,k) )
        {
          case Type::Elastic: a_epsp(e,k) = m_Elastic[m_index(e,k)].epsp(Eps); break;
          case Type::Cusp   : a_epsp(e,k) = m_Cusp   [m_index(e,k)].epsp(Eps); break;
          case Type::Smooth : a_epsp(e,k) = m_Smooth [m_index(e,k)].epsp(Eps); break;
          default: std::runtime_error("Unknown material");
        }
      }
    }
  }
}

// -------------------------------- compute stress for all entries ---------------------------------

inline xt::xtensor<double,4> Matrix::Sig(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,4> a_Sig = xt::empty<double>(a_Eps.shape());

  this->Sig(a_Eps, a_Sig);

  return a_Sig;
}


// -------------------------------- compute energy for all entries ---------------------------------

inline xt::xtensor<double,2> Matrix::energy(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,2> a_energy = xt::empty<double>(m_type.shape());

  this->energy(a_Eps, a_energy);

  return a_energy;
}

// ------------------------- find the current yield strain for all entries -------------------------

inline xt::xtensor<size_t,2> Matrix::find(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<size_t,2> a_idx = xt::empty<size_t>(m_type.shape());

  this->find(a_Eps, a_idx);

  return a_idx;
}

// ----------------------------- get the yield strain for all entries ------------------------------

inline xt::xtensor<double,2> Matrix::epsy(const xt::xtensor<size_t,2> &a_idx) const
{
  xt::xtensor<double,2> a_epsy = xt::empty<double>(m_type.shape());

  this->epsy(a_idx, a_epsy);

  return a_epsy;
}

// ----------------------- get the equivalent plastic strain for all entries -----------------------

inline xt::xtensor<double,2> Matrix::epsp(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,2> a_epsp = xt::empty<double>(m_type.shape());

  this->epsp(a_Eps, a_epsp);

  return a_epsp;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

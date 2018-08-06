/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_MATRIX_CPP
#define ELASTOPLASTICQPOT_CARTESIAN2D_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// =================================================================================================

namespace ElastoPlasticQPot {
namespace Cartesian2d {

// ------------------------------------------ constructor ------------------------------------------

inline Matrix::Matrix(const std::vector<size_t> &shape)
{
  // resize type and index look-up
  m_type .resize(shape);
  m_index.resize(shape);

  // initialize everything to unassigned
  m_type.setConstant(Type::Unset);
}

// --------------------------------------------- shape ---------------------------------------------

inline std::vector<size_t> Matrix::shape() const
{
  return m_type.shape();
}

// --------------------------------------------- shape ---------------------------------------------

inline size_t Matrix::shape(size_t i) const
{
  return m_type.shape(i);
}

// --------------------------------------------- type ----------------------------------------------

inline ArrS Matrix::type() const
{
  return m_type;
}

// ------------------------------------------ parameters -------------------------------------------

inline ArrD Matrix::K() const
{
  ArrD a_K(m_type.shape());

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

inline ArrD Matrix::G() const
{
  ArrD a_G(m_type.shape());

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

inline ArrS Matrix::isPlastic() const
{
  ArrS out = ArrS::Zero(m_type.shape());

  for ( size_t i = 0 ; i < m_type.size() ; ++i )
    if ( m_type[i] != Type::Unset and m_type[i] != Type::Elastic )
      out[i] = 1;

  return out;
}

// --------------------------- check that a type has been set everywhere ---------------------------

inline void Matrix::check() const
{
  for ( size_t i = 0 ; i < m_type.size() ; ++i )
    if ( m_type[i] == Type::Unset )
      throw std::runtime_error("No type set for: " + cppmat::to_string(m_type.decompress(i)));
}

// ---------------------------------- set Elastic material points ----------------------------------

inline void Matrix::setElastic(const ArrS &I, double K, double G)
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

inline void Matrix::setCusp(const ArrS &I,
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

inline void Matrix::setSmooth(const ArrS &I,
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

inline void Matrix::setElastic(const ArrS &I, const ArrS &idx, const ColD &K, const ColD &G)
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
    assert( idx.max() == K.size()-1 );
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

inline void Matrix::setCusp(const ArrS &I, const ArrS &idx,
  const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic)
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
    assert( idx.max() == K.size()-1 );
    // - consistency
    assert( K.size() == G.size()      );
    assert( K.size() == epsy.shape(0) );
  #endif

  // start index
  size_t N = m_Cusp.size();

  // store material definition
  for ( size_t i = 0 ; i < K.size() ; ++i ) {
    std::vector<double> y(epsy.item(i), epsy.item(i)+epsy.shape(1));
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

inline void Matrix::setSmooth(const ArrS &I, const ArrS &idx,
  const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic)
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
    assert( idx.max() == K.size()-1 );
    // - consistency
    assert( K.size() == G.size()      );
    assert( K.size() == epsy.shape(0) );
  #endif

  // start index
  size_t N = m_Smooth.size();

  // store material definition
  for ( size_t i = 0 ; i < K.size() ; ++i ) {
    std::vector<double> y(epsy.item(i), epsy.item(i)+epsy.shape(1));
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

inline void Matrix::Sig(const ArrD &a_Eps, ArrD &a_Sig) const
{
  assert( a_Eps.shape()   == a_Sig.shape() );
  assert( a_Eps.shape(-1) == m_ncomp );
  assert( cppmat::del(a_Eps.shape(),-1) == m_type.shape() );

  // iterators
  auto i_Eps = a_Eps.begin();
  auto i_Sig = a_Sig.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp);
      T2s Sig;
      // compute
      switch ( m_type[i] )
      {
        case Type::Elastic: Sig = m_Elastic[m_index[i]].Sig(Eps); break;
        case Type::Cusp   : Sig = m_Cusp   [m_index[i]].Sig(Eps); break;
        case Type::Smooth : Sig = m_Smooth [m_index[i]].Sig(Eps); break;
        default: std::runtime_error("Unknown material");
      }
      // store
      std::copy(Sig.begin(), Sig.end(), i_Sig+i*m_ncomp);
    }
  }
}

// -------------------------------- compute stress for all entries ---------------------------------

inline ArrD Matrix::Sig(const ArrD &a_Eps) const
{
  // check input
  assert( a_Eps.shape(-1) == m_ncomp );
  assert( cppmat::del(a_Eps.shape(),-1) == m_type.shape() );

  // allocate output: matrix of tensors
  ArrD a_Sig(a_Eps.shape());

  // iterators
  auto i_Eps = a_Eps.begin();
  auto i_Sig = a_Sig.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp);
      T2s Sig;
      // compute
      switch ( m_type[i] )
      {
        case Type::Elastic: Sig = m_Elastic[m_index[i]].Sig(Eps); break;
        case Type::Cusp   : Sig = m_Cusp   [m_index[i]].Sig(Eps); break;
        case Type::Smooth : Sig = m_Smooth [m_index[i]].Sig(Eps); break;
        default: std::runtime_error("Unknown material");
      }
      // store
      std::copy(Sig.begin(), Sig.end(), i_Sig+i*m_ncomp);
    }
  }

  return a_Sig;
}


// -------------------------------- compute energy for all entries ---------------------------------

inline ArrD Matrix::energy(const ArrD &a_Eps) const
{
  // check input
  assert( a_Eps.shape(-1) == m_ncomp );
  assert( cppmat::del(a_Eps.shape(),-1) == m_type.shape() );

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD a_energy(cppmat::del(a_Eps.shape(),-1));

  // iterators
  auto i_Eps = a_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp);
      // compute/store
      switch ( m_type[i] )
      {
        case Type::Elastic: a_energy[i] = m_Elastic[m_index[i]].energy(Eps); break;
        case Type::Cusp   : a_energy[i] = m_Cusp   [m_index[i]].energy(Eps); break;
        case Type::Smooth : a_energy[i] = m_Smooth [m_index[i]].energy(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return a_energy;
}

// ------------------------- find the current yield strain for all entries -------------------------

inline ArrS Matrix::find(const ArrD &a_Eps) const
{
  // check input
  assert( a_Eps.shape(-1) == m_ncomp );
  assert( cppmat::del(a_Eps.shape(),-1) == m_type.shape() );

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrS a_idx(cppmat::del(a_Eps.shape(),-1));

  // iterators
  auto i_Eps = a_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp);
      // compute/store
      switch ( m_type[i] )
      {
        case Type::Elastic: a_idx[i] = m_Elastic[m_index[i]].find(Eps); break;
        case Type::Cusp   : a_idx[i] = m_Cusp   [m_index[i]].find(Eps); break;
        case Type::Smooth : a_idx[i] = m_Smooth [m_index[i]].find(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return a_idx;
}

// ----------------------------- get the yield strain for all entries ------------------------------

inline ArrD Matrix::epsy(const ArrS &a_idx) const
{
  // check input
  assert( a_idx.shape() == m_type.shape() );

  // allocate output: matrix of scalars
  ArrD a_epsy(a_idx.shape());

  // loop over all points
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_type.size() ; ++i )
  {
    switch ( m_type[i] )
    {
      case Type::Elastic: a_epsy[i] = m_Elastic[m_index[i]].epsy(a_idx[i]); break;
      case Type::Cusp   : a_epsy[i] = m_Cusp   [m_index[i]].epsy(a_idx[i]); break;
      case Type::Smooth : a_epsy[i] = m_Smooth [m_index[i]].epsy(a_idx[i]); break;
      default: std::runtime_error("Unknown material");
    }
  }

  return a_epsy;
}

// ----------------------- get the equivalent plastic strain for all entries -----------------------

inline ArrD Matrix::epsp(const ArrD &a_Eps) const
{
  // check input
  assert( a_Eps.shape(-1) == m_ncomp );
  assert( cppmat::del(a_Eps.shape(),-1) == m_type.shape() );

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD a_epsp(cppmat::del(a_Eps.shape(),-1));

  // iterators
  auto i_Eps = a_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp);
      // compute/store
      switch ( m_type[i] )
      {
        case Type::Elastic: a_epsp[i] = m_Elastic[m_index[i]].epsp(Eps); break;
        case Type::Cusp   : a_epsp[i] = m_Cusp   [m_index[i]].epsp(Eps); break;
        case Type::Smooth : a_epsp[i] = m_Smooth [m_index[i]].epsp(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return a_epsp;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

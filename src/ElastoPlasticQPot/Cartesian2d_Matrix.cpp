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

Matrix::Matrix(const std::vector<size_t> &shape)
{
  // resize type and index look-up
  m_type .resize(shape);
  m_index.resize(shape);

  // initialize everything to unassigned
  m_type.setConstant(Type::Unset);
}

// --------------------------------------------- type ----------------------------------------------

inline ArrS Matrix::type() const
{
  return m_type;
}

// --------------------------- check that a type has been set everywhere ---------------------------

inline void Matrix::check() const
{
  for ( size_t i = 0 ; i < m_type.size() ; ++i )
    if ( m_type[i] == Type::Unset )
      throw std::runtime_error("No type set for: " + cppmat::to_string(m_type.decompress(i)));
}

// ---------------------------------- set Elastic material points ----------------------------------

inline void Matrix::setElastic(
  const ArrS &I, double K, double G)
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

inline void Matrix::setCusp(
  const ArrS &I, double K, double G, const std::vector<double> &epsy, bool init_elastic)
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

inline void Matrix::setSmooth(
  const ArrS &I, double K, double G, const std::vector<double> &epsy, bool init_elastic)
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

// ---------------------------------- add Elastic material points ----------------------------------

inline void Matrix::addElastic(
  const MatS &index, const ColD &K, const ColD &G)
{
  // check input
  assert( index.cols() == m_type.ndim() );
  assert( index.rows() == K.size()      );
  assert( index.rows() == G.size()      );

  // loop over entries
  for ( size_t i = 0 ; i < index.rows() ; ++i )
  {
    // check
    assert( m_type.at(index.beginRow(i), index.endRow(i)) == Type::Unset );
    // set type and position in material vector
    m_type .at(index.beginRow(i), index.endRow(i)) = Type::Elastic;
    m_index.at(index.beginRow(i), index.endRow(i)) = m_Elastic.size();
    // store material definition
    m_Elastic.push_back(Elastic(K[i], G[i]));
  }
}

// ----------------------------------- add Cusp material points ------------------------------------

inline void Matrix::addCusp(
  const MatS &index, const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic)
{
  // check input
  assert( index.cols() == m_type.ndim() );
  assert( index.rows() == K.size()      );
  assert( index.rows() == G.size()      );
  assert( index.rows() == index.rows()  );

  // loop over entries
  for ( size_t i = 0 ; i < index.rows() ; ++i )
  {
    // check
    assert( m_type.at(index.beginRow(i), index.endRow(i)) == Type::Unset );
    // set type and position in material vector
    m_type .at(index.beginRow(i), index.endRow(i)) = Type::Cusp;
    m_index.at(index.beginRow(i), index.endRow(i)) = m_Cusp.size();
    // get yield strains
    std::vector<double> y(epsy.item(i), epsy.item(i)+epsy.cols());
    // store material definition
    m_Cusp.push_back(Cusp(K[i], G[i], y, init_elastic));
  }
}

// ---------------------------------- add Smooth material points -----------------------------------

inline void Matrix::addSmooth(
  const MatS &index, const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic)
{
  // check input
  assert( index.cols() == m_type.ndim() );
  assert( index.rows() == K.size()      );
  assert( index.rows() == G.size()      );
  assert( index.rows() == index.rows()  );

  // loop over entries
  for ( size_t i = 0 ; i < index.rows() ; ++i )
  {
    // check
    assert( m_type.at(index.beginRow(i), index.endRow(i)) == Type::Unset );
    // set type and position in material vector
    m_type .at(index.beginRow(i), index.endRow(i)) = Type::Smooth;
    m_index.at(index.beginRow(i), index.endRow(i)) = m_Smooth.size();
    // get yield strains
    std::vector<double> y(epsy.item(i), epsy.item(i)+epsy.cols());
    // store material definition
    m_Smooth.push_back(Smooth(K[i], G[i], y, init_elastic));
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
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp, i_Eps+(i+1)*m_ncomp);
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
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp, i_Eps+(i+1)*m_ncomp);
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
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp, i_Eps+(i+1)*m_ncomp);
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
      T2s Eps = T2s::Copy(i_Eps+i*m_ncomp, i_Eps+(i+1)*m_ncomp);
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

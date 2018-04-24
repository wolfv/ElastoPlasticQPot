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

// ------------------------------------------ mean strain ------------------------------------------

inline ArrD epsm(const ArrD &mat_Eps)
{
  // number of tensor components
  static const size_t ncomp = 3;

  // check input
  assert( mat_Eps.shape(-1) == ncomp );
  assert( mat_Eps.ndim()    >= 2     );

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD mat_epsm(cppmat::del(mat_Eps.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < mat_epsm.size() ; ++i )
    {
      // map from matrix of strains
      vT2s Eps = vT2s::Map(&mat_Eps[i*ncomp]);
      // compute the mean strain
      mat_epsm[i] = Eps.trace()/2.;
    }
  }

  return mat_epsm;
}

// ---------------------------------- equivalent strain deviator -----------------------------------

inline ArrD epsd(const ArrD &mat_Eps)
{
  // number of tensor components
  static const size_t ncomp = 3;

  // check input
  assert( mat_Eps.shape(-1) == ncomp );
  assert( mat_Eps.ndim()    >= 2     );

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD mat_epsd(cppmat::del(mat_Eps.shape(),-1));

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < mat_epsd.size() ; ++i )
    {
      // map from matrix of strains
      vT2s Eps  = vT2s::Map(&mat_Eps[i*ncomp]);
      // compute the strain deviator
      T2s  Epsd = Eps - Eps.trace()/2. * T2d::I();
      // compute the equivalent strain deviator
      mat_epsd[i] = std::sqrt(.5*Epsd.ddot(Epsd));
    }
  }

  return mat_epsd;
}

// ----------------------------------- mean & equivalent stress ------------------------------------

inline ArrD sigm(const ArrD &mat_Sig) { return epsm(mat_Sig); }
inline ArrD sigd(const ArrD &mat_Sig) { return epsd(mat_Sig); }

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
    // - number of dimensions
    assert( I.ndim() == m_type.ndim() );
    // - each shape
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( I.shape(i) == m_type(i) );
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
    // - number of dimensions
    assert( I.ndim() == m_type.ndim() );
    // - each shape
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( I.shape(i) == m_type(i) );
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
    // - number of dimensions
    assert( I.ndim() == m_type.ndim() );
    // - each shape
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( I.shape(i) == m_type(i) );
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

inline ArrD Matrix::stress(const ArrD &mat_Eps) const
{
  // check input
  #ifndef NDEBUG
    // number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // number of indices
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( mat_Eps.shape(i) == m_type.shape(i) );
  #endif

  // allocate output: matrix of tensors
  ArrD mat_Sig(mat_Eps.shape());

  // iterator to beginning of the matrices containing the strain and stress
  auto itr_Eps = mat_Eps.begin();
  auto itr_Sig = mat_Sig.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // temporary tensors
    T2s Sig, Eps;
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      std::copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp, Eps.data());
      // compute stress
      switch ( m_type[i] )
      {
        case Type::Elastic: Sig = m_Elastic[m_index[i]].stress(Eps); break;
        case Type::Cusp   : Sig = m_Cusp   [m_index[i]].stress(Eps); break;
        case Type::Smooth : Sig = m_Smooth [m_index[i]].stress(Eps); break;
        default: std::runtime_error("Unknown material");
      }
      // store stress to matrix
      std::copy(Sig.begin(), Sig.end(), itr_Sig+i*m_ncomp);
    }
  }

  return mat_Sig;
}

// -------------------------------- compute energy for all entries ---------------------------------

inline ArrD Matrix::energy(const ArrD &mat_Eps) const
{
  // check input
  #ifndef NDEBUG
    // number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // number of indices
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( mat_Eps.shape(i) == m_type.shape(i) );
  #endif

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD mat_energy(cppmat::del(mat_Eps.shape(),-1));

  // iterator to beginning of the matrices containing the strain
  auto itr_Eps = mat_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp);
      // compute energy and store to matrix
      switch ( m_type[i] )
      {
        case Type::Elastic: mat_energy[i] = m_Elastic[m_index[i]].energy(Eps); break;
        case Type::Cusp   : mat_energy[i] = m_Cusp   [m_index[i]].energy(Eps); break;
        case Type::Smooth : mat_energy[i] = m_Smooth [m_index[i]].energy(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return mat_energy;
}

// ------------------------- find the current yield strain for all entries -------------------------

inline ArrS Matrix::find(const ArrD &mat_Eps) const
{
  // check input
  #ifndef NDEBUG
    // number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // number of indices
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( mat_Eps.shape(i) == m_type.shape(i) );
  #endif

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrS mat_idx(cppmat::del(mat_Eps.shape(),-1));

  // iterator to beginning of the matrices containing the strain
  auto itr_Eps = mat_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp);
      // find index and store to matrix
      switch ( m_type[i] )
      {
        case Type::Elastic: mat_idx[i] = m_Elastic[m_index[i]].find(Eps); break;
        case Type::Cusp   : mat_idx[i] = m_Cusp   [m_index[i]].find(Eps); break;
        case Type::Smooth : mat_idx[i] = m_Smooth [m_index[i]].find(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return mat_idx;
}

// ----------------------------- get the yield strain for all entries ------------------------------

inline ArrD Matrix::epsy(const ArrS &mat_idx) const
{
  // check input
  #ifndef NDEBUG
    // number of dimensions
    assert( mat_idx.ndim() == m_type.ndim() );
    // number of indices
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( mat_idx.shape(i) == m_type.shape(i) );
  #endif

  // allocate output: matrix of scalars
  ArrD mat_epsy(mat_idx.shape());

  // loop over all points
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_type.size() ; ++i )
  {
    switch ( m_type[i] )
    {
      case Type::Elastic: mat_epsy[i] = m_Elastic[m_index[i]].epsy(mat_idx[i]); break;
      case Type::Cusp   : mat_epsy[i] = m_Cusp   [m_index[i]].epsy(mat_idx[i]); break;
      case Type::Smooth : mat_epsy[i] = m_Smooth [m_index[i]].epsy(mat_idx[i]); break;
      default: std::runtime_error("Unknown material");
    }
  }

  return mat_epsy;
}

// ----------------------- get the equivalent plastic strain for all entries -----------------------

inline ArrD Matrix::epsp(const ArrD &mat_Eps) const
{
  // check input
  #ifndef NDEBUG
    // number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // number of indices
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( mat_Eps.shape(i) == m_type.shape(i) );
  #endif

  // allocate output: matrix of scalars (shape of the input matrix, without index)
  ArrD mat_epsp(cppmat::del(mat_Eps.shape(),-1));

  // iterator to beginning of the matrices containing the strain
  auto itr_Eps = mat_Eps.begin();

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // copy strain from matrix
      T2s Eps = T2s::Copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp);
      // compute plastic strain and store to matrix
      switch ( m_type[i] )
      {
        case Type::Elastic: mat_epsp[i] = m_Elastic[m_index[i]].epsp(Eps); break;
        case Type::Cusp   : mat_epsp[i] = m_Cusp   [m_index[i]].epsp(Eps); break;
        case Type::Smooth : mat_epsp[i] = m_Smooth [m_index[i]].epsp(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return mat_epsp;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

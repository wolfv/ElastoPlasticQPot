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
    // - temporary tensor
    vT2s Eps;

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < mat_epsm.size() ; ++i )
    {
      // -- map from matrix of strains
      Eps.map(&mat_Eps[i*ncomp]);
      // -- compute the mean strain
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
    // - temporary tensors
    vT2s Eps;
    T2s  Epsd;
    T2d  I = cm::identity2<double>();

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < mat_epsd.size() ; ++i )
    {
      // -- map from matrix of strains
      Eps.map(&mat_Eps[i*ncomp]);
      // -- compute the strain deviator
      Epsd = Eps - Eps.trace()/2. * I;
      // -- compute the equivalent strain deviator
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

// ---------------------------------- add Elastic material points ----------------------------------

inline void Matrix::addElastic(
  const ArrS &index, const ArrD &K, const ArrD &G)
{
  // check input
  assert( index.ndim()   == 2              );
  assert( index.shape(1) == m_type.ndim()  );
  assert( K.ndim()       == 1              );
  assert( K.size()       == index.shape(0) );
  assert( G.ndim()       == 1              );
  assert( G.size()       == index.shape(0) );

  // number of entries, number matrix dimensions
  size_t rows = index.shape(0);
  size_t cols = index.shape(1);

  // loop over entries
  for ( size_t i = 0 ; i < rows ; ++i )
  {
    // - check
    assert( m_type.at(index.item(i), index.item(i)+cols) == Type::Unset );
    // - set type and position in material vector
    m_type .at(index.item(i), index.item(i)+cols) = Type::Elastic;
    m_index.at(index.item(i), index.item(i)+cols) = m_Elastic.size();
    // - store material definition
    m_Elastic.push_back(Elastic(K[i], G[i]));
  }
}

// ----------------------------------- add Cusp material points ------------------------------------

inline void Matrix::addCusp(
  const ArrS &index, const ArrD &K, const ArrD &G, const ArrD &epsy, bool init_elastic)
{
  // check input
  assert( index.ndim()   == 2              );
  assert( index.shape(1) == m_type.ndim()  );
  assert( K.ndim()       == 1              );
  assert( K.size()       == index.shape(0) );
  assert( G.ndim()       == 1              );
  assert( G.size()       == index.shape(0) );
  assert( epsy.ndim()    == 2              );
  assert( epsy.shape(0)  == index.shape(0) );

  // number of entries, number matrix dimensions, number of yield strains per entry
  size_t rows = index.shape(0);
  size_t cols = index.shape(1);
  size_t ny   = epsy .shape(1);

  // loop over entries
  for ( size_t i = 0 ; i < rows ; ++i )
  {
    // - check
    assert( m_type.at(index.item(i), index.item(i)+cols) == Type::Unset );
    // - set type and position in material vector
    m_type .at(index.item(i), index.item(i)+cols) = Type::Cusp;
    m_index.at(index.item(i), index.item(i)+cols) = m_Cusp.size();
    // - get yield strains
    std::vector<double> y(epsy.item(i), epsy.item(i)+ny);
    // - store material definition
    m_Cusp.push_back(Cusp(K[i], G[i], y, init_elastic));
  }
}

// ---------------------------------- add Smooth material points -----------------------------------

inline void Matrix::addSmooth(
  const ArrS &index, const ArrD &K, const ArrD &G, const ArrD &epsy, bool init_elastic)
{
  // check input
  assert( index.ndim()   == 2              );
  assert( index.shape(1) == m_type.ndim()  );
  assert( K.ndim()       == 1              );
  assert( K.size()       == index.shape(0) );
  assert( G.ndim()       == 1              );
  assert( G.size()       == index.shape(0) );
  assert( epsy.ndim()    == 2              );
  assert( epsy.shape(0)  == index.shape(0) );

  // number of entries, number matrix dimensions, number of yield strains per entry
  size_t rows = index.shape(0);
  size_t cols = index.shape(1);
  size_t ny   = epsy .shape(1);

  // loop over entries
  for ( size_t i = 0 ; i < rows ; ++i )
  {
    // - check
    assert( m_type.at(index.item(i), index.item(i)+cols) == Type::Unset );
    // - set type and position in material vector
    m_type .at(index.item(i), index.item(i)+cols) = Type::Smooth;
    m_index.at(index.item(i), index.item(i)+cols) = m_Smooth.size();
    // - get yield strains
    std::vector<double> y(epsy.item(i), epsy.item(i)+ny);
    // - store material definition
    m_Smooth.push_back(Smooth(K[i], G[i], y, init_elastic));
  }
}

// -------------------------------- compute stress for all entries ---------------------------------

inline ArrD Matrix::stress(const ArrD &mat_Eps) const
{
  // check input
  #ifndef NDEBUG
    // - number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // - number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // - number of indices
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
    // - temporary tensors
    T2s Sig, Eps;

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp, Eps.data());
      // -- compute stress
      switch ( m_type[i] )
      {
        case Type::Elastic: Sig = m_Elastic[m_index[i]].stress(Eps); break;
        case Type::Cusp   : Sig = m_Cusp   [m_index[i]].stress(Eps); break;
        case Type::Smooth : Sig = m_Smooth [m_index[i]].stress(Eps); break;
        default: std::runtime_error("Unknown material");
      }
      // -- store stress to matrix
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
    // - number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // - number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // - number of indices
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
    // - temporary tensor
    T2s Eps;

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp, Eps.data());
      // -- compute energy and store to matrix
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
    // - number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // - number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // - number of indices
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
    // - temporary tensor
    T2s Eps;

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp, Eps.data());
      // -- find index and store to matrix
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
    // - number of dimensions
    assert( mat_idx.ndim() == m_type.ndim() );
    // - number of indices
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
    // - number of tensor-components
    assert( mat_Eps.shape(-1) == m_ncomp );
    // - number of dimensions
    assert( mat_Eps.ndim()-1 == m_type.ndim() );
    // - number of indices
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
    // - temporary tensor
    T2s Eps;

    // - loop over all points
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(itr_Eps+i*m_ncomp, itr_Eps+(i+1)*m_ncomp, Eps.data());
      // -- compute plastic strain and store to matrix
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

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

inline ArrD epsm(const ArrD &strain)
{
  // number of tensor components
  static const size_t ncomp = 3;

  // check input
  assert( strain.shape(-1) == ncomp );
  assert( strain.ndim()    >= 2     );

  // get shape of the input matrix
  std::vector<size_t> shape = strain.shape();

  // remove tensor-components from the shape
  shape.erase(shape.end()-1, shape.end());

  // allocate output as matrix of scalars
  ArrD out(shape);

  // start threads
  #pragma omp parallel
  {
    // - per thread; temporary variables
    vT2s Eps;

    // - per thread; compute
    #pragma omp for
    for ( size_t i = 0 ; i < out.size() ; ++i )
    {
      // -- map from matrix of strains
      Eps.map(&strain[i*ncomp]);
      // -- compute the mean strain
      out[i] = Eps.trace()/2.;
    }
  }

  return out;
}

// ---------------------------------- equivalent strain deviator -----------------------------------

inline ArrD epsd(const ArrD &strain)
{
  // number of tensor components
  static const size_t ncomp = 3;

  // check input
  assert( strain.shape(-1) == ncomp );

  // get shape of the input matrix
  std::vector<size_t> shape = strain.shape();

  // remove tensor-components from the shape
  shape.erase(shape.end()-1, shape.end());

  // allocate output as matrix of scalars
  ArrD out(shape);

  // start threads
  #pragma omp parallel
  {
    // - per thread; temporary variables
    double epsm;
    vT2s   Eps;
    T2s    Epsd;
    T2d    I = cm::identity2<double>();

    // - per thread; compute
    #pragma omp for
    for ( size_t i = 0 ; i < out.size() ; ++i )
    {
      // -- map from matrix of strains
      Eps.map(&strain[i*ncomp]);
      // -- compute the strain deviator
      epsm = Eps.trace()/2.;
      Epsd = Eps - epsm*I;
      // -- compute the equivalent strain deviator
      out[i] = std::sqrt(.5*Epsd.ddot(Epsd));
    }
  }

  return out;
}

// ----------------------------------- mean & equivalent stress ------------------------------------

inline ArrD sigm(const ArrD &stress) { return epsm(stress); }
inline ArrD sigd(const ArrD &stress) { return epsd(stress); }

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

inline ArrD Matrix::stress(const ArrD &matrix) const
{
  // check input
  #ifndef NDEBUG
    // - #tensor-components
    assert( matrix.shape(-1) == m_ncomp );
    // - matrix shape: number of indices
    assert( matrix.ndim()-1 == m_type.ndim() );
    // - matrix shape: per index
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( matrix.shape(i) == m_type.shape(i) );
  #endif

  // allocate output
  ArrD out(matrix.shape());

  // iterator to beginning of the matrices containing the strain and stress
  auto strain = matrix.begin();
  auto stress = out   .begin();

  // start threads
  #pragma omp parallel
  {
    // - per thread; stress/strain tensor
    T2s Sig, Eps;

    // - per thread; constitutive response
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(strain+i*m_ncomp, strain+(i+1)*m_ncomp, Eps.data());
      // -- compute stress
      switch ( m_type[i] )
      {
        case Type::Elastic: Sig = m_Elastic[m_index[i]].stress(Eps); break;
        case Type::Cusp   : Sig = m_Cusp   [m_index[i]].stress(Eps); break;
        case Type::Smooth : Sig = m_Smooth [m_index[i]].stress(Eps); break;
        default: std::runtime_error("Unknown material");
      }
      // -- store stress to matrix
      std::copy(Sig.begin(), Sig.end(), stress+i*m_ncomp);
    }
  }

  return out;
}

// -------------------------------- compute energy for all entries ---------------------------------

inline ArrD Matrix::energy(const ArrD &matrix) const
{
  // check input
  #ifndef NDEBUG
    // - #tensor-components
    assert( matrix.shape(-1) == m_ncomp );
    // - matrix shape: number of indices
    assert( matrix.ndim()-1 == m_type.ndim() );
    // - matrix shape: per index
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( matrix.shape(i) == m_type.shape(i) );
  #endif

  // get shape of the input matrix
  std::vector<size_t> shape = matrix.shape();

  // convert shape to corresponding shape without the tensor-components
  shape.erase(shape.end()-1, shape.end());

  // allocate output
  ArrD out(shape);

  // iterator to beginning of the matrices containing the strain
  auto strain = matrix.begin();

  // start threads
  #pragma omp parallel
  {
    // - per thread; strain tensor
    T2s Eps;

    // - per thread; compute energy
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(strain+i*m_ncomp, strain+(i+1)*m_ncomp, Eps.data());
      // -- compute energy and store to matrix
      switch ( m_type[i] )
      {
        case Type::Elastic: out[i] = m_Elastic[m_index[i]].energy(Eps); break;
        case Type::Cusp   : out[i] = m_Cusp   [m_index[i]].energy(Eps); break;
        case Type::Smooth : out[i] = m_Smooth [m_index[i]].energy(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return out;
}

// ------------------------- find the current yield strain for all entries -------------------------

inline ArrS Matrix::find(const ArrD &matrix) const
{
  // check input
  #ifndef NDEBUG
    // - #tensor-components
    assert( matrix.shape(-1) == m_ncomp );
    // - matrix shape: number of indices
    assert( matrix.ndim()-1 == m_type.ndim() );
    // - matrix shape: per index
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( matrix.shape(i) == m_type.shape(i) );
  #endif

  // get shape of the input matrix
  std::vector<size_t> shape = matrix.shape();

  // convert shape to corresponding shape without the tensor-components
  shape.erase(shape.end()-1, shape.end());

  // allocate output
  ArrS out(shape);

  // iterator to beginning of the matrices containing the strain
  auto strain = matrix.begin();

  // start threads
  #pragma omp parallel
  {
    // - per thread; strain tensor
    T2s Eps;

    // - per thread; find index
    #pragma omp for
    for ( size_t i = 0 ; i < m_type.size() ; ++i )
    {
      // -- copy strain from matrix
      std::copy(strain+i*m_ncomp, strain+(i+1)*m_ncomp, Eps.data());
      // -- find index and store to matrix
      switch ( m_type[i] )
      {
        case Type::Elastic: out[i] = m_Elastic[m_index[i]].find(Eps); break;
        case Type::Cusp   : out[i] = m_Cusp   [m_index[i]].find(Eps); break;
        case Type::Smooth : out[i] = m_Smooth [m_index[i]].find(Eps); break;
        default: std::runtime_error("Unknown material");
      }
    }
  }

  return out;
}

// ----------------------------- get the yield strain for all entries ------------------------------

inline ArrD Matrix::eps_y(const ArrS &matrix) const
{
  // check input
  #ifndef NDEBUG
    // - #tensor-components
    assert( matrix.shape(-1) == m_ncomp );
    // - matrix shape: number of indices
    assert( matrix.ndim()-1 == m_type.ndim() );
    // - matrix shape: per index
    for ( size_t i = 0 ; i < m_type.ndim() ; ++i )
      assert( matrix.shape(i) == m_type.shape(i) );
  #endif

  // allocate output
  ArrD out(matrix.shape());

  // start threads
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_type.size() ; ++i )
  {
    // -- get yield strain and store to matrix
    switch ( m_type[i] )
    {
      case Type::Elastic: out[i] = m_Elastic[m_index[i]].eps_y(matrix[i]); break;
      case Type::Cusp   : out[i] = m_Cusp   [m_index[i]].eps_y(matrix[i]); break;
      case Type::Smooth : out[i] = m_Smooth [m_index[i]].eps_y(matrix[i]); break;
      default: std::runtime_error("Unknown material");
    }
  }

  return out;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

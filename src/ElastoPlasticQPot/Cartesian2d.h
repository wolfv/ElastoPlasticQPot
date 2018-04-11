/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef ELASTOPLASTICQPOT_CARTESIAN2D_H
#define ELASTOPLASTICQPOT_CARTESIAN2D_H

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// ================================ ElastoPlasticQPot::Cartesian2d =================================

namespace ElastoPlasticQPot {
namespace Cartesian2d {

// -------------------------------------------------------------------------------------------------

namespace cm = cppmat::cartesian2d;

typedef cppmat::cartesian2d::tensor2 <double> T2;
typedef cppmat::cartesian2d::tensor2s<double> T2s;
typedef cppmat::cartesian2d::tensor2d<double> T2d;

// -------------------------- equivalent stress/strain (Cartesian2d.cpp) ---------------------------

inline double epsm(const T2s &Eps);
inline double epsd(const T2s &Eps);

// ----------------------- equivalent stress/strain (Cartesian2d_Matrix.cpp) -----------------------

inline ArrD epsm(const ArrD &Eps);
inline ArrD epsd(const ArrD &Eps);

// ---------------------- material point - elastic (Cartesian2d_Elastic.cpp) -----------------------

class Elastic
{
private:

  double m_K; // bulk  modulus
  double m_G; // shear modulus

public:

  // constructor
  Elastic(){};
  Elastic(double K, double G);

  // compute stress or energy at "Eps"
  T2s    stress(const T2s &Eps) const;
  double energy(const T2s &Eps) const;

  // find index of the current yield strain, return yield strain at this index
  size_t find (const T2s &Eps ) const;
  size_t find (double     epsd) const;
  double eps_y(size_t     i   ) const;
};

// -------------------- material point - cusp potential (Cartesian2d_Cusp.cpp) ---------------------

class Cusp
{
private:

  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:

  // constructor
  Cusp(){};
  Cusp(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // compute stress or energy at "Eps"
  T2s    stress(const T2s &Eps) const;
  double energy(const T2s &Eps) const;

  // find index of the current yield strain, return yield strain at this index
  size_t find (const T2s &Eps ) const;
  size_t find (double     epsd) const;
  double eps_y(size_t     i   ) const;
};

// ------------------ material point - smooth potential (Cartesian2d_Smooth.cpp) -------------------

class Smooth
{
private:

  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:

  // constructor
  Smooth(){};
  Smooth(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // compute stress or energy at "Eps"
  T2s    stress(const T2s &Eps) const;
  double energy(const T2s &Eps) const;

  // find index of the current yield strain, return yield strain at this index
  size_t find (const T2s &Eps ) const;
  size_t find (double     epsd) const;
  double eps_y(size_t     i   ) const;
};

// ----------------------- enumerator to switch between material definitions -----------------------

struct Type {
  enum Value {
    Unset,
    Elastic,
    Cusp,
    Smooth,
    PlanarCusp,
    PlanarSmooth,
  };
};

// ---------------------- Matrix of material points (Cartesian2d_Matrix.cpp) -----------------------

class Matrix
{
private:

  // material vectors
  std::vector<Elastic> m_Elastic;
  std::vector<Cusp>    m_Cusp;
  std::vector<Smooth>  m_Smooth;

  // identifiers for each matrix entry
  ArrS m_type;  // type (e.g. "Type::Elastic")
  ArrS m_index; // index from the relevant material vector (e.g. "m_Elastic")

  // dimensions
  static const size_t m_ndim=2;   // number of dimensions
  static const size_t m_ncomp=3;  // number of tensor components

public:

  // constructor
  Matrix(){};
  Matrix(const std::vector<size_t> &shape);

  // return type
  ArrS type() const;

  // add material definitions
  void addElastic(const ArrS &index, const ArrD &K, const ArrD &G);
  // -
  void addCusp(
    const ArrS &index, const ArrD &K, const ArrD &G, const ArrD &epsy, bool init_elastic=true
  );
  // -
  void addSmooth(
    const ArrS &index, const ArrD &K, const ArrD &G, const ArrD &epsy, bool init_elastic=true
  );

  // compute stress or energy at "Eps"
  ArrD stress(const ArrD &Eps) const;
  ArrD energy(const ArrD &Eps) const;

  // find index of the current yield strain, return yield strain at this index
  ArrS find  (const ArrD &Eps) const;
  ArrD eps_y (const ArrS &i  ) const;

};

// -------------------------------------------------------------------------------------------------

}}

// -------------------------------------------------------------------------------------------------

#endif

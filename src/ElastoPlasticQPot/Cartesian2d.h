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

// --------------------------------------------- alias ---------------------------------------------

typedef cppmat::array <size_t> ArrS;
typedef cppmat::array <double> ArrD;
typedef cppmat::matrix<size_t> MatS;
typedef cppmat::matrix<double> MatD;
typedef cppmat::vector<size_t> ColS;
typedef cppmat::vector<double> ColD;

typedef cppmat::tiny::cartesian::tensor2 <double,2> T2;
typedef cppmat::tiny::cartesian::tensor2s<double,2> T2s;
typedef cppmat::tiny::cartesian::tensor2d<double,2> T2d;

// -------------------------- equivalent stress/strain (Cartesian2d.cpp) ---------------------------

// mean
inline double sigm(const T2s &Sig);
inline double epsm(const T2s &Eps);

// equivalent deviator
inline double sigd(const T2s &Sig);
inline double epsd(const T2s &Eps);

// deviator
inline T2s Sigd(const T2s &Sig);
inline T2s Epsd(const T2s &Eps);

// ----------------------- equivalent stress/strain (Cartesian2d_Matrix.cpp) -----------------------

// mean
inline ArrD sigm(const ArrD &a_Sig);
inline ArrD epsm(const ArrD &a_Eps);

// equivalent deviator
inline ArrD sigd(const ArrD &a_Sig);
inline ArrD epsd(const ArrD &a_Eps);

// deviator
inline ArrD Sigd(const ArrD &a_Sig);
inline ArrD Epsd(const ArrD &a_Eps);

// ---------------------- material point - elastic (Cartesian2d_Elastic.cpp) -----------------------

class Elastic
{
private:

  // parameters
  double m_K; // bulk  modulus
  double m_G; // shear modulus

public:

  // constructor
  Elastic() = default;
  Elastic(double K, double G);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double K() const;
  double G() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
};

// -------------------- material point - cusp potential (Cartesian2d_Cusp.cpp) ---------------------

class Cusp
{
private:

  // parameters
  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:

  // constructor
  Cusp() = default;
  Cusp(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double K() const;
  double G() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
};

// ------------------ material point - smooth potential (Cartesian2d_Smooth.cpp) -------------------

class Smooth
{
private:

  // parameters
  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:

  // constructor
  Smooth() = default;
  Smooth(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double K() const;
  double G() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
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

// ------------- matrix of material points of different types (Cartesian2d_Matrix.cpp) -------------

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
  Matrix() = default;
  Matrix(const std::vector<size_t> &shape);

  // return shape
  std::vector<size_t> shape() const;
  size_t shape(size_t i) const;

  // return type
  ArrS type() const;

  // return plastic yes/no
  ArrS isPlastic() const;

  // parameters
  ArrD K() const;
  ArrD G() const;

  // check that a type has been set everywhere
  void check() const;

  // set material definition for a batch of points
  // -
  void setElastic(const ArrS &I,
    double K, double G);
  // -
  void setCusp(const ArrS &I,
    double K, double G, const std::vector<double> &epsy, bool init_elastic=true);
  // -
  void setSmooth(const ArrS &I,
    double K, double G, const std::vector<double> &epsy, bool init_elastic=true);

  // set material definition for a batch of points
  // -
  void setElastic(const ArrS &I, const ArrS &idx,
    const ColD &K, const ColD &G);
  // -
  void setCusp(const ArrS &I, const ArrS &idx,
    const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic=true);
  // -
  void setSmooth(const ArrS &I, const ArrS &idx,
    const ColD &K, const ColD &G, const MatD &epsy, bool init_elastic=true);

  // stress
  ArrD Sig(const ArrD &a_Eps) const;

  // energy
  ArrD energy(const ArrD &a_Eps) const;

  // index of the current yield strain
  ArrS find(const ArrD &a_Eps) const;

  // certain yield strain
  ArrD epsy(const ArrS &a_idx) const;

  // equivalent plastic strain
  ArrD epsp(const ArrD &a_Eps) const;

};

// -------------------------------------------------------------------------------------------------

}}

// -------------------------------------------------------------------------------------------------

#endif

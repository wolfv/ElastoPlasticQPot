/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef XELASTOPLASTICQPOT_CARTESIAN2D_H
#define XELASTOPLASTICQPOT_CARTESIAN2D_H

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot.h"

// ================================ ElastoPlasticQPot::Cartesian2d =================================

namespace xElastoPlasticQPot {
namespace Cartesian2d {

// --------------------------------------------- alias ---------------------------------------------

using T2s = xt::xtensor_fixed<double, xt::xshape<2,2>>;


// ---------------------------------------- tensor algebra -----------------------------------------

template <class T>
inline double trace(const T &A);

template <class T>
inline double ddot (const T &A, const T &B);

inline T2s I();

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

// no allocation
inline void sigm(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigm);
inline void epsm(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsm);
inline void sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigd);
inline void epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsd);
inline void Sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,4> &a_Sigd);
inline void Epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Epsd);

// return allocated result
inline xt::xtensor<double,2> sigm(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,2> epsm(const xt::xtensor<double,4> &a_Eps);
inline xt::xtensor<double,2> sigd(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,2> epsd(const xt::xtensor<double,4> &a_Eps);
inline xt::xtensor<double,4> Sigd(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,4> Epsd(const xt::xtensor<double,4> &a_Eps);

// compute maximum, avoiding allocation
inline double sigm_max(const xt::xtensor<double,4> &a_Sig);
inline double epsm_max(const xt::xtensor<double,4> &a_Eps);
inline double sigd_max(const xt::xtensor<double,4> &a_Sig);
inline double epsd_max(const xt::xtensor<double,4> &a_Eps);

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
  xt::xtensor<size_t,2> m_type;  // type (e.g. "Type::Elastic")
  xt::xtensor<size_t,2> m_index; // index from the relevant material vector (e.g. "m_Elastic")

public:

  // constructor
  Matrix() = default;
  Matrix(const std::vector<size_t> &shape);

  // return shape
  std::vector<size_t> shape() const;
  size_t shape(size_t i) const;

  // return type
  xt::xtensor<size_t,2> type() const;

  // return plastic yes/no
  xt::xtensor<size_t,2> isPlastic() const;

  // parameters
  xt::xtensor<double,2> K() const;
  xt::xtensor<double,2> G() const;

  // check that a type has been set everywhere
  void check() const;

  // set material definition for a batch of points
  // -
  void setElastic(const xt::xtensor<size_t,2> &I,
    double K, double G);
  // -
  void setCusp(const xt::xtensor<size_t,2> &I,
    double K, double G, const std::vector<double> &epsy, bool init_elastic=true);
  // -
  void setSmooth(const xt::xtensor<size_t,2> &I,
    double K, double G, const std::vector<double> &epsy, bool init_elastic=true);

  // set material definition for a batch of points
  // -
  void setElastic(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G);
  // -
  void setCusp(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G,
    const xt::xtensor<double,2> &epsy, bool init_elastic=true);
  // -
  void setSmooth(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &K, const xt::xtensor<double,1> &G,
    const xt::xtensor<double,2> &epsy, bool init_elastic=true);

  // compute (no allocation)
  void Sig   (const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Sig   ) const;
  void energy(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_energy) const;
  void find  (const xt::xtensor<double,4> &a_Eps, xt::xtensor<size_t,2> &a_find  ) const;
  void epsy  (const xt::xtensor<size_t,2> &a_idx, xt::xtensor<double,2> &a_epsy  ) const;
  void epsp  (const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsp  ) const;

  // compute (return allocated result)
  xt::xtensor<double,4> Sig   (const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<double,2> energy(const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<size_t,2> find  (const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<double,2> epsy  (const xt::xtensor<size_t,2> &a_idx) const;
  xt::xtensor<double,2> epsp  (const xt::xtensor<double,4> &a_Eps) const;

};

// -------------------------------------------------------------------------------------------------

}}

// -------------------------------------------------------------------------------------------------

#endif

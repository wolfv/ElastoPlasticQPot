/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>
// #include <xtensor-python/pytensor.hpp>

#include "ElastoPlasticQPot.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;

// ======================================= ElastoPlasticQPot =======================================

PYBIND11_MODULE(ElastoPlasticQPot, m) {

m.doc() = "Elasto-plastic material models";

// ================================ ElastoPlasticQPot::Cartesian2d =================================

{

// create submodule
py::module sm = m.def_submodule("Cartesian2d", "2d Cartesian coordinates");

// abbreviate name-space
namespace SM = ElastoPlasticQPot::Cartesian2d;

// abbreviate types(s)
typedef SM::T2s T2s;

// -------------------------------------------------------------------------------------------------

sm.def("epsm", py::overload_cast<const T2s &>(&SM::epsm), "Hydrostatic strain" , py::arg("Eps"));
sm.def("epsd", py::overload_cast<const T2s &>(&SM::epsd), "Eq. strain deviator", py::arg("Eps"));
sm.def("Epsd", py::overload_cast<const T2s &>(&SM::Epsd), "Strain deviator"    , py::arg("Eps"));

sm.def("sigm", py::overload_cast<const T2s &>(&SM::sigm), "Hydrostatic stress" , py::arg("Sig"));
sm.def("sigd", py::overload_cast<const T2s &>(&SM::sigd), "Eq. stress deviator", py::arg("Sig"));
sm.def("Sigd", py::overload_cast<const T2s &>(&SM::Sigd), "Stress deviator"    , py::arg("Sig"));

sm.def("epsm", py::overload_cast<const xt::xtensor<double,4> &>(&SM::epsm), "Hydrostatic strain" , py::arg("a_Eps"));
sm.def("epsd", py::overload_cast<const xt::xtensor<double,4> &>(&SM::epsd), "Eq. strain deviator", py::arg("a_Eps"));
sm.def("Epsd", py::overload_cast<const xt::xtensor<double,4> &>(&SM::Epsd), "Strain deviator"    , py::arg("a_Eps"));

sm.def("sigm", py::overload_cast<const xt::xtensor<double,4> &>(&SM::sigm), "Hydrostatic stress" , py::arg("a_Sig"));
sm.def("sigd", py::overload_cast<const xt::xtensor<double,4> &>(&SM::sigd), "Eq. stress deviator", py::arg("a_Sig"));
sm.def("Sigd", py::overload_cast<const xt::xtensor<double,4> &>(&SM::Sigd), "Stress deviator"    , py::arg("a_Sig"));

// -------------------------------------------------------------------------------------------------

py::class_<SM::Elastic>(sm, "Elastic")
  // constructor
  .def(
    py::init<double,double>(),
    "Elastic material",
    py::arg("K"),
    py::arg("G")
  )
  // methods
  .def("Sig"   , &SM::Elastic::Sig   , py::arg("Eps"))
  .def("energy", &SM::Elastic::energy, py::arg("Eps"))
  .def("epsy"  , &SM::Elastic::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Elastic::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&SM::Elastic::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Elastic::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&SM::Elastic::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const SM::Elastic &){
    return "<ElastoPlasticQPot.Cartesian2d.Elastic>"; });

// -------------------------------------------------------------------------------------------------

py::class_<SM::Cusp>(sm, "Cusp")
  // constructor
  .def(
    py::init<double,double,const std::vector<double>&, bool>(),
    "Cusp material",
    py::arg("K"),
    py::arg("G"),
    py::arg("epsy"),
    py::arg("init_elastic")=true
  )
  // methods
  .def("Sig"   , &SM::Cusp::Sig   , py::arg("Eps"))
  .def("energy", &SM::Cusp::energy, py::arg("Eps"))
  .def("epsy"  , &SM::Cusp::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Cusp::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&SM::Cusp::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Cusp::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&SM::Cusp::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const SM::Cusp &){
    return "<ElastoPlasticQPot.Cartesian2d.Cusp>"; });

// -------------------------------------------------------------------------------------------------

py::class_<SM::Smooth>(sm, "Smooth")
  // constructor
  .def(
    py::init<double,double,const std::vector<double>&, bool>(),
    "Smooth material",
    py::arg("K"),
    py::arg("G"),
    py::arg("epsy"),
    py::arg("init_elastic")=true
  )
  // methods
  .def("Sig"   , &SM::Smooth::Sig   , py::arg("Eps"))
  .def("energy", &SM::Smooth::energy, py::arg("Eps"))
  .def("epsy"  , &SM::Smooth::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Smooth::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&SM::Smooth::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Smooth::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&SM::Smooth::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const SM::Smooth &){
    return "<ElastoPlasticQPot.Cartesian2d.Smooth>"; });

// -------------------------------------------------------------------------------------------------

py::module smm = sm.def_submodule("Type", "Type enumerator");

py::enum_<SM::Type::Value>(smm, "Type")
    .value("Unset"       , SM::Type::Unset)
    .value("Elastic"     , SM::Type::Elastic)
    .value("Cusp"        , SM::Type::Cusp)
    .value("Smooth"      , SM::Type::Smooth)
    .value("PlanarCusp"  , SM::Type::PlanarCusp)
    .value("PlanarSmooth", SM::Type::PlanarSmooth)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<SM::Matrix>(sm, "Matrix")
  // constructor
  .def(
    py::init<const std::vector<size_t>&>(),
    "Matrix of materials",
    py::arg("shape")
  )
  // methods
  .def("setElastic", py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double                                                            >(&SM::Matrix::setElastic),py::arg("I"),               py::arg("K"),py::arg("G"))
  .def("setCusp"   , py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double,                        const std::vector<double>   &, bool>(&SM::Matrix::setCusp   ),py::arg("I"),               py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double,                        const std::vector<double>   &, bool>(&SM::Matrix::setSmooth ),py::arg("I"),               py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setElastic", py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &                                     >(&SM::Matrix::setElastic),py::arg("I"),py::arg("idx"),py::arg("K"),py::arg("G"))
  .def("setCusp"   , py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &, const xt::xtensor<double,2> &, bool>(&SM::Matrix::setCusp   ),py::arg("I"),py::arg("idx"),py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &, const xt::xtensor<double,2> &, bool>(&SM::Matrix::setSmooth ),py::arg("I"),py::arg("idx"),py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("shape"     , py::overload_cast<size_t>(&SM::Matrix::shape, py::const_))
  .def("shape"     , py::overload_cast<>      (&SM::Matrix::shape, py::const_))
  .def("type"      ,                           &SM::Matrix::type)
  .def("Sig"       , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Sig   , py::const_), py::arg("a_Eps"))
  .def("energy"    , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::energy, py::const_), py::arg("a_Eps"))
  .def("find"      , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::find  , py::const_), py::arg("a_Eps"))
  .def("epsy"      , py::overload_cast<const xt::xtensor<size_t,2> &>(&SM::Matrix::epsy  , py::const_), py::arg("a_idx"))
  .def("epsp"      , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::epsp  , py::const_), py::arg("a_Eps"))
  // print to screen
  .def("__repr__", [](const SM::Matrix &){
    return "<ElastoPlasticQPot.Cartesian2d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}


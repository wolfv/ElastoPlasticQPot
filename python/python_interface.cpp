/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "../src/ElastoPlasticQPot/ElastoPlasticQPot.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;

// abbreviate types(s)
typedef ElastoPlasticQPot::ArrD ArrD;

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
typedef ElastoPlasticQPot::Cartesian2d::T2s T2s;

// -------------------------------------------------------------------------------------------------

sm.def("epsm", py::overload_cast<const T2s & >(&SM::epsm), "Mean strain"        , py::arg("Eps"));
sm.def("epsd", py::overload_cast<const T2s & >(&SM::epsd), "Eq. strain deviator", py::arg("Eps"));
sm.def("epsm", py::overload_cast<const ArrD &>(&SM::epsm), "Mean strain"        , py::arg("Eps"));
sm.def("epsd", py::overload_cast<const ArrD &>(&SM::epsd), "Eq. strain deviator", py::arg("Eps"));
sm.def("sigm", py::overload_cast<const T2s & >(&SM::sigm), "Mean stress"        , py::arg("Sig"));
sm.def("sigd", py::overload_cast<const T2s & >(&SM::sigd), "Eq. stress deviator", py::arg("Sig"));
sm.def("sigm", py::overload_cast<const ArrD &>(&SM::sigm), "Mean stress"        , py::arg("Sig"));
sm.def("sigd", py::overload_cast<const ArrD &>(&SM::sigd), "Eq. stress deviator", py::arg("Sig"));

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
  .def("stress", &SM::Elastic::stress)
  .def("energy", &SM::Elastic::energy)
  .def("epsy"  , &SM::Elastic::epsy  )
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Elastic::epsp, py::const_))
  .def("epsp"  , py::overload_cast<double     >(&SM::Elastic::epsp, py::const_))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Elastic::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&SM::Elastic::find, py::const_))
  // print to screen
  .def("__repr__", [](const SM::Elastic &a){
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
  .def("stress", &SM::Cusp::stress)
  .def("energy", &SM::Cusp::energy)
  .def("epsy"  , &SM::Cusp::epsy  )
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Cusp::epsp, py::const_))
  .def("epsp"  , py::overload_cast<double     >(&SM::Cusp::epsp, py::const_))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Cusp::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&SM::Cusp::find, py::const_))
  // print to screen
  .def("__repr__", [](const SM::Cusp &a){
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
  .def("stress", &SM::Smooth::stress)
  .def("energy", &SM::Smooth::energy)
  .def("epsy"  , &SM::Smooth::epsy  )
  .def("epsp"  , py::overload_cast<const T2s &>(&SM::Smooth::epsp, py::const_))
  .def("epsp"  , py::overload_cast<double     >(&SM::Smooth::epsp, py::const_))
  .def("find"  , py::overload_cast<const T2s &>(&SM::Smooth::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&SM::Smooth::find, py::const_))
  // print to screen
  .def("__repr__", [](const SM::Smooth &a){
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
  .def("addElastic", &SM::Matrix::addElastic,py::arg("index"),py::arg("K"),py::arg("G"))
  .def("addCusp"   , &SM::Matrix::addCusp   ,py::arg("index"),py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("addSmooth" , &SM::Matrix::addSmooth ,py::arg("index"),py::arg("K"),py::arg("G"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("type"      , &SM::Matrix::type)
  .def("stress"    , &SM::Matrix::stress)
  .def("energy"    , &SM::Matrix::energy)
  .def("epsy"      , &SM::Matrix::epsy)
  .def("epsp"      , &SM::Matrix::epsp)
  .def("find"      , &SM::Matrix::find)
  // print to screen
  .def("__repr__", [](const SM::Matrix &a){
    return "<ElastoPlasticQPot.Cartesian2d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}


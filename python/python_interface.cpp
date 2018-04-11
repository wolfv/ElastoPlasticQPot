/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "../src/ElastoPlasticQPot/ElastoPlasticQPot.h"

namespace py = pybind11;

// =================================================================================================

typedef ElastoPlasticQPot::ArrD ArrD;

// =================================================================================================

PYBIND11_MODULE(ElastoPlasticQPot, m) {

// =================================================================================================

// set doc-string
m.doc() = "Elasto-plastic material model";

// ================================ ElastoPlasticQPot::Cartesian2d =================================

py::module m2d = m.def_submodule("Cartesian2d", "2d Cartesian coordinates");

{

typedef ElastoPlasticQPot::Cartesian2d::T2s T2s;

namespace M = ElastoPlasticQPot::Cartesian2d;

// -------------------------------------------------------------------------------------------------

m2d.def("epsm",
  py::overload_cast<const T2s &>(&M::epsm),
  "Compute mean strain",
  py::arg("Eps")
);

// -------------------------------------------------------------------------------------------------

m2d.def("epsd",
  py::overload_cast<const T2s &>(&M::epsd),
  "Compute equivalent strain deviator",
  py::arg("Eps")
);

// -------------------------------------------------------------------------------------------------

m2d.def("epsm",
  py::overload_cast<const ArrD &>(&M::epsm),
  "Compute mean strain",
  py::arg("Eps")
);

// -------------------------------------------------------------------------------------------------

m2d.def("epsd",
  py::overload_cast<const ArrD &>(&M::epsd),
  "Compute equivalent strain deviator",
  py::arg("Eps")
);

// -------------------------------------------------------------------------------------------------

py::class_<M::Elastic>(m2d, "Elastic")
  // constructor
  .def(
    py::init<double,double>(),
    "Elastic material",
    py::arg("K"),
    py::arg("G")
  )
  // methods
  .def("stress", &M::Elastic::stress)
  .def("energy", &M::Elastic::energy)
  .def("eps_y" , &M::Elastic::eps_y )
  .def("find"  , py::overload_cast<const T2s &>(&M::Elastic::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&M::Elastic::find, py::const_))
  // print to screen
  .def("__repr__", [](const M::Elastic &a){
    return "<ElastoPlasticQPot.Cartesian2d.Elastic>"; });

// -------------------------------------------------------------------------------------------------

py::class_<M::Cusp>(m2d, "Cusp")
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
  .def("stress", &M::Cusp::stress)
  .def("energy", &M::Cusp::energy)
  .def("eps_y" , &M::Cusp::eps_y )
  .def("find"  , py::overload_cast<const T2s &>(&M::Cusp::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&M::Cusp::find, py::const_))
  // print to screen
  .def("__repr__", [](const M::Cusp &a){
    return "<ElastoPlasticQPot.Cartesian2d.Cusp>"; });

// -------------------------------------------------------------------------------------------------

py::class_<M::Smooth>(m2d, "Smooth")
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
  .def("stress", &M::Smooth::stress)
  .def("energy", &M::Smooth::energy)
  .def("eps_y" , &M::Smooth::eps_y )
  .def("find"  , py::overload_cast<const T2s &>(&M::Smooth::find, py::const_))
  .def("find"  , py::overload_cast<double     >(&M::Smooth::find, py::const_))
  // print to screen
  .def("__repr__", [](const M::Smooth &a){
    return "<ElastoPlasticQPot.Cartesian2d.Smooth>"; });

// -------------------------------------------------------------------------------------------------

py::enum_<M::Type::Value>(m2d, "Type")
    .value("Unset"       , M::Type::Unset)
    .value("Elastic"     , M::Type::Elastic)
    .value("Cusp"        , M::Type::Cusp)
    .value("Smooth"      , M::Type::Smooth)
    .value("PlanarCusp"  , M::Type::PlanarCusp)
    .value("PlanarSmooth", M::Type::PlanarSmooth)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<M::Matrix>(m2d, "Matrix")
  // constructor
  .def(
    py::init<const std::vector<size_t>&>(),
    "Matrix of materials",
    py::arg("shape")
  )
  // methods
  .def("addElastic", &M::Matrix::addElastic)
  .def("addCusp"   , &M::Matrix::addCusp)
  .def("addSmooth" , &M::Matrix::addSmooth)
  .def("type"      , &M::Matrix::type)
  .def("stress"    , &M::Matrix::stress)
  .def("energy"    , &M::Matrix::energy)
  .def("find"      , &M::Matrix::find)
  .def("eps_y"     , &M::Matrix::eps_y)
  // print to screen
  .def("__repr__", [](const M::Matrix &a){
    return "<ElastoPlasticQPot.Cartesian2d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}


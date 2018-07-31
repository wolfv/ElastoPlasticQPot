
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include "../include/xElastoPlasticQPot/ElastoPlasticQPot.h"

namespace GM = xElastoPlasticQPot::Cartesian2d;

// =================================================================================================

TEST_CASE("ElastoPlasticQPot::Cartesian2d", "Cartesian2d.h")
{

// =================================================================================================

SECTION( "Elastic" )
{
  // material model
  // - parameters
  double K = 12.3;
  double G = 45.6;
  // - model
  GM::Elastic mat(K,G);

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm  );
  EQ( Sig(1,1), K * epsm  );
  EQ( Sig(0,1), G * gamma );
  EQ( Sig(1,0), G * gamma );
  // - plastic strain
  EQ( mat.epsp(Eps), 0 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 0 );
}

// =================================================================================================

SECTION( "Cusp" )
{
  // material model
  // - parameters
  double K = 12.3;
  double G = 45.6;
  // - model
  GM::Cusp mat(K,G,{0.01, 0.03, 0.10});

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm );
  EQ( Sig(1,1), K * epsm );
  EQ( Sig(0,1), 0.       );
  EQ( Sig(1,0), 0.       );
  // - plastic strain
  EQ( mat.epsp(Eps), 0.02 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 1 );
}

// =================================================================================================

SECTION( "Smooth" )
{
  // material model
  // - parameters
  double K = 12.3;
  double G = 45.6;
  // - model
  GM::Smooth mat(K,G,{0.01, 0.03, 0.10});

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm );
  EQ( Sig(1,1), K * epsm );
  EQ( Sig(0,1), 0.       );
  EQ( Sig(1,0), 0.       );
  // - plastic strain
  EQ( mat.epsp(Eps), 0.02 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 1 );
}

// =================================================================================================

SECTION( "Matrix" )
{
  // parameters
  double K = 12.3;
  double G = 45.6;

  // allocate matrix
  GM::Matrix mat({3,2});

  // row 0: elastic
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat.shape());

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(0,k) = 1;

    mat.setElastic(I,K,G);
  }

  // row 1: cups
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat.shape());

    std::vector<double> epsy = {0.01, 0.03, 0.10};

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(1,k) = 1;

    mat.setCusp(I,K,G,epsy);
  }

  // row 2: smooth
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat.shape());

    std::vector<double> epsy = {0.01, 0.03, 0.10};

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(2,k) = 1;

    mat.setCusp(I,K,G,epsy);
  }

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - strain/stress matrices
  xt::xtensor<double,4> eps = xt::empty<double>({3,2,2,2});
  xt::xtensor<double,4> sig;
  xt::xtensor<double,2> epsp;
  // - set strain
  for ( size_t e = 0 ; e < 3 ; ++e ) {
    for ( size_t k = 0 ; k < 2 ; ++k ) {
      auto eps_i = xt::view(eps, e, k, xt::all(), xt::all());
      eps_i = Eps;
    }
  }
  // - stress & plastic strain
  sig  = mat.Sig (eps);
  epsp = mat.epsp(eps);
  // - analytical solution
  EQ( sig(0,0,0,0), K * epsm ); EQ( sig(0,1,0,0), K * epsm );
  EQ( sig(0,0,1,1), K * epsm ); EQ( sig(0,1,1,1), K * epsm );
  EQ( sig(0,0,0,1), G * gamma); EQ( sig(0,1,0,1), G * gamma);
  EQ( sig(0,0,1,0), G * gamma); EQ( sig(0,1,1,0), G * gamma);
  EQ( sig(1,0,0,0), K * epsm ); EQ( sig(1,1,0,0), K * epsm );
  EQ( sig(1,0,1,1), K * epsm ); EQ( sig(1,1,1,1), K * epsm );
  EQ( sig(1,0,0,1), 0.       ); EQ( sig(1,1,0,1), 0.       );
  EQ( sig(1,0,1,0), 0.       ); EQ( sig(1,1,1,0), 0.       );
  EQ( sig(2,0,0,0), K * epsm ); EQ( sig(2,1,0,0), K * epsm );
  EQ( sig(2,0,1,1), K * epsm ); EQ( sig(2,1,1,1), K * epsm );
  EQ( sig(2,0,0,1), 0.       ); EQ( sig(2,1,0,1), 0.       );
  EQ( sig(2,0,1,0), 0.       ); EQ( sig(2,1,1,0), 0.       );
  // - plastic strain
  EQ( epsp(0,0), 0    ); EQ( epsp(0,1), 0    );
  EQ( epsp(1,0), gamma); EQ( epsp(1,1), gamma);
  EQ( epsp(2,0), gamma); EQ( epsp(2,1), gamma);
}

// =================================================================================================

}

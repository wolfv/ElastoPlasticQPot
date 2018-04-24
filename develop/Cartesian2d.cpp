
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include <ElastoPlasticQPot/ElastoPlasticQPot.h>

namespace GMat = ElastoPlasticQPot::Cartesian2d;

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
  GMat::Elastic mat(K,G);

  // allocate tensors
  GMat::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.stress(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm  );
  EQ( Sig(1,1), K * epsm  );
  EQ( Sig(0,1), G * gamma );
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
  GMat::Cusp mat(K,G,{0.01, 0.03, 0.10});

  // allocate tensors
  GMat::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.stress(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm );
  EQ( Sig(1,1), K * epsm );
  EQ( Sig(0,1), G * 0.   );
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
  GMat::Smooth mat(K,G,{0.01, 0.03, 0.10});

  // allocate tensors
  GMat::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.stress(Eps);
  // - analytical solution
  EQ( Sig(0,0), K * epsm );
  EQ( Sig(1,1), K * epsm );
  EQ( Sig(0,1), G * 0.   );
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
  GMat::Matrix mat({3,2});

  // row 0: elastic
  {
    cppmat::matrix2<size_t> index(2,2);
    cppmat::vector <double> k    (2  );
    cppmat::vector <double> g    (2  );

    k.setConstant(K);
    g.setConstant(G);

    index(0,0) = 0;  index(0,1) = 0;
    index(1,0) = 0;  index(1,1) = 1;

    mat.addElastic(index,k,g);
  }

  // row 1: cups
  {
    cppmat::matrix2<size_t> index(2,2);
    cppmat::vector <double> k    (2  );
    cppmat::vector <double> g    (2  );
    cppmat::matrix2<double> epsy (2,3);

    k.setConstant(K);
    g.setConstant(G);

    epsy(0,0) = 0.01; epsy(0,1) = 0.03; epsy(0,2) = 0.10;
    epsy(1,0) = 0.01; epsy(1,1) = 0.03; epsy(1,2) = 0.10;

    index(0,0) = 1;  index(0,1) = 0;
    index(1,0) = 1;  index(1,1) = 1;

    mat.addCusp(index,k,g,epsy);
  }

  // row 2: smooth
  {
    cppmat::matrix2<size_t> index(2,2);
    cppmat::vector <double> k    (2  );
    cppmat::vector <double> g    (2  );
    cppmat::matrix2<double> epsy (2,3);

    k.setConstant(K);
    g.setConstant(G);

    epsy(0,0) = 0.01; epsy(0,1) = 0.03; epsy(0,2) = 0.10;
    epsy(1,0) = 0.01; epsy(1,1) = 0.03; epsy(1,2) = 0.10;

    index(0,0) = 2;  index(0,1) = 0;
    index(1,0) = 2;  index(1,1) = 1;

    mat.addSmooth(index,k,g,epsy);
  }

  // allocate tensors
  GMat::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = epsm; Eps(0,1) = gamma;
  // - strain/stress matrices
  cppmat::matrix<double> eps({3,2,Eps.size()}), sig, epsp;
  // - set strain
  for ( size_t e = 0 ; e < 3 ; ++e )
    for ( size_t k = 0 ; k < 2 ; ++k )
      std::copy(Eps.begin(), Eps.end(), eps.item(e,k));
  // - stress & plastic strain
  sig  = mat.stress(eps);
  epsp = mat.epsp  (eps);
  // - analytical solution
  EQ( sig(0,0,0), K * epsm ); EQ( sig(0,1,0), K * epsm );
  EQ( sig(0,0,2), K * epsm ); EQ( sig(0,1,2), K * epsm );
  EQ( sig(0,0,1), G * gamma); EQ( sig(0,1,1), G * gamma);
  EQ( sig(1,0,0), K * epsm ); EQ( sig(1,1,0), K * epsm );
  EQ( sig(1,0,2), K * epsm ); EQ( sig(1,1,2), K * epsm );
  EQ( sig(1,0,1), G * 0.   ); EQ( sig(1,1,1), G * 0.   );
  EQ( sig(2,0,0), K * epsm ); EQ( sig(2,1,0), K * epsm );
  EQ( sig(2,0,2), K * epsm ); EQ( sig(2,1,2), K * epsm );
  EQ( sig(2,0,1), G * 0.   ); EQ( sig(2,1,1), G * 0.   );
  // - plastic strain
  EQ( epsp(0,0), 0    ); EQ( epsp(0,1), 0    );
  EQ( epsp(1,0), gamma); EQ( epsp(1,1), gamma);
  EQ( epsp(2,0), gamma); EQ( epsp(2,1), gamma);
}

// =================================================================================================

}
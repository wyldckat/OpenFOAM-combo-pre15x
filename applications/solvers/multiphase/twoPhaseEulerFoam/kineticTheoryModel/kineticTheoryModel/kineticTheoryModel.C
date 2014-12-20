/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Everything in this module is taken from the Ph.D. thesis of
    Berend van Wachem:

    'Derivation, Implementation and Validation 
                    of 
          Computer Simulation Models
         for Gas-Solid Fluidized Beds'
    
\*---------------------------------------------------------------------------*/

#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phasea,
    const Foam::volVectorField& Ub,
    const Foam::volScalarField& alpha,
    const Foam::dragModel& draga
)
:
    phasea_(phasea),
    Ua_(phasea.U()),
    Ub_(Ub),
    alpha_(alpha),
    phia_(phasea.phi()),
    draga_(draga),

    rhoa_(phasea.rho()),
    da_(phasea.d()),
    nua_(phasea.nu()),

    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            Ua_.time().constant(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    equilibrium_(kineticTheoryProperties_.lookup("equilibrium")),
    viscosityModel_
    (
        viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    alphaMinFriction_(kineticTheoryProperties_.lookup("alphaMinFriction")),
    Fr_(kineticTheoryProperties_.lookup("Fr")),
    eta_(kineticTheoryProperties_.lookup("eta")),
    p_(kineticTheoryProperties_.lookup("p")),
    phi_(dimensionedScalar(kineticTheoryProperties_.lookup("phi"))*M_PI/180.0),
    Theta_
    (
        IOobject
        (
            "Theta",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh()
    ),
    mua_
    (
        IOobject
        (
            "mua",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModel::solve()
{
    if (!kineticTheory_)
    {
        return;
    }

    word scheme("div(phi, Theta)");

    volScalarField alpha = alpha_;
    alpha.max(1.0e-5);
    const scalar piSqrt = pow(M_PI, 0.5);
        
    surfaceScalarField phi = 1.5*rhoa_*phia_*fvc::interpolate(alpha_);
        
    volTensorField dU = fvc::grad(Ua_);
    volTensorField dUT = dU.T();
    volTensorField D = 0.5*(dU + dUT);
        
    dimensionedScalar time1("1", dimensionSet(0,0,1,0,0,0,0), 1.0);
        
    // NB, drag = K*alpha*beta,
    // (the alpha and beta has been extracted from the drag function for
    // numerical reasons)
    volScalarField Ur = mag(Ua_ - Ub_);
    volScalarField betaPrim = alpha_*(scalar(1) - alpha_)*draga_.K(Ur);
        
    // Gidaspow, p.53 Table 3.4
    volScalarField gs0 = g0();
        
    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff =
        rhoa_*alpha_*(scalar(1) + 2.0*(1.0 + e_)*gs0*alpha_);
        
    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0,2,-2,0,0,0,0),
        1.0e-6
    );

    dimensionedScalar TsmallSqrt = pow(Tsmall, 0.5);
    volScalarField ThetaSqrt = pow(Theta_, 0.5);

    // 'thermal' conductivity (Table 3.3, p. 49)
    volScalarField kappa = 
        conductivityModel_->kappa(alpha, Theta_, gs0, rhoa_, da_, e_);

    // particle viscosity (Table 3.2, p.47)
    mua_ = viscosityModel_->mua(alpha, Theta_, gs0, rhoa_, da_, e_);

    // dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff =
        12.0*(1.0 - e_*e_)*pow(alpha, 2.0)*rhoa_*gs0*(1.0/da_)*ThetaSqrt/piSqrt;
        
    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1 = 3.0*betaPrim;
    volScalarField J2 = 
        0.25*pow(betaPrim, 2.0)*da_*pow(Ur, 2.0)
       /(alpha*rhoa_*piSqrt*(ThetaSqrt + TsmallSqrt));
        
    // bulk viscosity  p. 45 (Lun et al. 1984).
    volScalarField lambda = 
        (4.0/3.0)*pow(alpha_, 2.0)*rhoa_*da_*gs0*(1.0+e_)*ThetaSqrt/piSqrt;
        
    // frictional normal stress Eq. 3.29, p.52 (Johnson and Jackson, 1987)
    // NB! pf = 0.0 for alpha < alphaMinFriction
    // can become VERY large close to max packing hence quite a large
    // lower limit in the denominator
    volScalarField pf = 
        Fr_*pow(max(alpha-alphaMinFriction_, scalar(0)), eta_)
       /pow(max(alphaMax_-alpha, scalar(5.0e-2)), p_);
        
    // add frictional stress for alpha > alphaMin
    PsCoeff += pf/(Theta_+Tsmall);

    PsCoeff.min(1.0e+10);
    PsCoeff.max(-1.0e+10);

    // update particle pressure
    pa_ = PsCoeff*Theta_;
        
    // frictional shear stress, Eq. 3.30, p. 52
    // (with a correcttion of dimension)
    volScalarField muf = 0.5*pf*sin(phi_)*time1;

    // add frictional stress for alpha > alphaMin
    mua_ += muf;
        
    // stress tensor, Definitions, Table 3.1, p. 43
    volTensorField tau = 2.0*mua_*D + (lambda - (2.0/3.0)*mua_)*tr(D)*I;
        
    if (!equilibrium_)
    {
        // construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20
        // no grad infront of Ps
        // wrong sign infront of laplacian
        fvScalarMatrix ThetaEqn
        (
            fvm::ddt(1.5*alpha*rhoa_, Theta_)
          + fvm::div(phi, Theta_, scheme)
         ==
            fvm::SuSp(-((PsCoeff*I) && dU), Theta_)
          + (tau && dU)
          + fvm::laplacian(kappa, Theta_, "laplacian(kappa,Theta)")
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/(Theta_+Tsmall), Theta_)
        );
        
        ThetaEqn.relax();
        ThetaEqn.solve();
    }
    else
    {
        // equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1 = 2.0*(1.0 + e_)*rhoa_*gs0;
        volScalarField K3 = 0.5*da_*rhoa_*
            (
                (piSqrt/(3.0*(3.0-e_)))
               *(scalar(1) + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*gs0)
              + 1.6*alpha*gs0*(1.0 + e_)/piSqrt
            );
            
        volScalarField K2 = 
            4.0*da_*rhoa_*(1.0 + e_)*alpha*gs0/(3.0*piSqrt) - 2.0*K3/3.0;

        volScalarField K4 = 12.0*(1.0 - e_*e_)*rhoa_*gs0/(da_*piSqrt);
            
        volScalarField trD = tr(D);
        volTensorField D2 = D & D;
        volScalarField tr2D = trD*trD;
        volScalarField trD2 = tr(D2);
        
        volScalarField t1 = K1*alpha + rhoa_;
        volScalarField l1 = -t1*trD;
        volScalarField l2 = pow(t1, 2)*tr2D;
        volScalarField l3 = 4.0*K4*alpha*(2.0*K3*tr2D + K2*tr2D);

        Theta_ = 
            pow((l1 + sqrt(l2 + l3))/(2.0*(alpha + scalar(1.0e-4))*K4), 2.0);
    }

    Theta_.max(1.0e-15);
    Theta_.min(1.0e+3);

    // update properties
    mua_ = viscosityModel_->mua(alpha, Theta_, gs0, rhoa_, da_, e_) + muf;
    pa_ = PsCoeff*Theta_ + pf;

    Info<< "kinTheory: max(Theta) = " << max(Theta_).value() << endl;

    volScalarField ktn = mua_/rhoa_;
        
    Info<< "kinTheory: min(nua) = " << min(ktn).value()
        << ", max(nua) = " << max(ktn).value() << endl;

    Info<< "kinTheory: min(pa) = " << min(pa_).value()
        << ", max(pa) = " << max(pa_).value() << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModel::g0()
{
    // NN. should be made runtime selectable like viscosity and conductivity
    // models, but for now its Gidaspow...
    volScalarField alpha = alpha_;
    alpha.min(alphaMax_-1.0e-2);
    alpha.max(1.0e-10);

    // this term is quite problematic...
    // we should handle alpha > alphaMax_ better using some linear
    // extrapolation instead of just clipping alpha
    return (0.6/(scalar(1) - pow(alpha/alphaMax_, 1.0/3.0)));
}


// ************************************************************************* //

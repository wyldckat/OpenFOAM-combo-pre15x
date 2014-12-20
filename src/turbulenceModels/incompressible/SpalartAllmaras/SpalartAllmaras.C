/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    
\*---------------------------------------------------------------------------*/

#include "SpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(turbulenceModel, SpalartAllmaras, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::fv1() const
{
    volScalarField chi3 = pow3(nuTilda_/nu());

    return chi3/(chi3 + pow3(Cv1));
}


tmp<volScalarField> SpalartAllmaras::fv2() const
{
    volScalarField chi = nuTilda_/nu();
    return 1.0 - chi/(1.0 + chi*fv1());
}


tmp<volScalarField> SpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r = nuTilda_/(Stilda*sqr(kappa_*d_));
    volScalarField g = r + Cw2*(pow(r, 6) - r);

    return g*pow((1.0 + pow(Cw3, 6))/(pow(g, 6) + pow(Cw3, 6)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
SpalartAllmaras::SpalartAllmaras
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(typeName, U, phi, lamTransportModel),

    alphaNut(turbulenceModelCoeffs_.lookup("alphaNut")),

    Cb1(turbulenceModelCoeffs_.lookup("Cb1")),
    Cb2(turbulenceModelCoeffs_.lookup("Cb2")),
    Cw1(Cb1/sqr(kappa_) + alphaNut*(1.0 + Cb2)),
    Cw2(turbulenceModelCoeffs_.lookup("Cw2")),
    Cw3(turbulenceModelCoeffs_.lookup("Cw3")),
    Cv1(turbulenceModelCoeffs_.lookup("Cv1")),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    d_(mesh_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volTensorField> SpalartAllmaras::R() const
{
    return tmp<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k() - nut()*2*symm(fvc::grad(U_)),
            nuTilda_.boundaryField().types()
        )
    );
}


tmp<fvVectorMatrix> SpalartAllmaras::divR(volVectorField& U) const
{
    volScalarField nuEff_ = nuEff();

    return
    (
      - fvm::laplacian(nuEff_, U)
      - fvc::div(nuEff_*dev(fvc::grad(U)().T()))
    );
}


bool SpalartAllmaras::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.lookup("alphaNut") >> alphaNut;

        turbulenceModelCoeffs_.lookup("Cb1") >> Cb1;
        turbulenceModelCoeffs_.lookup("Cb2") >> Cb2;
        Cw1 = Cb1/sqr(kappa_) + alphaNut*(1.0 + Cb2);
        turbulenceModelCoeffs_.lookup("Cw2") >> Cw2;
        turbulenceModelCoeffs_.lookup("Cw3") >> Cw3;
        turbulenceModelCoeffs_.lookup("Cv1") >> Cv1;

        return true;
    }
    else
    {
        return false;
    }
}


void SpalartAllmaras::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    turbulenceModel::correct();

    if (mesh_.moving())
    {
        d_.correct();
    }

    volScalarField Stilda =
        ::sqrt(2.0)*mag(skew(fvc::grad(U_)))
      + nuTilda_*fv2()/sqr(kappa_*d_);

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(), nuTilda_)
      - alphaNut*Cb2*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1*Stilda*nuTilda_ - fvm::Sp(Cw1*fw(Stilda)*nuTilda_/sqr(d_), nuTilda_)
    );

    nuTildaEqn().relax();
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceModels
} // End namespace Foam

// ************************************************************************* //

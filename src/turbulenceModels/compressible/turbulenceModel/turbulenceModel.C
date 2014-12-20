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

#include "turbulenceModel.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulenceModel, 0);
defineRunTimeSelectionTable(turbulenceModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModel::turbulenceModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    IOdictionary
    (
        IOobject
        (
            "turbulenceProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),
    phi_(phi),
    thermophysicalModel_(thermophysicalModel),

    turbulence_(lookup("turbulence")),
    turbulenceModelCoeffs_(subDict(type + "Coeffs")),

    kappa_(subDict("wallFunctionCoeffs").lookup("kappa")),
    E_(subDict("wallFunctionCoeffs").lookup("E")),
    yPlusLam_(yPlusLam(kappa_.value(), E_.value())),

    k0_("k0", dimVelocity*dimVelocity, SMALL),
    epsilon0_("epsilon", k0_.dimensions()/dimTime, SMALL),
    epsilonSmall_("epsilonSmall", epsilon0_.dimensions(), SMALL),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar turbulenceModel::yPlusLam(const scalar kappa, const scalar E)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<scalarField> turbulenceModel::yPlus(const label patchNo) const
{
    const fvPatch& curPatch = mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp();

    if (typeid(curPatch) == typeid(wallFvPatch))
    {
        dimensionedScalar Cmu(turbulenceModelCoeffs_.lookup("Cmu"));

        Yp = pow(Cmu.value(), 0.25)
            *y_[patchNo]
            *sqrt(k()().boundaryField()[patchNo].patchInternalField())
           /(
                mu().boundaryField()[patchNo].patchInternalField()
                *rho_.boundaryField()[patchNo]
            );
    }
    else
    {
        Warning
            << "tmp<scalarField> turbulenceModel::yPlus(const label patchNo) "
            << "const : " << endl
            << "Patch " << patchNo << " is not a wall.  Returning blank field"
            << endl;

        Yp.setSize(0);
    }

    return tYp;
}


void turbulenceModel::correct()
{
    if (mesh_.moving())
    {
        y_.correct();
    }
}


bool turbulenceModel::read()
{
    if (regIOobject::read())
    {
        lookup("turbulence") >> turbulence_;
        turbulenceModelCoeffs_ = subDict(type() + "Coeffs");

        subDict("wallFunctionCoeffs").lookup("kappa") >> kappa_;
        subDict("wallFunctionCoeffs").lookup("E") >> E_;
        yPlusLam_ = yPlusLam(kappa_.value(), E_.value());

        if (found("k0"))
        {
            lookup("k0") >> k0_;
        }

        if (found("epsilon0"))
        {
            lookup("epsilon0") >> epsilon0_;
        }

        if (found("epsilonSmall"))
        {
            lookup("epsilonSmall") >> epsilonSmall_;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

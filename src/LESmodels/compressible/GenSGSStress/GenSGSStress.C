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

Description

\*---------------------------------------------------------------------------*/

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESmodels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
GenSGSStress::GenSGSStress
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel
)
:
    LESmodel
    (
        word("GenSGSStress"),
        rho,
        U,
        phi,
        thermoPhysicalModel
    ),

    ce_(LESmodelProperties().lookup("ce")),

    B_
    (
        IOobject
        (
            "B",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    muSgs_
    (
        IOobject
        (
            "muSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> GenSGSStress::divRhoB
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho()*B_ + 0.05*muSgs_*fvc::grad(U))
      + fvc::laplacian(0.95*muSgs_, U, "laplacian(muEff,U)")
      - fvm::laplacian(muEff(), U)
    );
}


void GenSGSStress::correct(const tmp<volTensorField>& gradU)
{
    LESmodel::correct(gradU);
}


bool GenSGSStress::read()
{
    if (LESmodel::read())
    {
        LESmodelProperties().lookup("ce") >> ce_;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

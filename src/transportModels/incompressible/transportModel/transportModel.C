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

Class
    transportModel

\*---------------------------------------------------------------------------*/

#include "transportModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(transportModel, 0);
defineRunTimeSelectionTable(transportModel, dictionary);


//- Lookup and return the phase transport properties dictionary
const dictionary& transportModel::lookupPhaseTransportProperties
(
    const dictionary& transportProperties,
    const word& phaseName
)
{
    if (phaseName.size())
    {
        return transportProperties.subDict(phaseName);
    }
    else
    {
        return transportProperties;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transportModel::transportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    U_(U),
    phi_(phi),
    phaseName_(phaseName),
    phaseTransportProperties_
    (
        lookupPhaseTransportProperties(*this, phaseName)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> transportModel::strainRate() const
{
    return mag(fvc::grad(U_));
}


bool transportModel::read()
{
    if (regIOobject::read())
    {
        phaseTransportProperties_ =
            lookupPhaseTransportProperties(*this, phaseName_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

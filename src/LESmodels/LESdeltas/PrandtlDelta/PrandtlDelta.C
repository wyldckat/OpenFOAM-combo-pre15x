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

#include "PrandtlDelta.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PrandtlDelta, 0);
addToRunTimeSelectionTable(LESdelta, PrandtlDelta, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void PrandtlDelta::calcDelta()
{
    delta_ = min
    (
        (const volScalarField&)(geometricDelta_()),
        (kappa_/Cdelta_)*(const volScalarField&)wallDist(mesh_)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PrandtlDelta::PrandtlDelta(const fvMesh& mesh, const dictionary& dd)
:
    LESdelta(mesh),
    geometricDelta_(LESdelta::New(mesh, dd.subDict(type() + "Coeffs"))),
    kappa_(dd.lookup("kappa")),
    Cdelta_(dd.subDict(type() + "Coeffs").lookup("Cdelta"))
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void PrandtlDelta::read(const dictionary& d)
{
    const dictionary& dd(d.subDict(type() + "Coeffs"));

    geometricDelta_().read(dd);
    d.lookup("kappa") >> kappa_;
    dd.lookup("Cdelta") >> Cdelta_;
    calcDelta();
}


void PrandtlDelta::correct()
{
    geometricDelta_().correct();

    if (mesh_.moving())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

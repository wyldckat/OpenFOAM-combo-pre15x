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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Toroidal coordinate system

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "toroidalCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(toroidalCS, 0);

addToRunTimeSelectionTable(coordinateSystem, toroidalCS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
toroidalCS::toroidalCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& direction,
    const scalar radius
)
:
    cartesianCS(name, origin, axis, direction),
    radius_(radius)
{}


// Construct from origin and a coordinate rotation
toroidalCS::toroidalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr,
    const scalar radius
)
:
    cartesianCS(name, origin, cr),
    radius_(radius)
{}


toroidalCS::toroidalCS
(
    const word& name,
    const dictionary& dict
)
:
    cartesianCS(name, dict),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Convert from local coordinate system to the global Cartesian system
vector toroidalCS::toGlobal(const vector& localV) const
{
    // Notation: r = localV.x()
    scalar theta = localV.y()*mathematicalConstant::pi/180.0;
    scalar phi = localV.z()*mathematicalConstant::pi/180.0;

    scalar rprime = radius_ + localV.x()*sin(phi);

    if ((localV.x()*sin(phi)) > (radius_))
    {
        FatalErrorIn("toroidalCS::toGlobal(vector) const")
            << "Badly defined toroidal coordinates"
            << abort(FatalError);
    }

    return cartesianCS::toGlobal
    (
        vector(rprime*cos(theta), rprime*sin(theta), localV.x()*cos(phi))
    );
}


tmp<vectorField> toroidalCS::toGlobal
(
    const vectorField& localV
) const
{
    const scalarField r = localV.component(vector::X);

    const scalarField theta =
        localV.component(vector::Y)*mathematicalConstant::pi/180.0;

    const scalarField phi =
        localV.component(vector::Z)*mathematicalConstant::pi/180.0;

    const scalarField rprime = radius_ + r*sin(phi);

    vectorField lc(localV.size());
    lc.replace(vector::X, rprime*cos(theta));
    lc.replace(vector::Y, rprime*sin(theta));
    lc.replace(vector::Z, r*cos(phi));

    return cartesianCS::toGlobal(lc);
}


// Convert from global Cartesian coordinate system to the local system
vector toroidalCS::toLocal(const vector& globalV) const
{
    notImplemented("vector toroidalCS::toLocal(const vector& globalV) const");

    return vector::zero;
}


tmp<vectorField> toroidalCS::toLocal
(
    const vectorField& globalV
) const
{
    notImplemented
    (
        "tmp<vectorField> toroidalCS::toLocal(const vectorField& globalV) const"
    );

    return tmp<vectorField>(&vectorField::null());
}


void toroidalCS::write(Ostream& os) const
{
    coordinateSystem::write(os);
    os << "radius: " << radius() << endl;
}


void toroidalCS::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT
        << "    origin " << origin() << token::END_STATEMENT << nl
        << "    axis " << axis() << token::END_STATEMENT << nl
        << "    direction " << direction() << token::END_STATEMENT << nl
        << "    radius " << radius() << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

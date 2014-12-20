/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "toroidalCS.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(toroidalCS, 0);

    addToRunTimeSelectionTable(coordinateSystem, toroidalCS, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::toroidalCS::toroidalCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& direction,
    const scalar radius
)
:
    coordinateSystem(name, origin, axis, direction),
    radius_(radius)
{}


Foam::toroidalCS::toroidalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr,
    const scalar radius
)
:
    coordinateSystem(name, origin, cr),
    radius_(radius)
{}


Foam::toroidalCS::toroidalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::toroidalCS::toGlobal(const vector& localV) const
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

    return coordinateSystem::toGlobal
    (
        vector(rprime*cos(theta), rprime*sin(theta), localV.x()*cos(phi))
    );
}


Foam::tmp<Foam::vectorField> Foam::toroidalCS::toGlobal
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

    return coordinateSystem::toGlobal(lc);
}


Foam::vector Foam::toroidalCS::toLocal(const vector& globalV) const
{
    notImplemented("vector toroidalCS::toLocal(const vector& globalV) const");

    return vector::zero;
}


Foam::tmp<Foam::vectorField> Foam::toroidalCS::toLocal
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


void Foam::toroidalCS::write(Ostream& os) const
{
    coordinateSystem::write(os);
    os << "radius: " << radius() << endl;
}


void Foam::toroidalCS::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
	os  << indent << name() << nl
	    << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    coordinateSystem::writeDict(os, false);
    os.writeKeyword("radius") << radius() << token::END_STATEMENT << nl;
   
    if (subDict)
    {
	os << decrIndent << indent << token::END_BLOCK << endl;
    }    
}


// ************************************************************************* //

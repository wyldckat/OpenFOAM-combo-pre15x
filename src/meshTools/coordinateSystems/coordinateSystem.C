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

#include "IOstream.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);

    defineRunTimeSelectionTable(coordinateSystem, origAxisDir);
    defineRunTimeSelectionTable(coordinateSystem, origRotation);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
}


const Foam::scalar Foam::coordinateSystem::nonOrthogonalError = 1.0e-8;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordinateSystem::calcTransformations()
{
    vector e1 = direction();
    vector e3 = axis()/mag(axis());

    // relax pickiness about nonorthogonality, but leave axis() alone
    e1 = e1 - (e1 & e3)*e3;
    e1 = e1/mag(e1);

    vector e2 = e3 ^ e1;

    if (mag(e1&e3)/(mag(e1)*mag(e3)) >= nonOrthogonalError)
    {
        FatalErrorIn("void coordinateSystem::calcTransformations()")
            << "coordinate system nonorthogonality " << endl
            << "mag(axis & direction) = " << mag(e1 & e3)
            << abort(FatalError);
    }

    // the global -> local transformation
    Rtrans_.x() = e1;
    Rtrans_.y() = e2;
    Rtrans_.z() = e3;
    
    // the local -> global transformation
    R_ = Rtrans_.T();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& dir
)
:
    name_(name),
    origin_(origin),
    axis_(axis),
    dir_(dir),
    R_(tensor::zero),
    Rtrans_(tensor::zero)
{
    calcTransformations();
}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    origin_(origin),
    axis_(cr.axis()),
    dir_(cr.direction()),
    R_(cr.R()),
    Rtrans_(cr.R().T())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    origin_(vector::zero),
    axis_(vector::zero),
    dir_(vector::zero),
    R_(tensor::zero),
    Rtrans_(tensor::zero)
{
    // unspecified origin is (0 0 0)
    if (dict.found("origin"))
    {
        dict.lookup("origin") >> origin_;
    }

    // specify via coordinateRotation or via axis/direction
    if (dict.found("coordinateRotation"))
    {
        autoPtr<coordinateRotation> cr =
            coordinateRotation::New(dict.subDict("coordinateRotation"));

        R_      = cr().R();
        Rtrans_ = R_.T();

        axis_   = cr().axis();
        dir_    = cr().direction();
    }
    else
    {
        dict.lookup("axis") >> axis_;
        dict.lookup("direction") >> dir_;
        calcTransformations();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystem::~coordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::coordinateSystem::toGlobal(const vector& localV) const
{
    return (R_ & localV) + origin();
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::toGlobal
(
    const vectorField& localV
) const
{
    return (R_ & localV) + origin();
}


Foam::vector Foam::coordinateSystem::toLocal(const vector& globalV) const
{
    return (Rtrans_ & (globalV - origin()));
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::toLocal
(
    const vectorField& globalV
) const
{
    return (Rtrans_ & (globalV - origin()));
}


void Foam::coordinateSystem::write(Ostream& os) const
{
    os  << nl << type()
        << " origin: " << origin()
        << " axis: " << axis() << " direction: " << direction() << nl;
}


void Foam::coordinateSystem::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << name() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    os.writeKeyword("type")      << type()      << token::END_STATEMENT << nl;
    os.writeKeyword("origin")    << origin()    << token::END_STATEMENT << nl;
    os.writeKeyword("axis")      << axis()      << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << direction() << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const coordinateSystem& p");
    return os;
}


// ************************************************************************* //

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
    A subset of mesh points.

\*---------------------------------------------------------------------------*/

#include "pointZone.H"
#include "addToRunTimeSelectionTable.H"
#include "pointZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointZone, 0);
    defineRunTimeSelectionTable(pointZone, dictionary);
    addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::pointZone::pointLookupMap() const
{
    if (!pointLookupMapPtr_)
    {
        calcPointLookupMap();
    }

    return *pointLookupMapPtr_;
}


void Foam::pointZone::calcPointLookupMap() const
{
    if (debug)
    {
        Info<< "void pointZone::calcPointLookupMap() const : "
            << "Calculating point lookup map"
            << endl;
    }

    if (pointLookupMapPtr_)
    {
        FatalErrorIn
        (
            "void pointZone::calcPointLookupMap() const"
        )   << "point lookup map already calculated"
            << abort(FatalError);
    }

    const labelList& addr = addressing();

    pointLookupMapPtr_ = new Map<label>(2*addr.size());
    Map<label>& plm = *pointLookupMapPtr_;

    forAll (addr, pointI)
    {
        plm.insert(addr[pointI], pointI);
    }

    if (debug)
    {
        Info<< "void pointZone::calcPointLookupMap() const : "
            << "Finished calculating point lookup map"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointZone::pointZone
(
    const word& name,
    const labelList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    indirectPointList(zm.mesh().allPoints(), addr),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// Construct from dictionary
Foam::pointZone::pointZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointZoneMesh& zm
)
:
    indirectPointList
    (
        zm.mesh().allPoints(),
        dict.lookup("pointLabels")
    ),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// Construct given the original zone and resetting the
// point list and zone mesh information
Foam::pointZone::pointZone
(
    const pointZone& pz,
    const pointZoneMesh& zm,
    const label index,
    const labelList& addr
)
:
    indirectPointList(zm.mesh().allPoints(), addr),
    name_(pz.name()),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointZone::~pointZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointZone::whichPoint(const label globalPointID) const
{
    const Map<label>& plm = pointLookupMap();

    Map<label>::const_iterator plmIter = plm.find(globalPointID);

    if (plmIter == plm.end())
    {
        return -1;
    }
    else
    {
        return plmIter();
    }
}


const Foam::pointZoneMesh& Foam::pointZone::zoneMesh() const
{
    return zoneMesh_;
}


void Foam::pointZone::clearAddressing()
{
    deleteDemandDrivenData(pointLookupMapPtr_);
}


void Foam::pointZone::resetAddressing(const labelList& addr)
{
    clearAddressing();
    indirectPointList::resetAddressing(addr);
}


void Foam::pointZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << addressing();
}


void Foam::pointZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    addressing().writeEntry("pointLabels", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointZone& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const pointZone& p");
    return os;
}


// ************************************************************************* //

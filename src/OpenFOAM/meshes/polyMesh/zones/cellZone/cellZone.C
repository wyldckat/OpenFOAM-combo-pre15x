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
    A subset of mesh cells.

\*---------------------------------------------------------------------------*/

#include "cellZone.H"
#include "addToRunTimeSelectionTable.H"
#include "cellZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "IOstream.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZone, 0);

    defineRunTimeSelectionTable(cellZone, dictionary);
    addToRunTimeSelectionTable(cellZone, cellZone, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::cellZone::cellLookupMap() const
{
    if (!cellLookupMapPtr_)
    {
        calcCellLookupMap();
    }

    return *cellLookupMapPtr_;
}


void Foam::cellZone::calcCellLookupMap() const
{
    if (debug)
    {
        Info<< "void cellZone::calcCellLookupMap() const : "
            << "Calculating cell lookup map"
            << endl;
    }

    if (cellLookupMapPtr_)
    {
        FatalErrorIn
        (
            "void cellZone::calcCellLookupMap() const"
        )   << "cell lookup map already calculated"
            << abort(FatalError);
    }

    const labelList& addr = addressing();

    cellLookupMapPtr_ = new Map<label>(2*addr.size());
    Map<label>& clm = *cellLookupMapPtr_;

    forAll (addr, cellI)
    {
        clm.insert(addr[cellI], cellI);
    }

    if (debug)
    {
        Info<< "void cellZone::calcCellLookupMap() const : "
            << "Finished calculating cell lookup map"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellZone::cellZone
(
    const word& name,
    const labelList& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    indirectCellList(zm.mesh().allCells(), addr),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    cellLookupMapPtr_(NULL)
{}


// Construct from dictionary
Foam::cellZone::cellZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const cellZoneMesh& zm
)
:
    indirectCellList
    (
        zm.mesh().allCells(),
        dict.lookup("cellLabels")
    ),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    cellLookupMapPtr_(NULL)
{}


// Construct given the original zone and resetting the
//  cell list and zone mesh information
Foam::cellZone::cellZone
(
    const cellZone& cz,
    const cellZoneMesh& zm,
            const label index,
    const labelList& addr
)
:
    indirectCellList(zm.mesh().allCells(), addr),
    name_(cz.name()),
    index_(index),
    zoneMesh_(zm),
    cellLookupMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellZone::~cellZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellZone::whichCell(const label globalCellID) const
{
    const Map<label>& clm = cellLookupMap();

    Map<label>::const_iterator clmIter = clm.find(globalCellID);

    if (clmIter == clm.end())
    {
        return -1;
    }
    else
    {
        return clmIter();
    }
}


const Foam::cellZoneMesh& Foam::cellZone::zoneMesh() const
{
    return zoneMesh_;
}


void Foam::cellZone::clearAddressing()
{
    deleteDemandDrivenData(cellLookupMapPtr_);
}


void Foam::cellZone::resetAddressing(const labelList& addr)
{
    clearAddressing();
    indirectCellList::resetAddressing(addr);
}


void Foam::cellZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << addressing();
}


void Foam::cellZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    addressing().writeEntry("cellLabels", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const cellZone& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const cellZone& p");
    return os;
}


// ************************************************************************* //

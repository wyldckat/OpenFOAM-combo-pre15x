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
    The class implements fast off-diagonal lookup for the lduAddressing class.

\*---------------------------------------------------------------------------*/

#include "fastLduTriLookup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fastLduTriLookup, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fastLduTriLookup::expandRow(const label d) const
{
    // Expand the row d into the next available slot into the main storage

    label endLabel = ldu_.ownerStartAddr()[d + 1];

    const unallocLabelList& neighbour = ldu_.upperAddr();

    for (label i = ldu_.ownerStartAddr()[d]; i < endLabel; i++)
    {
        mainAddr_[nextFree_*ldu_.size() + neighbour[i]] = i;
    }

    diagAddr_[d] = nextFree_;
    if (reverseAddr_[nextFree_] > -1)
    {
        diagAddr_[reverseAddr_[nextFree_]] = -1;
    }

    reverseAddr_[nextFree_] = d;
    ++nextFree_;
    nextFree_ %= band_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fastLduTriLookup::fastLduTriLookup
(
    const lduAddressing& ldu,
    const label band
)
:
    ldu_(ldu),
    band_(band),
    diagAddr_(ldu.size(), -1),
    mainAddr_(band_*ldu.size()),
    reverseAddr_(band_, -1),
    nextFree_(0)
{
    if (debug)
    {
        Info<< "fastLduTriLookup::fastLduTriLookup(const lduAddressing& ldu, "
            << "const label band) : "
            << "Creating fast ldu lookup.  Band = " << band_
            << " Matrix size: " << ldu_.size() << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fastLduTriLookup::~fastLduTriLookup()
{
    if (debug)
    {
        Info<< "fastLduTriLookup::~fastLduTriLookup() : "
            << "Deleting fast ldu lookup." << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fastLduTriLookup::triIndex
(
    const label a,
    const label b
) const
{
    low_ = min(a, b);
    high_ = max(a, b);

    if (diagAddr_[low_] == -1)
    {
        expandRow(low_);
    }

    return mainAddr_[diagAddr_[low_]*ldu_.size() + high_];
}


// ************************************************************************* //

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
    A bounding box defined in terms of the points at it's extremities.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "boundBox.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as the bounding box of the given pointField
boundBox::boundBox(const pointField& points)
:
    min_(vector::zero),
    max_(vector::zero)
{
    if (points.size() == 0)
    {
        Warning
            << "boundBox::boundBox(const pointField& points) : "
            << "cannot find bounding box for zero sized pointField"
            << "returning zero" << endl;

        return;
    }

    min_ = points[0];
    max_ = points[0];

    forAll(points, i)
    {
        min_ = ::Foam::min(min_, points[i]);
        max_ = ::Foam::max(max_, points[i]);
    }

    // Reduce parallel information
    reduce(min_, minOp<point>());
    reduce(max_, maxOp<point>());
}


// Construct from Istream
boundBox::boundBox(Istream& is)
:
    min_(is),
    max_(is)
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const boundBox& b)
{
    return os << b.min() << token::SPACE << b.max();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

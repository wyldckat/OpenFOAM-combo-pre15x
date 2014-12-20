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

\*---------------------------------------------------------------------------*/

#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(wedgePolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, word);
addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void wedgePolyPatch::initTransforms()
{
    const pointField& points = allPoints();

    vector n(static_cast<const faceList&>(*this)[0].normal(points));
    n /= mag(n);

    centreNormal_ =
        vector
        (
            sign(n.x())*(::Foam::max(mag(n.x()), 0.5) - 0.5),
            sign(n.y())*(::Foam::max(mag(n.y()), 0.5) - 0.5),
            sign(n.z())*(::Foam::max(mag(n.z()), 0.5) - 0.5)
        );
    centreNormal_ /= mag(centreNormal_);

    if
    (
        mag(centreNormal_.x() + centreNormal_.y() + centreNormal_.z())
      < (1 - SMALL)
    )
    {
        FatalErrorIn
        (
            "wedgePolyPatch::wedgePolyPatch(const polyPatch&, "
            "const fvBoundaryMesh&)"
        )   << "wedge does not align with a coordinate plane"
            << exit(FatalError);
    }

    axis_ = centreNormal_ ^ n;
    axis_ /= mag(axis_);

    faceT_ = transformationTensor(centreNormal_, n);
    cellT_ = faceT_ & faceT_;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

// Construct from components
wedgePolyPatch::wedgePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm)
{
    initTransforms();
}


// Construct from Istream
wedgePolyPatch::wedgePolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(is, index, bm)
{
    initTransforms();
}


// Construct from dictionary
wedgePolyPatch::wedgePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm)
{
    initTransforms();
}


//- Construct as copy, resetting the boundary mesh
wedgePolyPatch::wedgePolyPatch
(
    const wedgePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{
    initTransforms();
}


//- Construct as copy, resetting the face list and boundary mesh data
wedgePolyPatch::wedgePolyPatch
(
    const wedgePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{
    initTransforms();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

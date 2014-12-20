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

#include "constantPlane.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "interpolation.H"
#include "dictionary.H"
#include "plane.H"
#include "cuttingPlane.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(constantPlane, 0);

addToRunTimeSelectionTable(surface, constantPlane, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::constantPlane::makeTriangles(const cuttingPlane& cut)
{
    // Count triangles
    label nTris = 0;
    forAll(cut.faces(), cutFaceI)
    {
        const face& f = cut.faces()[cutFaceI];

        nTris += f.nTriangles(cut.points());
    }

    // Triangulation uses all points from cut
    points_ = cut.points();

    faces_.setSize(nTris);
    meshCells_.setSize(nTris);

    label triI = 0;
    label oldTriI = 0;

    forAll(cut.faces(), cutFaceI)
    {
        const face& f = cut.faces()[cutFaceI];

        f.triangles(cut.points(), triI, faces_);

        for(label i = oldTriI; i < triI; i++)
        {
            meshCells_[i] = cut.cells()[cutFaceI];
        }

        oldTriI = triI;
    }
}


void Foam::constantPlane::copyFaces(const cuttingPlane& cut)
{
    points_ = cut.points();
    faces_ = cut.faces();
    meshCells_ = cut.cells();
}


void Foam::constantPlane::createGeometry()
{
    cuttingPlane plane(mesh(), planeDesc_);

    if (triangulate_)
    {
        makeTriangles(plane);
    }
    else
    {
        copyFaces(plane);
    }

    Info<< "Created " << name() << " :"
        << "  base:" << planeDesc_.refPoint()
        << "  normal:" << planeDesc_.normal()
        << "  faces:" << faces_.size()
        << "  points:" << points_.size() << endl;
}


template <class T>
Foam::Field<T> Foam::constantPlane::doInterpolate
(
    const word& fieldName,
    const fieldsCache<T>& cache
) const
{
    // One value per face
    Field<T> values(meshCells_.size());

    const GeometricField<T, fvPatchField, volMesh>& field =
        *cache[fieldName];

    forAll(meshCells_, elemI)
    {
        values[elemI] = field[meshCells_[elemI]];
    }

    return values;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::constantPlane::constantPlane
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const plane& planeDesc,
    const bool triangulate
)
:
    surface(mesh, searchEngine, name),
    planeDesc_(planeDesc),
    triangulate_(triangulate),
    points_(0),
    faces_(0),
    meshCells_(0)
{
    createGeometry();
}


// Construct from dictionary. Note how we construct plane from point/normal
// to bypass the elaborate plane description format.
Foam::constantPlane::constantPlane
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    surface(mesh, searchEngine, dict),
    planeDesc_(dict.lookup("basePoint"), dict.lookup("normalVector")),
    triangulate_(getBool(dict, "triangulate", true)),
    points_(0),
    faces_(0),
    meshCells_(0)
{
    createGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantPlane::~constantPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantPlane::correct
(
    const bool meshChanged,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes,
    const fieldsCache<scalar>& scalarCache,
    const fieldsCache<vector>& vectorCache,
    const fieldsCache<tensor>& tensorCache
)
{
    if (meshChanged)
    {
        // Only change of mesh changes plane.
        createGeometry();
    }
}


Foam::scalarField Foam::constantPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<scalar>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return doInterpolate(fieldName, cache);
}


Foam::vectorField Foam::constantPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<vector>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return doInterpolate(fieldName, cache);
}


Foam::tensorField Foam::constantPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<tensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return doInterpolate(fieldName, cache);
}


// ************************************************************************* //

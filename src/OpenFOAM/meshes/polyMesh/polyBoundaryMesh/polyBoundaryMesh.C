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

#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyBoundaryMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make identity map.
labelList polyBoundaryMesh::ident(const label len)
{
    labelList elems(len);
    forAll(elems, elemI)
    {
        elems[elemI] = elemI;
    }
    return elems;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject and a polyMesh reference
polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if (readOpt() == IOobject::MUST_READ)
    {
        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        ptrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.hook
            (
                polyPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const polyMesh&)"
        );

        close();
    }
}


// Construct given size. Patches will be hooked later
polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    polyPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


// Calculate the geometry for the patches (transformation tensors etc.)
void polyBoundaryMesh::calcGeometry()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initGeometry();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).calcGeometry();
    }
}


// Return a list of patch names
wordList polyBoundaryMesh::names() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


// Return a list of patch types
wordList polyBoundaryMesh::types() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


// Return a list of physical types
wordList polyBoundaryMesh::physicalTypes() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].physicalType();
    }

    return t;
}


label polyBoundaryMesh::findPatchID(const word& patchName) const
{
    const polyPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    if (debug)
    {
        Info<< "label polyBoundaryMesh::findPatchID(const word& "
            << "patchName) const"
            << "Patch named " << patchName << " not found.  "
            << "List of available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


// Return patch index for a given face label
label polyBoundaryMesh::whichPatch(const label faceIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (faceIndex >= mesh().nFaces())
    {
        FatalErrorIn
        (
            "polyBoundaryMesh::whichPatch(const label faceIndex) const"
        )   << "given label greater than the number of geometric faces"
            << abort(FatalError);
    }

    if (faceIndex < mesh().nInternalFaces())
    {
        return -1;
    }

    forAll (*this, patchI)
    {
        const polyPatch& bp = operator[](patchI);

        if
        (
            faceIndex >= bp.start()
         && faceIndex < bp.start() + bp.size()
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it's trouble!
    FatalErrorIn
    (
        "label polyBoundaryMesh::whichPatch(const label faceIndex) const"
    )   << "Cannot find face " << faceIndex << " in any of the patches "
        << names() << nl
        << "It seems your patches are not consistent with the mesh :"
        << " internalFaces:" << mesh().nInternalFaces()
        << "  total number of faces:" << mesh().nFaces()
        << abort(FatalError);

    return -1;
}


bool polyBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalFaces();
    const polyBoundaryMesh& bm = *this;

    bool boundaryError = false;

    forAll (bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart)
        {
            boundaryError = true;

            Info<< "bool polyBoundaryMesh::checkDefinition("
                << "const bool report) const : "
                << "Problem with boundary patch " << patchI
                << " named " << bm[patchI].name()
                << " of type " <<  bm[patchI].type()
                << ".\nThe patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchI].start()
                << "." << nl << endl;
        }

        nextPatchStart += bm[patchI].size();
    }

    if (boundaryError)
    {
        SeriousError
            << "bool polyBoundaryMesh::checkDefinition("
            << "const bool report) const : "
            << "This mesh is not valid: boundary definition is in error."
            << endl;
    }
    else
    {
        if (debug || report)
        {
            Info << "Boundary definition OK." << endl;
        }
    }

    return boundaryError;
}


bool polyBoundaryMesh::checkFaceOrder(const bool report) const
{
    // Create dummy topology change.
    polyTopoChange noChange(mesh_);

    // Create map from mesh to itself

    const polyBoundaryMesh& patches = *this;

    labelList oldPatchStarts(patches.size());
    forAll(oldPatchStarts, patchI)
    {
        oldPatchStarts[patchI] = patches[patchI].start();
    }

    // Map from mesh to itself.
    // Note: a few arguments are not set since they are not used in sendOrder.
    // Just lazy ;-)
    mapPolyMesh identiMap
    (
        mesh_,                  // polyMesh
        mesh_.nPoints(),        // nOldPoints
        mesh_.nFaces(),         // nOldFaces
        mesh_.nCells(),         // nOldCells
        ident(mesh_.nPoints()), // pointMap
        ident(mesh_.nFaces()),  // faceMap
        List<objectMap>(),      // facesFromPoints
        List<objectMap>(),      // facesFromEdges
        ident(mesh_.nCells()),  // cellMap
        List<objectMap>(),      // cellsFromPoints
        List<objectMap>(),      // cellsFromEdges
        List<objectMap>(),      // cellsFromFaces
        ident(mesh_.nPoints()), // reversePointMap
        ident(mesh_.nFaces()),  // reverseFaceMap
        ident(mesh_.nCells()),  // reverseCellMap
        labelHashSet(),         // flipFaceFlux
        labelListList(),        // patchPointMap, NOT SET
        labelListList(),        // pointZoneMap, NOT SET
        labelListList(),        // faceZonePointMap, NOT SET
        labelListList(),        // faceZoneFaceMap, NOT SET
        labelListList(),        // cellZoneMap, NOT SET
        mesh_.allPoints(),      // preMotionPoints
        oldPatchStarts,         // oldPatchStarts
        labelList()             // oldPatchNMeshPoints, NOT SET
    );
    
    bool boundaryError = false;

    forAll(patches, patchI)
    {
        patches[patchI].sendOrder(noChange, identiMap);
    }
    
    forAll(patches, patchI)
    {
        labelList pfMap;
        labelList pfRot;

        if (patches[patchI].order(noChange, identiMap, pfMap, pfRot))
        {
            boundaryError = true;

            Sout<< "bool polyBoundaryMesh::checkFaceOrder("
                << "const bool report) const : "
                << "Problem with coupled boundary patch " << patchI
                << " named " << patches[patchI].name()
                << " of type " <<  patches[patchI].type()
                << ".\nThe patch faces are not ordered such that the faces"
                << " on the two sides are numbered identically"
                << "." << nl << endl;
        }
    }

    reduce(boundaryError, orOp<bool>());

    if (boundaryError)
    {
        SeriousError
            << "bool polyBoundaryMesh::checkFaceOrder("
            << "const bool report) const : "
            << "This mesh is not valid: coupled face ordering not consistent."
            << endl;
    }
    else
    {
        if (debug || report)
        {
            Info << "Patch face ordering OK." << endl;
        }
    }

    return boundaryError;
}


// Correct polyBoundaryMesh after moving points
void polyBoundaryMesh::movePoints(const pointField& p)
{
    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initMovePoints(p);
    }

    forAll(patches, patchi)
    {
        patches[patchi].movePoints(p);
    }
}


// writeData member function required by regIOobject
bool polyBoundaryMesh::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const polyBoundaryMesh& patches)
{
    os  << patches.size() << nl << token::BEGIN_LIST;

    forAll(patches, patchI)
    {
        patches[patchI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

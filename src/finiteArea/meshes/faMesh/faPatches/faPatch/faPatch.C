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

#include "faPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faPatch, 0);

defineRunTimeSelectionTable(faPatch, dictionary);

addToRunTimeSelectionTable(faPatch, faPatch, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
faPatch::faPatch
(
    const word& name,
    const labelList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm
)
:
    labelList(edgeLabels),
    patchIdentifier(name, index),
    ngbPolyPatchIndex_(-1),
    boundaryMesh_(bm)
{}


// Construct from dictionary
faPatch::faPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    labelList(dict.lookup("edgeLabels")),
    patchIdentifier(name, dict, index),
    ngbPolyPatchIndex_(readInt(dict.lookup("ngbPolyPatchIndex"))),
    boundaryMesh_(bm)
{}

faPatch::faPatch(const faPatch& p, const faBoundaryMesh& bm)
:
    labelList(p),
    patchIdentifier(p, p.index()),
    ngbPolyPatchIndex_(p.ngbPolyPatchIndex_),
    boundaryMesh_(bm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faPatch::~faPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label faPatch::ngbPolyPatchIndex() const
{
    return ngbPolyPatchIndex_;
}

const faBoundaryMesh& faPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


label faPatch::start() const
{
    return boundaryMesh().mesh().patchStarts()[index()];
}


labelList faPatch::pointLabels() const
{
    SLList<label> labels;

    UList<edge> edges =
        patchSlice(boundaryMesh().mesh().edges());

    forAll(edges, edgeI)
    {
        bool existStart = false;
        bool existEnd = false;

        for
        (
            SLList<label>::iterator iter = labels.begin();
            iter != labels.end();
            ++iter
        )
        {
            if(*iter == edges[edgeI].start())
            {
                existStart = true;
            }

            if(*iter == edges[edgeI].end())
            {
                existEnd = true;
            }
        }

        if(!existStart)
        {
            labels.append(edges[edgeI].start());
        }

        if(!existEnd)
        {
            labels.append(edges[edgeI].end());
        }
    }

    return labels;
}


labelListList faPatch::pointEdges() const
{
    labelList points = pointLabels();

    labelListList pEdges(points.size());

    UList<edge> edges =
        patchSlice(boundaryMesh().mesh().edges());

    forAll(points, pointI)
    {
        SLList<label> curEdges;

        forAll(edges, edgeI)
        {
            if(edges[edgeI].start() == points[pointI] ||
               edges[edgeI].end() == points[pointI])
            {
                curEdges.append(edgeI);
            }
        }

        pEdges[pointI] = curEdges;
    }

    return pEdges;
}


labelList faPatch::ngbPolyPatchFaces() const
{
    labelList ngbFaces;
    
    if(ngbPolyPatchIndex() == -1)
    {
        return ngbFaces;
    }

    ngbFaces.setSize(size());

    const faMesh& aMesh = boundaryMesh().mesh();
    const polyMesh& pMesh = aMesh();
    const indirectPrimitivePatch& patch = aMesh.patch();

    const labelListList& edgeFaces = pMesh.edgeFaces();

    labelList faceCells (patch.size(), -1);

    forAll (faceCells, faceI)
    {
        label faceID = aMesh.faceLabels()[faceI];

        faceCells[faceI] = pMesh.faceOwner()[faceID];
    }

    labelList meshEdges =
        patch.meshEdges
        (
            pMesh.edges(),
            pMesh.cellEdges(),
            faceCells
        );
    
    forAll(ngbFaces, edgeI)
    {
        ngbFaces[edgeI] = -1;

        label curEdge = (*this)[edgeI];

        label curPMeshEdge = meshEdges[curEdge];

        forAll(edgeFaces[curPMeshEdge], faceI)
        {
            label curFace = edgeFaces[curPMeshEdge][faceI];

            label curPatchID = pMesh.boundaryMesh().whichPatch(curFace);

            if (curPatchID == ngbPolyPatchIndex())
            {
                ngbFaces[edgeI] = curFace;
            }            
        }

        if(ngbFaces[edgeI] == -1)
        {
            Info<< "faPatch::edgeNgbPolyPatchFaces(): "
                << "Problem with determination of edge ngb faces!" << endl;
        }
    }

    return ngbFaces;
}


tmp<vectorField> faPatch::ngbPolyPatchFaceNormals() const
{
    tmp<vectorField> fN(new vectorField());

    if (ngbPolyPatchIndex() == -1)
    {
        return fN;
    }

    fN().setSize(size());

    labelList ngbFaces = ngbPolyPatchFaces();

    const polyMesh& pMesh = boundaryMesh().mesh()();

    const faceList& faces = pMesh.faces();
    const pointField& points = pMesh.points();

    forAll(fN(), faceI)
    {
        fN() = faces[ngbFaces[faceI]].normal(points)
            /faces[ngbFaces[faceI]].mag(points);
    }

    return fN;    
}


tmp<vectorField> faPatch::ngbPolyPatchPointNormals() const
{
    if (ngbPolyPatchIndex() == -1)
    {
        return tmp<vectorField>(new vectorField());
    }

    labelListList pntEdges = pointEdges();

    tmp<vectorField> pN(new vectorField(pntEdges.size(), vector::zero));

    vectorField faceNormals = ngbPolyPatchFaceNormals();

    forAll(pN(), pointI)
    {
        forAll(pntEdges[pointI], edgeI)
        {
            pN()[pointI] += faceNormals[pntEdges[pointI][edgeI]];
        }
    }

    pN() /= mag(pN());

    return pN;
}


labelList::subList faPatch::edgeFaces() const
{
//     return labelList::subList
//     (
//         boundaryMesh().mesh()().edgeOwner(),
//         size(),
//         start()
//     );

    return  patchSlice(boundaryMesh().mesh().edgeOwner());
}


// Return the patch edge centres
const vectorField& faPatch::edgeCentres() const
{
    return boundaryMesh().mesh().centres().boundaryField()[index()];
}


// Return the patch edges length vectors
const vectorField& faPatch::edgeLengths() const
{
    return boundaryMesh().mesh().Le().boundaryField()[index()];
}


// Return the patch edge length magnitudes
const scalarField& faPatch::magEdgeLengths() const
{
    return boundaryMesh().mesh().magLe().boundaryField()[index()];
}


// Return the patch edge unit normals
tmp<vectorField> faPatch::edgeNormals() const
{
    tmp<vectorField> eN(new vectorField(size()));

    eN() = edgeLengths()/magEdgeLengths();

    return eN;
}


// Return the patch edge neighbour face centres
tmp<vectorField> faPatch::edgeFaceCentres() const
{
    tmp<vectorField> tfc(new vectorField(size()));
    vectorField& fc = tfc();

    // get reference to global face centres
    const vectorField& gfc = boundaryMesh().mesh().centres().internalField();

    const labelList::subList faceLabels = edgeFaces();

    forAll (faceLabels, edgeI)
    {
        fc[edgeI] = gfc[faceLabels[edgeI]];
    }

    return tfc;
}


// Return cell-centre to face-centre vector
tmp<vectorField> faPatch::delta() const
{
    return edgeCentres() - edgeFaceCentres();
}


// Make delta coefficients as patch face - neighbour cell distances
void faPatch::makeDeltaCoeffs(scalarField& dc) const
{
    dc = 1.0/(edgeNormals() & delta());
}


// Return delta coefficients
const scalarField& faPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


void faPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


const scalarField& faPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


void faPatch::movePoints(const pointField& points)
{}


void faPatch::write(Ostream& os) const
{
    os  << nl << type()
        << static_cast<const patchIdentifier&>(*this) << endl
        << static_cast<const labelList&>(*this) << endl
        << ngbPolyPatchIndex_;
}


void faPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << ';' << nl;

    patchIdentifier::writeDict(os);

    os  << "    edgeLabels " << static_cast<const labelList&>(*this) << ';'
        << nl;

    os  << "    ngbPolyPatchIndex " << ngbPolyPatchIndex_ << ';' << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const faPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const faPatch& p");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

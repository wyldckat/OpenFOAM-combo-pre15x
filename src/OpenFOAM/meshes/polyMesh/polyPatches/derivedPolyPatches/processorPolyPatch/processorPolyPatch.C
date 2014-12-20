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

\*---------------------------------------------------------------------------*/

#include "processorPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorPolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, processorPolyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct from components
processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
)
:
    coupledPolyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


// Construct from Istream
processorPolyPatch::processorPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(is, index, bm),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is)),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


// Construct from dictionary
processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


//- Construct as copy, resetting the boundary mesh
processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


//- Construct as copy, resetting the face list and boundary mesh data
processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorPolyPatch::~processorPolyPatch()
{
    deleteDemandDrivenData(neighbPointsPtr_);
    deleteDemandDrivenData(neighbEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void processorPolyPatch::initGeometry()
{
    if (Pstream::parRun())
    {
        OPstream toNeighbProc
        (
            neighbProcNo(),
          + 3*(sizeof(label) + size()*sizeof(vector))
        );

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void processorPolyPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc
            (
                neighbProcNo(),
                3*(sizeof(label) + size()*sizeof(vector))
            );
            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }

        scalarField magSf = mag(faceAreas());

        forAll(magSf, facei)
        {
            scalar nmagSf = mag(neighbFaceAreas_[facei]);
            scalar avSf = (magSf[facei] + nmagSf)/2.0;

            if (avSf > VSMALL && mag(magSf[facei] - nmagSf)/avSf > 1e-6)
            {
                FatalErrorIn
                (
                    "processorPolyPatch::calcGeometry()"
                )   << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf[facei] - nmagSf)/avSf
                    << "% -- possible face ordering problem"
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            faceCentres(),
            neighbFaceCentres_,
            faceNormals(),
            neighbFaceAreas_/(mag(neighbFaceAreas_)+VSMALL)
        );
    }
}


void processorPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    initGeometry();
}


void processorPolyPatch::movePoints(const pointField&)
{
    calcGeometry();
}


void processorPolyPatch::initUpdateTopology()
{
    // For completeness
    polyPatch::initUpdateTopology();

    deleteDemandDrivenData(neighbPointsPtr_);
    deleteDemandDrivenData(neighbEdgesPtr_);

    if (Pstream::parRun() && boundaryMesh().extendedAddressing())
    {
        // Express all points as patch face and index in face.
        labelList patchFace(nPoints());
        labelList indexInFace(nPoints());

        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
        {
            label faceI = pointFaces()[patchPointI][0];

            patchFace[patchPointI] = faceI;

            const face& f = localFaces()[faceI];

            indexInFace[patchPointI] = findIndex(f, patchPointI);
        }

        OPstream toNeighbProc
        (
            neighbProcNo(),
            3*sizeof(label)
          + 2*nPoints()*sizeof(label)
          + nEdges()*sizeof(edge)
        );

        toNeighbProc
            << patchFace
            << indexInFace
            << edges();
    }
}


void processorPolyPatch::updateTopology()
{
    // For completeness
    polyPatch::updateTopology();

    if (Pstream::parRun() && boundaryMesh().extendedAddressing())
    {
        labelList nbrPatchFace(nPoints());
        labelList nbrIndexInFace(nPoints());
        edgeList nbrEdges(nEdges());

        {
            // Note cannot predict exact size since edgeList not (yet) sent as
            // binary entity but as List of edges.
            IPstream fromNeighbProc(neighbProcNo());

            fromNeighbProc
                >> nbrPatchFace
                >> nbrIndexInFace
                >> nbrEdges;
        }

        // Convert neighbour edges and indices into face back into my edges and
        // points.
        neighbPointsPtr_ = new labelList(nPoints());
        labelList& neighbPoints = *neighbPointsPtr_;

        // Inverse of neighbPoints so from neighbour point to current point.
        labelList nbrToThis(nPoints(), -1);

        forAll(nbrPatchFace, nbrPointI)
        {
            // Find face and index in face on this side.
            const face& f = localFaces()[nbrPatchFace[nbrPointI]];
            label index = (f.size() - nbrIndexInFace[nbrPointI]) % f.size();
            label patchPointI = f[index];

            neighbPoints[patchPointI] = nbrPointI;
            nbrToThis[nbrPointI] = patchPointI;
        }

        // Convert edges.
        neighbEdgesPtr_ = new labelList(nEdges());
        labelList& neighbEdges = *neighbEdgesPtr_;

        forAll(nbrEdges, nbrEdgeI)
        {
            const edge& nbrEdge = nbrEdges[nbrEdgeI];

            // Get edge in local point numbering
            edge e(nbrToThis[nbrEdge[0]], nbrToThis[nbrEdge[1]]);

            // Find the edge.
            const labelList& pEdges = pointEdges()[e[0]];

            label edgeI = -1;

            forAll(pEdges, i)
            {
                if (edges()[pEdges[i]] == e)
                {
                    edgeI = pEdges[i];
                    break;
                }
            }

            if (edgeI == -1)
            {
                FatalErrorIn("processorPolyPatch::calcGeometry()")
                    << "Cannot find patch edge with vertices " << e
                    << '.' << nl << "Can only find edges "
                    << IndirectList<edge>(edges(), pEdges)
                    << " connected to first vertex" << abort(FatalError);
            }

            neighbEdges[edgeI] = nbrEdgeI;
        }
    }
}


const labelList& processorPolyPatch::neighbPoints() const
{
    if (!neighbPointsPtr_)
    {
        FatalErrorIn("processorPolyPatch::neighbPoints() const")
            << "No extended addressing calculated. Please call"
            << " polyBoundaryMesh::setExtendedAddressing(true)"
            << " on all processors before calling this function"
            << abort(FatalError);
    }
    return *neighbPointsPtr_;
}


const labelList& processorPolyPatch::neighbEdges() const
{
    if (!neighbEdgesPtr_)
    {
        FatalErrorIn("processorPolyPatch::neighbEdges() const")
            << "No extended addressing calculated. Please call"
            << " polyBoundaryMesh::setExtendedAddressing(true)"
            << " on all processors before calling this function"
            << abort(FatalError);
    }
    return *neighbEdgesPtr_;
}


// Write
void processorPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);

    os  << nl << myProcNo_ << token::SPACE << neighbProcNo_;
}

void processorPolyPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;
    patchIdentifier::writeDict(os);
    os  << "    nFaces " << size() << token::END_STATEMENT << nl
        << "    startFace " << start() << token::END_STATEMENT << nl
        << "    myProcNo " << myProcNo_ << token::END_STATEMENT << nl
        << "    neighbProcNo " << neighbProcNo_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

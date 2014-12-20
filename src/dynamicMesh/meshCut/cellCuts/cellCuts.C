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

#include "cellCuts.H"
#include "polyMesh.H"
#include "Time.H"
#include "ListOps.H"
#include "cellLooper.H"
#include "refineCell.H"
#include "meshTools.H"
#include "geomCellLooper.H"
#include "OFstream.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cellCuts, 0);

}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

// Find val in first nElems elements of list.
Foam::label Foam::cellCuts::findPartIndex
(
    const labelList& elems,
    const label nElems,
    const label val
)
{
    for (label i = 0; i < nElems; i++)
    {
        if (elems[i] == val)
        {
            return i;
        }
    }
    return -1;
}


Foam::labelList Foam::cellCuts::makeIdent(const label size)
{
    labelList result(size);

    forAll(result, i)
    {
        result[i] = i;
    }
    return result;
}


Foam::boolList Foam::cellCuts::expand
(
    const label size,
    const labelList& labels
)
{
    boolList result(size, false);

    forAll(labels, labelI)
    {
        result[labels[labelI]] = true;
    }
    return result;
}


Foam::scalarField Foam::cellCuts::expand
(
    const label size,
    const labelList& labels,
    const scalarField& weights
)
{
    scalarField result(size, -GREAT);

    forAll(labels, labelI)
    {
        result[labels[labelI]] = weights[labelI];
    }
    return result;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCuts::writeOBJ
(
    const fileName& dir,
    const label cellI,
    const pointField& loopPoints,
    const labelList& anchors
) const
{
    //- Cell edges
    fileName cutsFile(dir / "cell_" + name(cellI) + ".obj");

    Info<< "Writing cell for time " <<  mesh().time().timeName()
        << " to " << cutsFile << nl;

    OFstream cutsStream(cutsFile);

    meshTools::writeOBJ
    (
        cutsStream,
        mesh().cells(),
        mesh().faces(),
        mesh().points(),
        labelList(1, cellI)
    );


    //- Loop cutting cell in two
    fileName loopFile(dir / "cellLoop_" + name(cellI) + ".obj");

    Info<< "Writing loop for time " <<  mesh().time().timeName()
        << " to " << loopFile << nl;

    OFstream loopStream(loopFile);

    label vertI = 0;

    writeOBJ(loopStream, loopPoints, vertI);


    //- Anchors for cell
    fileName anchorFile(dir / "anchors_" + name(cellI) + ".obj");

    Info<< "Writing anchors for time " <<  mesh().time().timeName()
        << " to " << anchorFile << endl;

    OFstream anchorStream(anchorFile);

    forAll(anchors, i)
    {
        meshTools::writeOBJ(anchorStream, mesh().points()[anchors[i]]);
    }
}


// Find face on cell using the two edges.
Foam::label Foam::cellCuts::edgeEdgeToFace
(
    const label cellI,
    const label edgeA,
    const label edgeB
) const
{
    const labelList& cFaces = mesh().cells()[cellI];

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        const labelList& fEdges = mesh().faceEdges()[faceI];

        if
        (
            findIndex(fEdges, edgeA) != -1
         && findIndex(fEdges, edgeB) != -1
        )
        {
           return faceI;
        }
    }

    // Coming here means the loop is illegal since the two edges
    // are not shared by a face. We just mark loop as invalid and continue.

    WarningIn
    (
        "Foam::cellCuts::edgeEdgeToFace"
        "(const label cellI, const label edgeA,"
        "const label edgeB) const"
    )   << "cellCuts : Cannot find face on cell "
        << cellI << " that has both edges " << edgeA << ' ' << edgeB << endl
        << "faces : " << cFaces << endl
        << "edgeA : " << mesh().edges()[edgeA] << endl
        << "edgeB : " << mesh().edges()[edgeB] << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


// Find face on cell using an edge and a vertex.
Foam::label Foam::cellCuts::edgeVertexToFace
(
    const label cellI,
    const label edgeI,
    const label vertI
) const
{
    const labelList& cFaces = mesh().cells()[cellI];

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        const face& f = mesh().faces()[faceI];

        const labelList& fEdges = mesh().faceEdges()[faceI];

        if
        (
            findIndex(fEdges, edgeI) != -1
         && findIndex(f, vertI) != -1
        )
        {
           return faceI;
        }
    }

    WarningIn
    (
        "Foam::cellCuts::edgeVertexToFace"
        "(const label cellI, const label edgeI, "
        "const label vertI) const"
    )   << "cellCuts : Cannot find face on cell "
        << cellI << " that has both edge " << edgeI << " and vertex "
        << vertI << endl
        << "faces : " << cFaces << endl
        << "edge : " << mesh().edges()[edgeI] << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


// Find face using two vertices (guaranteed not to be along edge)
Foam::label Foam::cellCuts::vertexVertexToFace
(
    const label cellI,
    const label vertA,
    const label vertB
) const
{
    const labelList& cFaces = mesh().cells()[cellI];

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        const face& f = mesh().faces()[faceI];

        if
        (
            findIndex(f, vertA) != -1
         && findIndex(f, vertB) != -1
        )
        {
           return faceI;
        }
    }

    WarningIn("Foam::cellCuts::vertexVertexToFace")
        << "cellCuts : Cannot find face on cell "
        << cellI << " that has vertex " << vertA << " and vertex "
        << vertB << endl
        << "faces : " << cFaces << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


void Foam::cellCuts::calcFaceCuts() const
{
    if (faceCutsPtr_)
    {
        FatalErrorIn("cellCuts::calcFaceCuts()")
            << "faceCuts already calculated" << abort(FatalError);
    }

    const faceList& faces = mesh().faces();


    faceCutsPtr_ = new labelListList(mesh().nFaces());
    labelListList& faceCuts = *faceCutsPtr_;

    for (label faceI = 0; faceI < mesh().nFaces(); faceI++)
    {
        const face& f = faces[faceI];

        // Big enough storage. Shrink later on.
        labelList& cuts = faceCuts[faceI];

        cuts.setSize(2*f.size());

        label cutI = 0;

        // Do point-edge-point walk over face and collect all cuts.
        // Problem is that we want to start from one of the endpoints of a 
        // string of connected cuts; we don't want to start somewhere in the
        // middle.

        // Pass1: find first point cut not preceeded by a cut.
        label startFp = -1;

        forAll(f, fp)
        {
            label v0 = f[fp];

            if (pointIsCut_[v0])
            {
                label fpMin1 = (fp == 0 ? f.size()-1 : fp-1);
                label vMin1 = f[fpMin1];

                if 
                (
                    !pointIsCut_[vMin1]
                 && !edgeIsCut_[findEdge(faceI, v0, vMin1)]
                )
                {
                    cuts[cutI++] = vertToEVert(v0);
                    startFp = (fp+1) % f.size();
                    break;
                }
            }
        }

        // Pass2: first edge cut not preceeded by point cut
        if (startFp == -1)
        {
            forAll(f, fp)
            {
                label fp1 = (fp+1) % f.size();

                label v0 = f[fp];
                label v1 = f[fp1];

                label edgeI = findEdge(faceI, v0, v1);

                if (edgeIsCut_[edgeI] && !pointIsCut_[v0])
                {
                    cuts[cutI++] = edgeToEVert(edgeI);
                    startFp = fp1;
                    break;
                }
            }
        }

        // Pass3: nothing found so far. Either face is not cut at all or all
        // vertices are cut. Start from 0.
        if (startFp == -1)
        {
            startFp = 0;
        }

        // Store all cuts starting from startFp;
        label fp = startFp;

        bool allVerticesCut = true;

        forAll(f, i)
        {
            label fp1 = (fp+1) % f.size();

            // Get the three items: current vertex, next vertex and edge
            // inbetween
            label v0 = f[fp];
            label v1 = f[fp1];
            label edgeI = findEdge(faceI, v0, v1);

            if (pointIsCut_[v0])
            {
                cuts[cutI++] = vertToEVert(v0);
            }
            else
            {
                // For check if all vertices have been cut (= illegal)
                allVerticesCut = false;
            }

            if (edgeIsCut_[edgeI])
            {
                cuts[cutI++] = edgeToEVert(edgeI);
            }

            fp = fp1;
        }


        if (allVerticesCut)
        {
            WarningIn("Foam::cellCuts::calcFaceCuts() const")
                << "Face " << faceI << " vertices " << f
                << " has all its vertices cut. Not cutting face." << endl;

            cutI = 0;
        }

        // Remove duplicate starting point
        if (cutI > 1 && cuts[cutI-1] == cuts[0])
        {
            cutI--;
        }
        cuts.setSize(cutI);
    }
}


// Find edge on face using two vertices
Foam::label Foam::cellCuts::findEdge
(
    const label faceI,
    const label v0,
    const label v1
) const
{
    const edgeList& edges = mesh().edges();

    const labelList& fEdges = mesh().faceEdges()[faceI];

    forAll(fEdges, i)
    {
        const edge& e = edges[fEdges[i]];

        if
        (
            (e[0] == v0 && e[1] == v1)
         || (e[0] == v1 && e[1] == v0)
        )
        {
            return fEdges[i];
        }
    }

    return -1;
}


// Check if there is a face on the cell on which all cuts are.
Foam::label Foam::cellCuts::loopFace
(
    const label cellI,
    const labelList& loop
) const
{
    const cell& cFaces = mesh().cells()[cellI];

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        const labelList& fEdges = mesh().faceEdges()[faceI];
        const face& f = mesh().faces()[faceI];

        bool allOnFace = true;

        forAll(loop, i)
        {
            label cut = loop[i];

            if (isEdge(cut))
            {
                if (findIndex(fEdges, getEdge(cut)) == -1)
                {
                    // Edge not on face. Skip face.
                    allOnFace = false;
                    break;
                }
            }
            else
            {
                if (findIndex(f, getVertex(cut)) == -1)
                {
                    // Vertex not on face. Skip face.
                    allOnFace = false;
                    break;
                }
            }
        }

        if (allOnFace)
        {
            // Found face where all elements of loop are on the face.
            return faceI;
        }
    }
    return -1;
}


// From point go into connected face
bool Foam::cellCuts::walkPoint
(
    const label cellI,
    const label startCut,

    const label exclude0,
    const label exclude1,

    const label otherCut,

    label& nVisited,
    labelList& visited
) const
{
    label vertI = getVertex(otherCut);

    const labelList& pFaces = mesh().pointFaces()[vertI];

    forAll(pFaces, pFaceI)
    {
        label otherFaceI = pFaces[pFaceI];

        if
        (
            otherFaceI != exclude0
         && otherFaceI != exclude1
         && meshTools::faceOnCell(mesh(), cellI, otherFaceI)
        )
        {
            label oldNVisited = nVisited;

            bool foundLoop =
                walkCell
                (
                    cellI,
                    startCut,
                    otherFaceI,
                    otherCut,
                    nVisited,
                    visited
                );

            if (foundLoop)
            {
                return true;
            }

            // No success. Restore state and continue
            nVisited = oldNVisited;
        }
    }
    return false;
}


// Cross cut (which is edge on faceI) onto next face
bool Foam::cellCuts::crossEdge
(
    const label cellI,
    const label startCut,
    const label faceI,
    const label otherCut,

    label& nVisited,
    labelList& visited
) const
{
    // Cross edge to other face
    label edgeI = getEdge(otherCut);

    label otherFaceI = meshTools::otherFace(mesh(), cellI, faceI, edgeI);

    // Store old state
    label oldNVisited = nVisited;

    // Recurse to otherCut
    bool foundLoop =
        walkCell
        (
            cellI,
            startCut,
            otherFaceI,
            otherCut,
            nVisited,
            visited
        );

    if (foundLoop)
    {
        return true;
    }
    else
    {    
        // No success. Restore state (i.e. backtrack)
        nVisited = oldNVisited;

        return false;
    }
}


bool Foam::cellCuts::addCut
(
    const label cellI,
    const label cut,
    label& nVisited,
    labelList& visited
) const
{
    // Check for duplicate cuts.
    if (findPartIndex(visited, nVisited, cut) != -1)
    {
        // Truncate (copy of) visited for ease of printing.
        labelList truncVisited(visited);
        truncVisited.setSize(nVisited);

        Info<< "For cell " << cellI << " : found duplicate cut " << cut
            << " in path:";
        writeCuts(Info, truncVisited);
        Info<< endl;

        return false;
    }
    else
    {
        visited[nVisited++] = cut;

        return true;
    }
}


// Walk across faceI, storing cuts as you go. Returns last two cuts visisted.
// Returns true if valid walk.
bool Foam::cellCuts::walkFace
(
    const label cellI,
    const label startCut,
    const label faceI,
    const label cut,

    label& lastCut,
    label& beforeLastCut,
    label& nVisited,
    labelList& visited
) const
{
    const labelList& fCuts = faceCuts()[faceI];

    if (fCuts.size() < 2)
    {
        return false;
    }

    // Easy case : two cuts.
    if (fCuts.size() == 2)
    {
        if (fCuts[0] == cut)
        {
            if (!addCut(cellI, cut, nVisited, visited))
            {
                return false;
            }

            beforeLastCut = cut;
            lastCut = fCuts[1];

            return true;
        }
        else
        {
            if (!addCut(cellI, cut, nVisited, visited))
            {
                return false;
            }

            beforeLastCut = cut;
            lastCut = fCuts[0];

            return true;
        }
    }

    // Harder case: more than 2 cuts on face.
    // Should be start or end of string of cuts. Store all but last cut.
    if (fCuts[0] == cut)
    {
        // Walk forward
        for (label i = 0; i < fCuts.size()-1; i++)
        {
            if (!addCut(cellI, fCuts[i], nVisited, visited))
            {
                return false;
            }
        }
        beforeLastCut = fCuts[fCuts.size()-2];
        lastCut = fCuts[fCuts.size()-1];
    }
    else if (fCuts[fCuts.size()-1] == cut)
    {
        for (label i = fCuts.size()-1; i >= 1; --i)
        {
            if (!addCut(cellI, fCuts[i], nVisited, visited))
            {
                return false;
            }
        }
        beforeLastCut = fCuts[1];
        lastCut = fCuts[0];
    }
    else
    {
        WarningIn("Foam::cellCuts::walkFace")
            << "In middle of cut. cell:" << cellI << " face:" << faceI
            << " cuts:" << fCuts << " current cut:" << cut << endl;

        return false;
    }

    return true;
}



// Walk across cuts (cut edges or cut vertices) of cell. Stops when hit cut 
// already visited. Returns true when loop of 3 or more vertices found.
bool Foam::cellCuts::walkCell
(
    const label cellI,
    const label startCut,
    const label faceI,
    const label cut,

    label& nVisited,
    labelList& visited
) const
{
    // Walk across face, storing cuts. Return last two cuts visited.
    label lastCut = -1;
    label beforeLastCut = -1;

    bool validWalk = walkFace
    (
        cellI,
        startCut,
        faceI,
        cut,

        lastCut,
        beforeLastCut,
        nVisited,
        visited
    );

    if (!validWalk)
    {
        return false;
    }

    // Check if starting point reached.
    if (lastCut == startCut)
    {
        if (nVisited >= 3)
        {
            if (debug & 2)
            {
                // Truncate (copy of) visited for ease of printing.
                labelList truncVisited(visited);
                truncVisited.setSize(nVisited);

                Info<< "For cell " << cellI << " : found closed path:";
                writeCuts(Info, truncVisited);
                Info<< " closed by " << lastCut << endl;
            }

            return true;
        }
        else
        {
            // Loop but too short. Probably folded back on itself.
            return false;
        }
    }


    // Check last cut and one before that to determine type.

    // From beforeLastCut to lastCut:
    //  - from edge to edge
    //      (- always not along existing edge)
    //  - from edge to vertex
    //      - not along existing edge
    //      - along existing edge: not allowed?
    //  - from vertex to edge
    //      - not along existing edge
    //      - along existing edge. Not allowed. See above.
    //  - from vertex to vertex
    //      - not along existing edge
    //      - along existing edge

    if (isEdge(beforeLastCut))
    {
        if (isEdge(lastCut))
        {
            // beforeLastCut=edge, lastCut=edge.

            // Cross lastCut (=edge) into face not faceI
            return crossEdge
            (
                cellI,
                startCut,
                faceI,
                lastCut,
                nVisited,
                visited
            );
        }
        else
        {
            // beforeLastCut=edge, lastCut=vertex.

            // Go from lastCut to all connected faces (except faceI)
            return walkPoint
            (
                cellI,
                startCut,
                faceI,          // exclude0: don't cross back on current face
                -1,             // exclude1
                lastCut,
                nVisited,
                visited
            );
        }
    }
    else
    {
        if (isEdge(lastCut))
        {
            // beforeLastCut=vertex, lastCut=edge.
            return crossEdge
            (
                cellI,
                startCut,
                faceI,
                lastCut,
                nVisited,
                visited
            );
        }
        else
        {
            // beforeLastCut=vertex, lastCut=vertex. Check if along existing
            // edge.
            label edgeI =
                findEdge
                (
                    faceI,
                    getVertex(beforeLastCut),
                    getVertex(lastCut)
                );

            if (edgeI != -1)
            {
                // Cut along existing edge. So is in fact on two faces.
                // Get faces on both sides of the edge to make
                // sure we dont fold back on to those.

                label f0, f1;
                meshTools::getEdgeFaces(mesh(), cellI, edgeI, f0, f1);

                return walkPoint
                (
                    cellI,
                    startCut,
                    f0,
                    f1,
                    lastCut,
                    nVisited,
                    visited
                );

            }
            else
            {
                // Cross cut across face.
                return walkPoint
                (
                    cellI,
                    startCut,
                    faceI,      // face to exclude
                    -1,         // face to exclude
                    lastCut,
                    nVisited,
                    visited
                );
            }
        }
    }
}


// Determine for every cut cell the loop (= face) it is cut by. Done by starting
// from a cut edge or cut vertex and walking across faces, from cut to cut,
// until starting cut hit.
// If multiple loops are possible across a cell circumference takes the first
// one found.
void Foam::cellCuts::calcCellLoops(const labelList& cutCells)
{
    // Calculate cuts per face.
    const labelListList& allFaceCuts = faceCuts();

    // Per cell the number of faces with valid cuts. Is used as quick
    // rejection to see if cell can be cut.
    labelList nCutFaces(mesh().nCells(), 0);

    forAll(allFaceCuts, faceI)
    {
        const labelList& fCuts = allFaceCuts[faceI];

        if (fCuts.size() == mesh().faces()[faceI].size())
        {
            // Too many cuts on face. WalkCell would get very upset so disable.
            nCutFaces[mesh().faceOwner()[faceI]] = labelMin;

            if (mesh().isInternalFace(faceI))
            {
                nCutFaces[mesh().faceNeighbour()[faceI]] = labelMin;
            }
        }
        else if (fCuts.size() >= 2)
        {
            // Could be valid cut. Update count for owner and neighbour.
            nCutFaces[mesh().faceOwner()[faceI]]++;

            if (mesh().isInternalFace(faceI))
            {
                nCutFaces[mesh().faceNeighbour()[faceI]]++;
            }
        }
    }


    // Stack of visited cuts (nVisited used as stack pointer)
    // Size big enough.
    labelList visited(mesh().nPoints());

    forAll(cutCells, i)
    {
        label cellI = cutCells[i];

        // Quick rejection: has enough faces that are cut?
        if (nCutFaces[cellI] >= 3)
        {
            const labelList& cFaces = mesh().cells()[cellI];

            // Determine the first cut face to start walking from.
            label cutFaceI = -1;
            label firstCut = -1;

            forAll(cFaces, cFaceI)
            {
                const labelList& fCuts = allFaceCuts[cFaces[cFaceI]];

                // Take first cut of multiple on face. Note that in calcFaceCuts
                // we have already made sure this is the start or end of a
                // string of cuts.
                if (fCuts.size() >= 2)
                {
                    cutFaceI = cFaces[cFaceI];
                    firstCut = fCuts[0];

                    break;
                }
            }

            if (cutFaceI != -1)
            {
                // Got possibly cut cell with a cutface and an initial cut
                // on the face.

                // Reset stack of visited cuts.
                label nVisited = 0;

                bool validLoop =
                    walkCell
                    (
                        cellI,
                        firstCut,
                        cutFaceI,
                        firstCut,

                        nVisited,
                        visited
                    );

                if (validLoop)
                {
                    // Copy nVisited elements out of visited (since visited is
                    // never truncated for efficiency reasons)

                    labelList& loop = cellLoops_[cellI];

                    loop.setSize(nVisited);

                    for (label i = 0; i < nVisited; i++)
                    {
                        loop[i] = visited[i];
                    }
                }
                else
                {
                    // Invalid loop. Leave cellLoops_[cellI] zero size which
                    // flags this.
                }
            }
        }
    }
}


// Count number of faces used by anchorPoints since the anchor points are
// used (in meshCutter) to create a new cell which should have a minimum
// of 4 faces (so three without the face created from the loop itself)
bool Foam::cellCuts::checkFaces
(
    const label cellI,
    const labelList& anchorPoints
) const
{
    labelHashSet usedFaces(4*anchorPoints.size());

    forAll(anchorPoints, i)
    {
        label pointI = anchorPoints[i];

        const labelList& pFaces = mesh().pointFaces()[pointI];

        forAll(pFaces, pFaceI)
        {
            if (meshTools::faceOnCell(mesh(), cellI, pFaces[pFaceI]))
            {
                usedFaces.insert(pFaces[pFaceI]);
            }
        }
    }

    if (usedFaces.size() < 3)
    {
        return false;
    }
    else
    {
        return true;
    }
}


// On cell get a point which is not part of the cell loop.
Foam::label Foam::cellCuts::findWalkStart
(
    const labelList& loop,
    const label cellI
) const
{
    const labelList& cPoints = mesh().cellPoints()[cellI];

    forAll(cPoints, i)
    {
        label pointI = cPoints[i];

        if (findIndex(loop, vertToEVert(pointI)) == -1)
        {
            return pointI;
        }
    }

    FatalErrorIn("cellCuts::findWalkStart(const label cellI)")
        << "Can not find point on cell " << cellI
        << " which is not cut by loop " << loop << endl
        << "The cut pattern is illegal" << abort(FatalError);

    return -1;
}


// Walk edges of single cell from starting point and collects mesh vertices.
// Does not cross loop (cutEdge or cutPoint)
void Foam::cellCuts::walkEdges
(
    const label cellI,
    const labelList& loop,

    const label pointI,

    DynamicList<label>& connectedPoints
) const
{
    if
    (
        findIndex(loop, vertToEVert(pointI)) == -1
     && findPartIndex(connectedPoints, connectedPoints.size(), pointI) == -1
    )
    {
        // PointI not in loop and not yet visited.

        connectedPoints.append(pointI);

        const labelList& pEdges = mesh().pointEdges()[pointI];

        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];

            if
            (
                meshTools::edgeOnCell(mesh(), cellI, edgeI)
             && findIndex(loop, edgeToEVert(edgeI)) == -1
            )
            {
                // edgeI not cut by loop so recurse.

                label v2 = mesh().edges()[edgeI].otherVertex(pointI);

                walkEdges(cellI, loop, v2, connectedPoints);
            }
        }
    }
}


// Invert anchor point selection.
Foam::labelList Foam::cellCuts::nonAnchorPoints
(
    const labelList& cellPoints,
    const labelList& anchorPoints,
    const labelList& loop
) const
{
    labelList newElems(cellPoints.size());
    label newElemI = 0;

    forAll(cellPoints, i)
    {
        label pointI = cellPoints[i];

        if
        (
            findIndex(anchorPoints, pointI) == -1 
         && findIndex(loop, vertToEVert(pointI)) == -1
        )
        {
            newElems[newElemI++] = pointI;
        }
    }

    newElems.setSize(newElemI);

    return newElems;
}


//- Check anchor points on 'outside' of loop
bool Foam::cellCuts::loopAnchorConsistent
(
    const label cellI,
    const pointField& loopPts,
    const labelList& anchorPoints
) const
{
    // Create identity face for ease of calculation of normal etc.
    face f(makeIdent(loopPts.size()));

    vector normal = f.normal(loopPts);
    point ctr = f.centre(loopPts);


    // Get average position of anchor points.
    vector avg(vector::zero);

    forAll(anchorPoints, ptI)
    {
        avg += mesh().points()[anchorPoints[ptI]];
    }
    avg /= anchorPoints.size();


    if (((avg - ctr) & normal) > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::cellCuts::postCalcAnchorPoints()
{
    forAll(cellLoops_, cellI)
    {
        const labelList& loop = cellLoops_[cellI];

        if (loop.size() > 0 && cellAnchorPoints_[cellI].size() == 0)
        {
            const labelList& cPoints = mesh().cellPoints()[cellI];

            // Visit (and store in connectedPoints) all points reachable from
            // a point not on the loop without crossing any point/edge in loops
            DynamicList<label> connectedPoints(cPoints.size());
            walkEdges(cellI, loop, findWalkStart(loop, cellI), connectedPoints);
            connectedPoints.shrink();

            // Check if the anchorPoints have the correct orientation
            bool loopOk =
                loopAnchorConsistent
                (
                    cellI,
                    loopPoints(cellI),
                    connectedPoints
                );

            if (loopOk)
            {
                cellAnchorPoints_[cellI].transfer(connectedPoints);
                connectedPoints.clear();
            }
            else
            {
                cellAnchorPoints_[cellI] =
                    nonAnchorPoints
                    (
                        cPoints,
                        connectedPoints,
                        loop
                    );
            }
        }
    }
}


Foam::pointField Foam::cellCuts::loopPoints
(
    const labelList& loop,
    const scalarField& loopWeights
) const
{
    pointField loopPts(loop.size());

    forAll(loop, fp)
    {
        loopPts[fp] = coord(loop[fp], loopWeights[fp]);
    }
    return loopPts;
}


// Returns weights of loop. Inverse of loopPoints.
Foam::scalarField Foam::cellCuts::loopWeights(const labelList& loop) const
{
    scalarField weights(loop.size());

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            weights[fp] = edgeWeight_[edgeI];
        }
        else
        {
            weights[fp] = -GREAT;
        }
    }
    return weights;
}


// Check if cut edges in loop are compatible with ones in edgeIsCut_
bool Foam::cellCuts::validEdgeLoop
(
    const labelList& loop,
    const scalarField& loopWeights
) const
{
    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            // Check: cut compatible only if can be snapped to existing one.
            if (edgeIsCut_[edgeI])
            {
                scalar edgeLen = mesh().edges()[edgeI].mag(mesh().points());

                if
                (
                    mag(loopWeights[fp] - edgeWeight_[edgeI])
                  > geomCellLooper::snapTol()*edgeLen
                )
                {
                    // Positions differ too much->would create two cuts.
                    return false;
                }
            }
        }
    }
    return true;
}


// Counts cuts on face. Includes cuts through vertices and through edges.
// Assumes that if edge is cut both in edgeIsCut and in loop that the position
// of the cut is the same.
Foam::label Foam::cellCuts::countFaceCuts
(
    const label faceI,
    const labelList& loop
) const
{
    label nCuts = 0;

    // Count cut vertices
    const face& f = mesh().faces()[faceI];

    forAll(f, fp)
    {
        label vertI = f[fp];

        // Vertex already cut or mentioned in current loop.
        if
        (
            pointIsCut_[vertI]
         || (findIndex(loop, vertToEVert(vertI)) != -1)
        )
        {
            nCuts++;
        }
    }

    // Count cut edges.
    const labelList& fEdges = mesh().faceEdges()[faceI];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];

        // Edge already cut or mentioned in current loop.
        if
        (
            edgeIsCut_[edgeI]
         || (findIndex(loop, edgeToEVert(edgeI)) != -1)
        )
        {
            nCuts++;
        }
    }

    return nCuts;
}


// Determine compatibility of loop with existing cut pattern. Does not use
// cut-addressing (faceCuts_, cutCuts_)
bool Foam::cellCuts::conservativeValidLoop
(
    const label cellI,
    const labelList& loop
) const
{

    if (loop.size() < 2)
    {
        return false;
    }

    forAll(loop, cutI)
    {
        if (isEdge(loop[cutI]))
        {
            label edgeI = getEdge(loop[cutI]);

            if (edgeIsCut_[edgeI])
            {
                // edge compatibility already checked.
            }
            else
            {
                // Quick rejection: vertices of edge itself cannot be cut.
                const edge& e = mesh().edges()[edgeI];

                if (pointIsCut_[e.start()] || pointIsCut_[e.end()])
                {
                    return false;
                }


                // Check faces using this edge
                const labelList& eFaces = mesh().edgeFaces()[edgeI];

                forAll(eFaces, eFaceI)
                {
                    label nCuts = countFaceCuts(eFaces[eFaceI], loop);

                    if (nCuts > 2)
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            // Vertex cut

            label vertI = getVertex(loop[cutI]);

            if (!pointIsCut_[vertI])
            {
                // New cut through vertex.

                // Check edges using vertex.
                const labelList& pEdges = mesh().pointEdges()[vertI];

                forAll(pEdges, pEdgeI)
                {
                    label edgeI = pEdges[pEdgeI];

                    if (edgeIsCut_[edgeI])
                    {
                        return false;
                    }
                }

                // Check faces using vertex.
                const labelList& pFaces = mesh().pointFaces()[vertI];

                forAll(pFaces, pFaceI)
                {
                    label nCuts = countFaceCuts(pFaces[pFaceI], loop);

                    if (nCuts > 2)
                    {
                        return false;
                    }
                }
            }
        }
    }

    // All ok.
    return true;
}


// Determine compatibility of loop with existing cut pattern. Does not use
// derived cut-addressing (faceCuts), only pointIsCut, edgeIsCut.
// Adds any cross-cuts found to newFaceSplitCut and sets cell points on
// one side of the loop in anchorPoints.
bool Foam::cellCuts::validLoop
(
    const label cellI,
    const labelList& loop,
    const scalarField& loopWeights,

    Map<edge>& newFaceSplitCut,
    labelList& anchorPoints
) const
{
    if (loop.size() < 2)
    {
        return false;
    }

    if (debug & 4)
    {
        // Allow as fallback the 'old' loop checking where only a single
        // cut per face is allowed.
        if (!conservativeValidLoop(cellI, loop))
        {
            return  false;
        }
    }

    forAll(loop, fp)
    {
        label cut = loop[fp];
        label nextCut = loop[(fp+1) % loop.size()];

        // Label (if any) of face cut (so cut not along existing edge)
        label meshFaceI = -1;

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            // Look one cut ahead to find if it is along existing edge.

            if (isEdge(nextCut))
            {
                // From edge to edge -> cross cut
                label nextEdgeI = getEdge(nextCut);

                // Find face and mark as to be split.
                meshFaceI = edgeEdgeToFace(cellI, edgeI, nextEdgeI);

                if (meshFaceI == -1)
                {
                    // Can't find face using both cut edges.
                    return false;
                }
            }
            else
            {
                // From edge to vertex -> cross cut only if no existing edge.

                label nextVertI = getVertex(nextCut);
                const edge& e = mesh().edges()[edgeI];

                if (e.start() != nextVertI && e.end() != nextVertI)
                {
                    // New edge. Find face and mark as to be split.
                    meshFaceI = edgeVertexToFace(cellI, edgeI, nextVertI);

                    if (meshFaceI == -1)
                    {
                        // Can't find face. Ilegal.
                        return false;
                    }
                }
            }
        }
        else
        {
            // Cut is vertex.
            label vertI = getVertex(cut);

            if (isEdge(nextCut))
            {
                // From vertex to edge -> cross cut only if no existing edge.
                label nextEdgeI = getEdge(nextCut);

                const edge& nextE = mesh().edges()[nextEdgeI];

                if (nextE.start() != vertI && nextE.end() != vertI)
                {
                    // Should be cross cut. Find face.
                    meshFaceI = edgeVertexToFace(cellI, nextEdgeI, vertI);

                    if (meshFaceI == -1)
                    {
                        return false;
                    }
                }
            }
            else
            {
                // From vertex to vertex -> cross cut only if no existing edge.
                label nextVertI = getVertex(nextCut);

                if (meshTools::findEdge(mesh(), vertI, nextVertI) == -1)
                {
                    // New edge. Find face.
                    meshFaceI = vertexVertexToFace(cellI, vertI, nextVertI);

                    if (meshFaceI == -1)
                    {
                        return false;
                    }
                }
            }
        }

        if (meshFaceI != -1)
        {
            // meshFaceI is cut across along cut-nextCut (so not along existing
            // edge). Check if this is compatible with existing pattern.
            edge cutEdge(cut, nextCut);

            Map<edge>::const_iterator iter = faceSplitCut_.find(meshFaceI);

            if (iter == faceSplitCut_.end())
            {
                // Face not yet cut so insert.
                newFaceSplitCut.insert(meshFaceI, cutEdge);
            }
            else
            {
                // Face already cut. Ok if same edge.
                if (iter() != cutEdge)
                {
                    return false;
                }
            }
        }
    }

    // Is there a face on which all cuts are?
    label faceContainingLoop = loopFace(cellI, loop);

    if (faceContainingLoop != -1)
    {
        WarningIn("Foam::cellCuts::validLoop")
            << "Found loop on cell " << cellI << " with all points"
            << " on face " << faceContainingLoop << endl;

        writeOBJ(".", cellI, loopPoints(loop, loopWeights), labelList(0));

        return false;
    }


    // Final check: are more than 2 faces used in each of the two halves of the
    // cell

    // Calculate anchor points

    const labelList& cPoints = mesh().cellPoints()[cellI];

    // Visit (and store in connectedPoints) all points reachable from
    // a point not on the loop without crossing any point/edge in loops
    DynamicList<label> connectedPoints(cPoints.size());
    walkEdges(cellI, loop, findWalkStart(loop, cellI), connectedPoints);
    connectedPoints.shrink();

    // Determine all other non-loop points as well.
    labelList otherPoints(nonAnchorPoints(cPoints, connectedPoints, loop));

    // Points on loop
    pointField loopPts(loopPoints(loop, loopWeights));


    // Check that all cell points are either in loop, in anchor set or in rest.
    {
        label nLoopPoints = 0;

        forAll(loop, i)
        {
            if (!isEdge(loop[i]))
            {
                nLoopPoints++;
            }
        }

        if
        (
                nLoopPoints + connectedPoints.size() + otherPoints.size()
             != cPoints.size()
        )
        {
            writeOBJ(".", cellI, loopPts, connectedPoints);

            FatalErrorIn("validLoop")
                << "For cell:" << cellI
                << " points should be in one of -loop -anchor -nonanchor"
                << endl
                << "cellPoints:" << cPoints << endl
                << "loop:" << loop << endl
                << "anchors:" << connectedPoints << endl
                << "otherPoints:" << otherPoints << endl
                << abort(FatalError);
        }
    }


    // Check if faces split off will have enough vertices
    if (!checkFaces(cellI, connectedPoints))
    {
        WarningIn("Foam::cellCuts::validLoop")
            << "Invalid loop " << loop << " for cell " << cellI
            << " since would have too few faces on one side." << nl
            << "All faces:" << mesh().cells()[cellI] << endl;

        writeOBJ(".", cellI, loopPts, connectedPoints);

        return false;
    }

    if (!checkFaces(cellI, otherPoints))
    {
        WarningIn("Foam::cellCuts::validLoop")
            << "Invalid loop " << loop << " for cell " << cellI
            << " since would have too few faces on one side." << nl
            << "All faces:" << mesh().cells()[cellI] << endl;

        writeOBJ(".", cellI, loopPts, otherPoints);

        return false;
    }

    // Check which one of point sets to use.
    bool loopOk = loopAnchorConsistent(cellI, loopPts, connectedPoints);

    //if (debug)
    {
        // Additional check: are non-anchor points on other side?
        bool otherLoopOk = loopAnchorConsistent(cellI, loopPts, otherPoints);

        if (loopOk == otherLoopOk)
        {
            // Both sets of points are supposedly on the same side as the
            // loop normal. Oops.

            WarningIn("Foam::cellCuts::validLoop")
                << " For cell:" << cellI
                << " achorpoints and nonanchorpoints are geometrically"
                << " on same side!" << endl
                << "cellPoints:" << cPoints << endl
                << "loop:" << loop << endl
                << "anchors:" << connectedPoints << endl
                << "otherPoints:" << otherPoints << endl;
                //<< abort(FatalError);

            writeOBJ(".", cellI, loopPts, connectedPoints);
        }
    }


    if (loopOk)
    {
        // connectedPoints on 'outside' of loop so these are anchor points
        anchorPoints.transfer(connectedPoints);
        connectedPoints.clear();
    }
    else
    {
        anchorPoints.transfer(otherPoints);
    }

    // All ok.
    return true;
}


// Update basic cut information (pointIsCut, edgeIsCut) from cellLoops.
// Assumes cellLoops_ and edgeWeight_ already set and consistent.
// Does not use any other information.
void Foam::cellCuts::setFromCellLoops()
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    faceSplitCut_.clear();

    forAll(cellLoops_, cellI)
    {
        const labelList& loop = cellLoops_[cellI];

        if (loop.size() > 0)
        {
            // Storage for cross-face cuts
            Map<edge> faceSplitCuts(loop.size());
            // Storage for points on one side of cell.
            labelList anchorPoints;

            if
            (
               !validLoop
                (
                    cellI,
                    loop,
                    loopWeights(loop),
                    faceSplitCuts,
                    anchorPoints
                )
            )
            {
                writeOBJ(".", cellI, loopPoints(cellI), anchorPoints);

                //FatalErrorIn("cellCuts::setFromCellLoops()")
                WarningIn("cellCuts::setFromCellLoops")
                    << "Illegal loop " << loop
                    << " when recreating cut-addressing"
                    << " from existing cellLoops" << endl
                    << "cellI:" << cellI << endl;
                    //<< abort(FatalError);

                cellLoops_[cellI].setSize(0);
                cellAnchorPoints_[cellI].setSize(0);
            }
            else
            {
                // Copy anchor points.
                cellAnchorPoints_[cellI].transfer(anchorPoints);

                // Copy faceSplitCuts into overall faceSplit info.
                for
                (
                    Map<edge>::const_iterator iter = faceSplitCuts.begin();
                    iter != faceSplitCuts.end();
                    ++iter
                )
                {
                    faceSplitCut_.insert(iter.key(), iter());
                }

                // Update edgeIsCut, pointIsCut information
                forAll(loop, cutI)
                {
                    label cut = loop[cutI];

                    if (isEdge(cut))
                    {
                        edgeIsCut_[getEdge(cut)] = true;
                    }
                    else
                    {
                        pointIsCut_[getVertex(cut)] = true;
                    }
                }
            }
        }
    }

    // Reset edge weights
    forAll(edgeIsCut_, edgeI)
    {
        if (!edgeIsCut_[edgeI])
        {
            // Weight not used. Set to illegal value to make any use fall over.
            edgeWeight_[edgeI] = -GREAT;
        }
    }
}


// Upate basic cut information from single cellLoop. Returns true if loop
// was valid.
bool Foam::cellCuts::setFromCellLoop
(
    const label cellI,
    const labelList& loop,
    const scalarField& loopWeights
)
{
    // Dump loop for debugging.
    if (debug)
    {
        OFstream str("last_cell.obj");

        str<< "# edges of cell " << cellI << nl;

        meshTools::writeOBJ
        (
            str,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            labelList(1, cellI)
        );


        OFstream loopStr("last_loop.obj");

        loopStr<< "# looppoints for cell " << cellI << nl;

        pointField pointsOfLoop = loopPoints(loop, loopWeights);

        forAll(pointsOfLoop, i)
        {
            meshTools::writeOBJ(loopStr, pointsOfLoop[i]);
        }

        str << 'l';

        forAll(pointsOfLoop, i)
        {
            str << ' ' << i + 1;
        }
        str << ' ' << 1 << nl;
    }

    bool okLoop = false;

    if (validEdgeLoop(loop, loopWeights))
    {
        // Storage for cuts across face
        Map<edge> faceSplitCuts(loop.size());
        // Storage for points on one side of cell
        labelList anchorPoints;

        pointField pointsOfLoop(loopPoints(loop, loopWeights));

        okLoop = 
            validLoop(cellI, loop, loopWeights, faceSplitCuts, anchorPoints);

        if (okLoop)
        {
            // Valid loop. Copy cellLoops and anchorPoints
            cellLoops_[cellI] = loop;
            cellAnchorPoints_[cellI].transfer(anchorPoints);

            // Copy split cuts
            for
            (
                Map<edge>::const_iterator iter = faceSplitCuts.begin();
                iter != faceSplitCuts.end();
                ++iter
            )
            {
                faceSplitCut_.insert(iter.key(), iter());
            }


            // Update edgeIsCut, pointIsCut information
            forAll(loop, cutI)
            {
                label cut = loop[cutI];

                if (isEdge(cut))
                {
                    label edgeI = getEdge(cut);

                    edgeIsCut_[edgeI] = true;

                    edgeWeight_[edgeI] = loopWeights[cutI];
                }
                else
                {
                    label vertI = getVertex(cut);

                    pointIsCut_[vertI] = true;
                }
            }
        }
    }

    return okLoop;
}


// Update basic cut information from cellLoops. Checks for consistency with
// existing cut pattern.
void Foam::cellCuts::setFromCellLoops
(
    const labelList& cellLabels,
    const labelListList& cellLoops,
    const List<scalarField>& cellLoopWeights
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    forAll(cellLabels, cellLabelI)
    {
        label cellI = cellLabels[cellLabelI];

        const labelList& loop = cellLoops[cellLabelI];

        if (loop.size() > 0)
        {
            const scalarField& loopWeights = cellLoopWeights[cellLabelI];

            if (setFromCellLoop(cellI, loop, loopWeights))
            {
                // Valid loop. Call above will have upated all already.
            }
            else
            {
                // Clear cellLoops
                cellLoops_[cellI].setSize(0);
            }
        }
    }
}


// Cut cells and update basic cut information from cellLoops. Checks each loop
// for consistency with existing cut pattern.
void Foam::cellCuts::setFromCellCutter
(
    const cellLooper& cellCutter,
    const List<refineCell>& refCells
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    // storage for loop of cuts (cut vertices and/or cut edges)
    labelList cellLoop;
    scalarField cellLoopWeights;

    // For debugging purposes
    DynamicList<label> invalidCutCells(2);
    DynamicList<labelList> invalidCutLoops(2);
    DynamicList<scalarField> invalidCutLoopWeights(2);


    forAll(refCells, refCellI)
    {
        const refineCell& refCell = refCells[refCellI];

        label cellI = refCell.cellNo();

        const vector& refDir = refCell.direction();


        // Cut cell. Determines cellLoop and cellLoopWeights
        bool goodCut =
            cellCutter.cut
            (
                refDir,
                cellI,

                pointIsCut_,
                edgeIsCut_,
                edgeWeight_,

                cellLoop,
                cellLoopWeights
            );

        // Check whether edge refinement is on a per face basis compatible with
        // current pattern.
        if (goodCut)
        {
            if (setFromCellLoop(cellI, cellLoop, cellLoopWeights))
            {
                // Valid loop. Will have updated all info already.
            }
            else
            {
                cellLoops_[cellI].setSize(0);

                // Discarded by validLoop
                if (debug)
                {
                    invalidCutCells.append(cellI);
                    invalidCutLoops.append(cellLoop);
                    invalidCutLoopWeights.append(cellLoopWeights);
                }
            }
        }
        else
        {
            // Clear cellLoops
            cellLoops_[cellI].setSize(0);
        }
    }

    if (debug && invalidCutCells.size() > 0)
    {
        invalidCutCells.shrink();
        invalidCutLoops.shrink();
        invalidCutLoopWeights.shrink();

        fileName cutsFile("invalidLoopCells.obj");

        Info<< "cellCuts : writing inValidLoops cells to " << cutsFile << endl;

        OFstream cutsStream(cutsFile);

        meshTools::writeOBJ
        (
            cutsStream,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            invalidCutCells
        );

        fileName loopsFile("invalidLoops.obj");

        Info<< "cellCuts : writing inValidLoops loops to " << loopsFile << endl;

        OFstream loopsStream(loopsFile);

        label vertI = 0;

        forAll(invalidCutLoops, i)
        {
            writeOBJ
            (
                loopsStream,
                loopPoints(invalidCutLoops[i], invalidCutLoopWeights[i]),
                vertI
            );
        }
    }
}


// Same as one before but cut plane prescribed (instead of just normal)
void Foam::cellCuts::setFromCellCutter
(
    const cellLooper& cellCutter,
    const labelList& cellLabels,
    const List<plane>& cellCutPlanes
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    // storage for loop of cuts (cut vertices and/or cut edges)
    labelList cellLoop;
    scalarField cellLoopWeights;

    // For debugging purposes
    DynamicList<label> invalidCutCells(2);
    DynamicList<labelList> invalidCutLoops(2);
    DynamicList<scalarField> invalidCutLoopWeights(2);


    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        // Cut cell. Determines cellLoop and cellLoopWeights
        bool goodCut =
            cellCutter.cut
            (
                cellCutPlanes[i],
                cellI,

                pointIsCut_,
                edgeIsCut_,
                edgeWeight_,

                cellLoop,
                cellLoopWeights
            );

        // Check whether edge refinement is on a per face basis compatible with
        // current pattern.
        if (goodCut)
        {
            if (setFromCellLoop(cellI, cellLoop, cellLoopWeights))
            {
                // Valid loop. Will have updated all info already.
            }
            else
            {
                cellLoops_[cellI].setSize(0);

                // Discarded by validLoop
                if (debug)
                {
                    invalidCutCells.append(cellI);
                    invalidCutLoops.append(cellLoop);
                    invalidCutLoopWeights.append(cellLoopWeights);
                }
            }
        }
        else
        {
            // Clear cellLoops
            cellLoops_[cellI].setSize(0);
        }
    }

    if (debug && invalidCutCells.size() > 0)
    {
        invalidCutCells.shrink();
        invalidCutLoops.shrink();
        invalidCutLoopWeights.shrink();

        fileName cutsFile("invalidLoopCells.obj");

        Info<< "cellCuts : writing inValidLoops cells to " << cutsFile << endl;

        OFstream cutsStream(cutsFile);

        meshTools::writeOBJ
        (
            cutsStream,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            invalidCutCells
        );

        fileName loopsFile("invalidLoops.obj");

        Info<< "cellCuts : writing inValidLoops loops to " << loopsFile << endl;

        OFstream loopsStream(loopsFile);

        label vertI = 0;

        forAll(invalidCutLoops, i)
        {
            writeOBJ
            (
                loopsStream,
                loopPoints(invalidCutLoops[i], invalidCutLoopWeights[i]),
                vertI
            );
        }
    }
}


// Set orientation of loops
void Foam::cellCuts::orientPlanesAndLoops()
{
    // Determine anchorPoints if not yet done by validLoop.
    postCalcAnchorPoints();

    if (debug & 2)
    {
        Info<< "cellAnchorPoints:" << endl;
    }
    forAll(cellAnchorPoints_, cellI)
    {
        if (cellLoops_[cellI].size() > 0)
        {
            if (cellAnchorPoints_[cellI].size() == 0)
            {
                FatalErrorIn("orientPlanesAndLoops()")
                    << "No anchor points for cut cell " << cellI << endl
                    << "cellLoop:" << cellLoops_[cellI] << abort(FatalError);
            }

            if (debug & 2)
            {
                Info<< "    cell:" << cellI << " anchored at "
                    << cellAnchorPoints_[cellI] << endl;
            }
        }
    }

    // Calculate number of valid cellLoops
    nLoops_ = 0;

    forAll(cellLoops_, cellI)
    {
        if (cellLoops_[cellI].size() > 0)
        {
            nLoops_++;
        }
    }
}


// Do all: calculate addressing, calculate loops splitting cells
void Foam::cellCuts::calcLoopsAndAddressing(const labelList& cutCells)
{
    // Sanity check on weights
    forAll(edgeIsCut_, edgeI)
    {
        if (edgeIsCut_[edgeI])
        {
            scalar weight = edgeWeight_[edgeI];

            if (weight < 0 || weight > 1)
            {
                FatalErrorIn
                (
                    "cellCuts::calcLoopsAndAddressing(const labelList&)"
                )   << "Weight out of range [0,1]. Edge " << edgeI
                    << " verts:" << mesh().edges()[edgeI]
                    << " weight:" << weight << abort(FatalError);
            }
        }
        else
        {
            // Weight not used. Set to illegal value to make any use fall over.
            edgeWeight_[edgeI] = -GREAT;
        }
    }


    // Calculate faces that split cells in two
    calcCellLoops(cutCells);

    if (debug & 2)
    {
        Info<< "-- cellLoops --" << endl;
        forAll(cellLoops_, cellI)
        {
            const labelList& loop = cellLoops_[cellI];

            if (loop.size() > 0)
            {
                Info<< "cell:" << cellI << "  ";
                writeCuts(Info, loop);
                Info<< endl;
            }
        }
    }

    // Redo basic cut information (pointIsCut, edgeIsCut, faceSplitCut)
    // using cellLoop only.
    setFromCellLoops();
}


void Foam::cellCuts::check() const
{
    // Check weights for unsnapped values
    forAll(edgeIsCut_, edgeI)
    {
        if (edgeIsCut_[edgeI])
        {
            if
            (
                edgeWeight_[edgeI] < geomCellLooper::snapTol()
             || edgeWeight_[edgeI] > (1 - geomCellLooper::snapTol())
            )
            {
                // Should have been snapped.
                //FatalErrorIn("cellCuts::check()")
                WarningIn("cellCuts::check()")
                    << "edge:" << edgeI << " vertices:"
                    << mesh().edges()[edgeI]
                    << " weight:" << edgeWeight_[edgeI] << " should have been"
                    << " snapped to one of its endpoints"
                    //<< abort(FatalError);
                    << endl;
            }
        }
        else
        {
            if (edgeWeight_[edgeI] > - 1)
            {
                FatalErrorIn("cellCuts::check()")
                    << "edge:" << edgeI << " vertices:"
                    << mesh().edges()[edgeI]
                    << " weight:" << edgeWeight_[edgeI] << " is not cut but"
                    << " its weight is not marked invalid"
                    << abort(FatalError);
            }
        }
    }

    // Check that all elements of cellloop are registered
    forAll(cellLoops_, cellI)
    {
        const labelList& loop = cellLoops_[cellI];

        forAll(loop, i)
        {
            label cut = loop[i];

            if
            (
                (isEdge(cut) && !edgeIsCut_[getEdge(cut)])
             || (!isEdge(cut) && !pointIsCut_[getVertex(cut)])
            )
            {
                writeCut(Info, cut);

                FatalErrorIn("cellCuts::check()")
                    << "cell:" << cellI << " loop:"
                    << loop
                    << " cut:" << cut << " is not marked as cut"
                    << abort(FatalError);
            }
        }
    }

    // Check that no elements of cell loop are anchor point.
    forAll(cellLoops_, cellI)
    {
        const labelList& anchors = cellAnchorPoints_[cellI];

        const labelList& loop = cellLoops_[cellI];

        if (loop.size() > 0 && anchors.size() == 0)
        {
            FatalErrorIn("cellCuts::check()")
                << "cell:" << cellI << " loop:" << loop
                << " has no anchor points"
                << abort(FatalError);
        }


        forAll(loop, i)
        {
            label cut = loop[i];

            if
            (
                !isEdge(cut)
             && findIndex(anchors, getVertex(cut)) != -1
            )
            {
                FatalErrorIn("cellCuts::check()")
                    << "cell:" << cellI << " loop:" << loop
                    << " anchor points:" << anchors
                    << " anchor:" << getVertex(cut) << " is part of loop"
                    << abort(FatalError);
            }
        }
    }


    // Check that cut faces have a neighbour that is cut.
    for
    (
        Map<edge>::const_iterator iter = faceSplitCut_.begin();
        iter != faceSplitCut_.end();
        ++iter
    )
    {
        label faceI = iter.key();

        if (mesh().isInternalFace(faceI))
        {
            label own = mesh().faceOwner()[faceI];
            label nei = mesh().faceNeighbour()[faceI];

            if
            (
                cellLoops_[own].size() == 0
             && cellLoops_[nei].size() == 0
            )
            {
                FatalErrorIn("cellCuts::check()")
                    << "Internal face:" << faceI << " cut by " << iter()
                    << " has owner:" << own
                    << " and neighbour:" << nei
                    << " that are both uncut"
                    << abort(FatalError);
            }
        }
        else
        {
            label own = mesh().faceOwner()[faceI];

            if (cellLoops_[own].size() == 0)
            {
                FatalErrorIn("cellCuts::check()")
                    << "Boundary face:" << faceI << " cut by " << iter()
                    << " has owner:" << own
                    << " that is uncut"
                    << abort(FatalError);
            }
        }
    }
}


void Foam::cellCuts::clearOut()
{
    deleteDemandDrivenData(faceCutsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from cells to cut and pattern of cuts
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& cutCells,
    const labelList& meshVerts,
    const labelList& meshEdges,
    const scalarField& meshEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(expand(mesh.nPoints(), meshVerts)),
    edgeIsCut_(expand(mesh.nEdges(), meshEdges)),
    edgeWeight_(expand(mesh.nEdges(), meshEdges, meshEdgeWeights)),
    faceCutsPtr_(NULL),
    faceSplitCut_(cutCells.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Info<< "cellCuts : constructor from cut verts and edges" << endl;
    }

    calcLoopsAndAddressing(cutCells);

    // Calculate planes and flip cellLoops if nessecary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Info<< "cellCuts : leaving constructor from cut verts and edges"
            << endl;
    }
}


// Construct from pattern of cuts. Finds out itself which cells are cut.
// (can go wrong if e.g. all neighbours of cell are refined)
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& meshVerts,
    const labelList& meshEdges,
    const scalarField& meshEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(expand(mesh.nPoints(), meshVerts)),
    edgeIsCut_(expand(mesh.nEdges(), meshEdges)),
    edgeWeight_(expand(mesh.nEdges(), meshEdges, meshEdgeWeights)),
    faceCutsPtr_(NULL),
    faceSplitCut_(mesh.nFaces()/10 + 1),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Info<< "cellCuts : constructor from cellLoops" << endl;
    }

    labelList cutCells(mesh.nCells());

    forAll(cutCells, cellI)
    {
        cutCells[cellI] = cellI;
    }
    calcLoopsAndAddressing(cutCells);

    // Calculate planes and flip cellLoops if nessecary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Info<< "cellCuts : leaving constructor from cellLoops" << endl;
    }
}


// Construct from complete cellLoops. Assumes correct cut pattern.
// Only constructs cut-cut addressing and cellAnchorPoints
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelListList& cellLoops,
    const List<scalarField>& cellEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -GREAT),
    faceCutsPtr_(NULL),
    faceSplitCut_(cellLabels.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Info<< "cellCuts : constructor from cellLoops" << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellLoops(cellLabels, cellLoops, cellEdgeWeights);

    // Calculate planes and flip cellLoops if nessecary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Info<< "cellCuts : leaving constructor from cellLoops" << endl;
    }
}


// Construct from list of cells to cut and cell cutter.
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const cellLooper& cellCutter,
    const List<refineCell>& refCells
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -GREAT),
    faceCutsPtr_(NULL),
    faceSplitCut_(refCells.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Info<< "cellCuts : constructor from cellCutter" << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellCutter(cellCutter, refCells);

    // Calculate planes and flip cellLoops if nessecary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Info<< "cellCuts : leaving constructor from cellCutter" << endl;
    }
}


// Construct from list of cells to cut, plane to cut with and cell cutter.
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const cellLooper& cellCutter,
    const labelList& cellLabels,
    const List<plane>& cutPlanes
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -GREAT),
    faceCutsPtr_(NULL),
    faceSplitCut_(cellLabels.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Info<< "cellCuts : constructor from cellCutter with prescribed plane"
            << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellCutter(cellCutter, cellLabels, cutPlanes);

    // Calculate planes and flip cellLoops if nessecary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Info<< "cellCuts : leaving constructor from cellCutter with prescribed"
            << " plane" << endl;
    }
}


// Construct from components
Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const boolList& pointIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,
    const Map<edge>& faceSplitCut,
    const labelListList& cellLoops,
    const label nLoops,
    const labelListList& cellAnchorPoints
)
:
    edgeVertex(mesh),
    pointIsCut_(pointIsCut),
    edgeIsCut_(edgeIsCut),
    edgeWeight_(edgeWeight),
    faceCutsPtr_(NULL),
    faceSplitCut_(faceSplitCut),
    cellLoops_(cellLoops),
    nLoops_(nLoops),
    cellAnchorPoints_(cellAnchorPoints)
{
    if (debug)
    {
        Info<< "cellCuts : constructor from components" << endl;
        Info<< "cellCuts : leaving constructor from components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCuts::~cellCuts()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::cellCuts::loopPoints(const label cellI) const
{
    const labelList& loop = cellLoops_[cellI];

    pointField loopPts(loop.size());

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            loopPts[fp] = coord(cut, edgeWeight_[edgeI]);
        }
        else
        {
            loopPts[fp] = coord(cut, -GREAT);
        }
    }
    return loopPts;
}


// Flip loop for cell
void Foam::cellCuts::flip(const label cellI)
{
    labelList& loop = cellLoops_[cellI];

    reverse(loop);

    // Reverse anchor point set.
    cellAnchorPoints_[cellI] =
        nonAnchorPoints
        (
            mesh().cellPoints()[cellI],
            cellAnchorPoints_[cellI],
            loop
        );  
}


// Flip loop only
void Foam::cellCuts::flipLoopOnly(const label cellI)
{
    labelList& loop = cellLoops_[cellI];

    reverse(loop);
}


void Foam::cellCuts::writeOBJ
(
    Ostream& os,
    const pointField& loopPts,
    label& vertI
) const
{
    label startVertI = vertI;

    forAll(loopPts, fp)
    {
        const point& pt = loopPts[fp];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;

        vertI++;
    }

    os << 'f';
    forAll(loopPts, fp)
    {
        os  << ' ' << startVertI + fp + 1;
    }
    os<< endl;
}


void Foam::cellCuts::writeOBJ(Ostream& os) const
{
    label vertI = 0;

    forAll(cellLoops_, cellI)
    {
        writeOBJ(os, loopPoints(cellI), vertI);
    }
}


void Foam::cellCuts::writeCellOBJ(const fileName& dir, const label cellI) const
{
    const labelList& anchors = cellAnchorPoints_[cellI];

    writeOBJ(dir, cellI, loopPoints(cellI), anchors);
}


// ************************************************************************* //

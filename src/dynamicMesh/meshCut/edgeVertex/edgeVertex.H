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

Class
    edgeVertex

Description
    Combines edge or vertex in single label. Used to specify cuts across
    cell circumference.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef edgeVertex_H
#define edgeVertex_H

#include "label.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class refineCell;

/*---------------------------------------------------------------------------*\
                           Class edgeVertex Declaration
\*---------------------------------------------------------------------------*/

class edgeVertex
{
    // Private data

        //- Reference to mesh. (could be primitive mesh but keeping polyMesh
        //  here saves storing reference at higher levels where we do need it)
        const polyMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        edgeVertex(const edgeVertex&);

        //- Disallow default bitwise assignment
        void operator=(const edgeVertex&);


public:

    // Static Functions

        //- Update refine list from map. Used to update cell/face labels
        //  after morphing
        static void updateLabels(const labelList& map, List<refineCell>&);

        //- Update map from map. Used to update cell/face labels
        //  after morphing
        static void updateLabels(const labelList& map, Map<label>&);

        //- Update map from map. Used to update cell/face labels
        //  after morphing
        static void updateLabels(const labelList& map, labelHashSet&);



    // Constructors

        //- Construct from mesh
        edgeVertex(const polyMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }


    // EdgeVertex handling

        //- is eVert an edge?
        bool isEdge(const label eVert) const
        {
            if (eVert < 0 || eVert >= (mesh_.nPoints() + mesh_.nEdges()))
            {
                FatalErrorIn("edgeVertex::isEdge(const label)")
                    << "EdgeVertex " << eVert << " out of range "
                    << mesh_.nPoints() << " to "
                    << (mesh_.nPoints() + mesh_.nEdges() - 1)
                    << abort(FatalError);
            }
            
            return eVert >= mesh_.nPoints();
        }

        //- convert eVert to edge label
        label getEdge(const label eVert) const
        {
            if (!isEdge(eVert))
            {
                FatalErrorIn("edgeVertex::getEdge(const label)")
                    << "EdgeVertex " << eVert << " not an edge"
                    << abort(FatalError);
            }
            return eVert - mesh_.nPoints();
        }

        //- convert eVert to vertex label
        label getVertex(const label eVert) const
        {
            if (isEdge(eVert) || (eVert < 0))
            {
                FatalErrorIn("edgeVertex::getVertex(const label)")
                    << "EdgeVertex " << eVert << " not a vertex"
                    << abort(FatalError);
            }            
            return eVert;
        }

        //- Convert pointI to eVert
        label vertToEVert(const label vertI) const
        {
            if ((vertI < 0) || (vertI >= mesh_.nPoints()))
            {
                FatalErrorIn("edgeVertex::vertToEVert(const label)")
                    << "Illegal vertex number " << vertI
                    << abort(FatalError);
            }
            return vertI;
        }

        //- Convert edgeI to eVert
        label edgeToEVert(const label edgeI) const
        {
            if ((edgeI < 0) || (edgeI >= mesh_.nEdges()))
            {
                FatalErrorIn("edgeVertex::edgeToEVert(const label)")
                    << "Illegal edge number " << edgeI
                    << abort(FatalError);
            }
            return mesh_.nPoints() + edgeI;
        }

        //- Return coordinate of cut (uses weight if edgeCut)
        point coord(const label cut, const scalar weight) const;

        //- Find mesh edge (or -1) between two cuts. 
        label cutPairToEdge(const label cut0, const label cut1) const;

        //- Write cut description to Ostream
        Ostream& writeCut(Ostream& os, const label cut, const scalar) const;

        //- Write cut descriptions to Ostream
        Ostream& writeCuts(Ostream& os, const labelList&, const scalarField&)
         const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

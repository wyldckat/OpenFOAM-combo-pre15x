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

Class
    twoDPointCorrector

Description
    Class applies a two-dimensional correction to mesh motion point field.
    The correction guarantees that the mesh does not get twisted during motion
    and thus introduce a third dimension into a 2-D problem.

    The operation is performed by looping through all edges approximately
    normal to the plane and enforcing their orthoginality onto the plane by
    adjusting points on their ends.

SourceFiles
    twoDPointCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointCorrector_H
#define twoDPointCorrector_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class twoDPointCorrector Declaration
\*---------------------------------------------------------------------------*/

class twoDPointCorrector
{
    // Private data

        //- Reference to moving mesh
        const polyMesh& mesh_;

        //- 2-D plane unit normal
        mutable vector* planeNormalPtr_;

        //- Indices of edges normal to plane
        mutable labelList* normalEdgeIndicesPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        twoDPointCorrector(const twoDPointCorrector&);

        //- Disallow default bitwise assignment
        void operator=(const twoDPointCorrector&);


        //- Calculate addressing
        void calcAddressing() const;

        //- Clear addressing
        void clearAddressing() const;


    // Static data members

        //- Edge orthogonality tolerance
        static const scalar edgeOrthogonalityTol;


public:

    // Constructors

        //- Construct from components
        twoDPointCorrector(const polyMesh& mesh);


    // Member Functions

        //- Return plane normal
        const vector& planeNormal() const;

        //- Return indices of normal edges.
        const labelList& normalEdgeIndices() const;

        //- Return direction normal to plane
        direction normalDir() const;

        //- Correct motion points
        void correctPoints(pointField& p) const;

        //- Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

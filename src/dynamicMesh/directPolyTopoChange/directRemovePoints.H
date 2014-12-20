/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    directRemovePoints

Description
    Removes selected points from mesh and updates faces using these
    points.

SourceFiles
    directRemovePoints.C

\*---------------------------------------------------------------------------*/

#ifndef directRemovePoints_H
#define directRemovePoints_H

#include "typeInfo.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class directPolyTopoChange;
class mapPolyMesh;
class face;

/*---------------------------------------------------------------------------*\
                           Class directRemovePoints Declaration
\*---------------------------------------------------------------------------*/

class directRemovePoints
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Change the vertices of the face whilst keeping everything else
        //  (patch, zone) the same.
        void modifyFace
        (
            const label faceI,
            const face&,
            directPolyTopoChange&
        ) const;


        //- Disallow default bitwise copy construct
        directRemovePoints(const directRemovePoints&);

        //- Disallow default bitwise assignment
        void operator=(const directRemovePoints&);

public:

    //- Runtime type information
    ClassName("directRemovePoints");


    // Constructors

        //- Construct from mesh
        directRemovePoints(const polyMesh& mesh);


    // Member Functions

        // Helper functions

            //- Mark in pointCanBeDeleted the points that can be deleted
            //  (parallel synchronised) and returns the global number of these
            //  points. (this number is the global number before synchronisation
            //  so might be off!)
            //  A point can be deleted if
            //  - it is not used by any edge.
            //  or
            //  - is not used by an internal edge
            //  - is used by only two boundary edges. (note that these two
            //    edges will always be boundary ones!)
            //  - these two edges are sufficiently in line (cos > minCos)
            //  - all processors agree that point can be deleted.
            label countPointUsage
            (
                const scalar minCos,
                boolList& pointCanBeDeleted
            ) const;

        // Topology changes

            //- Play commands into polyTopoChange to remove points. Gets
            //  boolList (output of countPointUsage) as input.
            //  Does no check for whether resulting face is legal.
            void setRefinement(const boolList&, directPolyTopoChange&);

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

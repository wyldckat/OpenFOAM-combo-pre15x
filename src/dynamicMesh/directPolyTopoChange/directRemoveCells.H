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
    directRemoveCells

Description
    Given list of cells to remove insert all the topology changes. Works in
    two passes:
    - get faces that will become boundary faces
    - given these faces and the patches they have to go into make the
      changes.

SourceFiles
    directRemoveCells.C

\*---------------------------------------------------------------------------*/

#ifndef directRemoveCells_H
#define directRemoveCells_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class directPolyTopoChange;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class directRemoveCells Declaration
\*---------------------------------------------------------------------------*/

class directRemoveCells
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether or not to synchronize parallel case.
        const bool syncPar_;


    // Private Member Functions

        //- Decrease count of elements of f
        static void uncount
        (
            const labelList& f,
            labelList& nUsage
        );


public:

    //- Runtime type information
    ClassName("directRemoveCells");


    // Constructors

        //- Construct from mesh. syncPar: do parallel synchronization.
        directRemoveCells(const polyMesh& mesh, const bool syncPar = true);


    // Member Functions

        // Topology changes

            //- Get labels of exposed faces. These are
            //  - internal faces that become boundary faces
            //  - coupled faces that become uncoupled (since one of the sides
            //    gets deleted)
            labelList getExposedFaces(const labelList& cellsToRemove) const;

            //- Play commands into polyTopoChange to remove cells.
            //  patchIDs is for every element in facesToExpose (see above) the
            //  patch it has to go into. This cannot be a coupled patch!
            void setRefinement
            (
                const labelList& cellsToRemove,
                const labelList& facesToExpose,
                const labelList& patchIDs,
                directPolyTopoChange&
            ) const;

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&)
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

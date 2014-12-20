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
    directRemoveCells

Description
    Given list of cells to remove insert all the topology changes.

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

        //- Patch to give to new external faces
        word patchName_;

public:

    //- Runtime type information
    ClassName("directRemoveCells");


    // Constructors

        //- Construct from mesh and patch to assign to newly exposed faces.
        directRemoveCells(const polyMesh& mesh, const word& patchName);


    // Member Functions

        // Topology changes

            //- Play commands into polyTopoChange to remove cells.
            void setRefinement(const labelList&, directPolyTopoChange&) const;

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&)
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

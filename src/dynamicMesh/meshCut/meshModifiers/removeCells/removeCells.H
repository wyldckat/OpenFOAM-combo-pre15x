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
    removeCells

Description
    Given list of cells to remove insert all the topology changes.

SourceFiles
    removeCells.C

\*---------------------------------------------------------------------------*/

#ifndef removeCells_H
#define removeCells_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class polyTopoChange;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class removeCells Declaration
\*---------------------------------------------------------------------------*/

class removeCells
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Patch to give to new external faces
        word patchName_;

public:

    //- Runtime type information
    ClassName("removeCells");


    // Constructors

        //- Construct from mesh and patch to assign to newly exposed faces.
        removeCells(const polyMesh& mesh, const word& patchName);


    // Member Functions

        // Topology changes

            //- Play commands into polyTopoChange to remove cells.
            void setRefinement(const labelList&, polyTopoChange&) const;

            //- Force recalculation of locally stored data on topological change
            void updateTopology(const mapPolyMesh&)
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    addPatchCellLayer

Description
    Adds layer of cells to outside of polyPatch. Call setRefinement with
    offset vector for every patch point.

SourceFiles
    addPatchCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef addPatchCellLayer_H
#define addPatchCellLayer_H

#include "labelList.H"
#include "typeInfo.H"
#include "pointField.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class polyTopoChange;

/*---------------------------------------------------------------------------*\
                           Class addPatchCellLayer Declaration
\*---------------------------------------------------------------------------*/

class addPatchCellLayer
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        addPatchCellLayer(const addPatchCellLayer&);

        //- Disallow default bitwise assignment
        void operator=(const addPatchCellLayer&);


public:

    //- Runtime type information
    ClassName("addPatchCellLayer");


    // Constructors

        //- Construct from mesh.
        addPatchCellLayer(const polyMesh& mesh);


    // Member Functions

        // Edit

            //- Play commands into polyTopoChange to create add points.
            void setRefinement
            (
                const indirectPrimitivePatch& pp,
                const pointField& points,
                polyTopoChange&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

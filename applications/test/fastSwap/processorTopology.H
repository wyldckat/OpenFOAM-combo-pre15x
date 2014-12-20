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
    processorTopology

Description
    Determines processor-processor connection. After instantiation contains
    on all processors the processor-processor connection table.

    *this[procI] gives the list of neighbouring processors.

SourceFiles
    processorTopology.C

\*---------------------------------------------------------------------------*/

#ifndef processorTopology_H
#define processorTopology_H

#include "labelList.H"
#include "polyPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class processorTopology Declaration
\*---------------------------------------------------------------------------*/

class processorTopology
:
    public labelListList
{
    // Private data

        //- Local map from neighbourproc to patchI. Different per processor!
        //  -1 or patchI for connection to procID
        labelList procPatchMap_;

        //- From cell (=processor) to face (processor-processor communication)
        mutable labelListList* cellFacesPtr_;
    
        //- From face (processor-processor communication) to cell (=processor)
        mutable labelListList* faceCellsPtr_;


    // Private Member Functions

        //- Return all neighbouring processors of this processor. Set
        //  procPatchMap_.
        labelList procNeighbours(const polyPatchList&);

        //- Calculate cell-face and face-cell addressing.
        void calcAddressing() const;

public:


    // Constructors

        //- Construct from boundaryMesh
        processorTopology(const polyPatchList& patches);


    // Destructor

        ~processorTopology();


    // Member Functions

        //- From processor to unique inter-processor 'face'.
        const labelListList& cellFaces() const
        {
            if (!cellFacesPtr_)
            {
                calcAddressing();
            }
            return *cellFacesPtr_;
        }


        //- From unique inter-processor 'face' to the two processors using it.
        const labelListList& faceCells() const
        {
            if (!faceCellsPtr_)
            {
                calcAddressing();
            }
            return *faceCellsPtr_;
        }


        //- from neighbour processor to index in boundaryMesh. Local information
        //  (so not same over all processors)
        const labelList& procPatchMap() const
        {
            return procPatchMap_;
        }


        //- Get label of 'face' for connection to neighbour procID on current
        //  processor.
        label getFace(const label nbrProcID) const;

        //- Get label of 'face' for connection to neighbour procID 
        label getFace(const label myProcID, const label nbrProcID) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

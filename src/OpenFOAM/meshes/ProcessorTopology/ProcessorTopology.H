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
    ProcessorTopology

Description
    Determines processor-processor connection. After instantiation contains
    on all processors the processor-processor connection table.

    *this[procI] gives the list of neighbouring processors.

SourceFiles
    ProcessorTopology.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorTopology_H
#define ProcessorTopology_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ProcessorTopology Declaration
\*---------------------------------------------------------------------------*/

//- Struct to hold the patch index and the initialisation flag
//  for the patch schedule
struct labool
{
    label patch;
    bool init;
    bool bufferedTransfer;

    friend bool operator!=(const labool& lb1, const labool& lb2)
    {
        return true;
    }

    friend Ostream& operator<<(Ostream& os, const labool& lb)
    {
        os  << lb.patch << token::SPACE
            << lb.init << token::SPACE
            << lb.bufferedTransfer;
        return os;
    }
};

typedef List<labool> patchScheduleList;


template<class Patch, class ProcPatch>
class ProcessorTopology
:
    public labelListList
{

private:

    // Private data

        //- Local map from neighbour proc to patchI. Different per processor!
        //  -1 or patchI for connection to procID
        labelList procPatchMap_;

        //- From cell (=processor) to face (processor-processor communication)
        labelListList cellFaces_;
    
        //- From face (processor-processor communication) to cell (=processor)
        labelListList faceCells_;

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        patchScheduleList patchSchedule_;


    // Private Member Functions

        //- Return all neighbouring processors of this processor. Set
        //  procPatchMap_.
        labelList procNeighbours(const PtrList<Patch>&);

        //- Calculate cell-face and face-cell addressing.
        void calcAddressing();


public:

    // Constructors

        //- Construct from boundaryMesh
        ProcessorTopology(const PtrList<Patch>& patches);


    // Member Functions

        //- From processor to unique inter-processor 'face'.
        const labelListList& cellFaces() const
        {
            return cellFaces_;
        }

        //- From unique inter-processor 'face' to the two processors using it.
        const labelListList& faceCells() const
        {
            return faceCells_;
        }

        //- from neighbour processor to index in boundaryMesh. Local information
        //  (so not same over all processors)
        const labelList& procPatchMap() const
        {
            return procPatchMap_;
        }

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        const patchScheduleList& patchSchedule() const
        {
            return patchSchedule_;
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

#ifdef NoRepository
#   include "ProcessorTopology.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

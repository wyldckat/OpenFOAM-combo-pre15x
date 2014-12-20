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
    regionSplit

Description
    This class separates the mesh into distint unconnected regions,
    each of which is then given a label.

SourceFiles
    regionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit_H
#define regionSplit_H

#include "polyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class regionInfo;

/*---------------------------------------------------------------------------*\
                           Class regionSplit Declaration
\*---------------------------------------------------------------------------*/

class regionSplit
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Number of regions
        mutable label nRegions_;

        //- Cell-to-region list
        mutable labelList*  cellToRegionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionSplit(const regionSplit&);

        //- Disallow default bitwise assignment
        void operator=(const regionSplit&);


        //- Calculate region split
        void calcRegionSplit() const;

        //- Given a seed cell label, fill mask with 1 for contiguous
        //  region around it
        void fillSeedMask
        (
            labelList& region,
            const label seedCellID,
            const label mark
        ) const;


public:

    // Constructors

        //- Construct from mesh
        regionSplit(const polyMesh& m)
        :
            mesh_(m),
            nRegions_(-1),
            cellToRegionPtr_(NULL)
        {}

    // Destructor

        ~regionSplit()
        {
            deleteDemandDrivenData(cellToRegionPtr_);
        }


    // Member Functions

        //- Return number of regions
        label nRegions() const;

        //- Return region split
        const labelList& cellToRegion() const;

        //- Return mask covering contiguous cells from the given seed cell
        labelList seedMask(const label seedCellID) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    channelIndex

Description
    does indexing for a standard (Christer-style) channel.
Assumes that the blocks are aranged in the y direction.

SourceFiles
    channelIndex.C
    channelIndexIO.C

\*---------------------------------------------------------------------------*/

#ifndef channelIndex_H
#define channelIndex_H

#include "fvCFD.H"


/*---------------------------------------------------------------------------*\
                           Class channelIndex Declaration
\*---------------------------------------------------------------------------*/

class channelIndex
{

    // Private data

    IOdictionary indexingDict;

    const label Nx;
    const labelList Ny;
    const label Nz;

    labelList cumNy;

    label Nlayers;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        channelIndex(const channelIndex&);
        void operator=(const channelIndex&);


public:

    // Constructors

        //- from components
        channelIndex(const fvMesh& m);

    // Destructor

        ~channelIndex();


    // Member Functions

        // Access

            //- return number of layers
            label nLayers() const
            {
                return Nlayers;
            }

            //- number of cells in X direction
            label nx() const
            {
                return Nx;
            }

            //- number of cells in Z direction
            label nz() const
            {
                return Nz;
            }

            //- collapse a field to a line
            scalarField collapse
            (
                const volScalarField& vsf
            ) const;

            //- return the field of Y locations from the cell centres
            scalarField y
            (
                const volVectorField& cellCentres
            ) const;


    // Member Operators

        label operator()(const label, const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

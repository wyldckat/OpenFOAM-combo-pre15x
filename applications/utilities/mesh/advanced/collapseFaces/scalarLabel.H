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
    scalarLabel

Description
    Holds scalar and label and nothing more.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef scalarLabel_H
#define scalarLabel_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarLabel Declaration
\*---------------------------------------------------------------------------*/

class scalarLabel
{
    // Private data

        scalar val_;

        label index_;

public:

    // Constructors

        //- Construct empty
        scalarLabel()
        :
            val_(-GREAT),
            index_(-1)
        {}

        //- Construct from components
        scalarLabel(const scalar val, const label index)
        :
            val_(val),
            index_(index)
        {}

    // Member Functions

        // Access

            scalar val() const
            {
                return val_;
            }

            label index() const
            {
                return index_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

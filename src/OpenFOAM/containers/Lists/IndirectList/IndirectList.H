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
    IndirectList

Description

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef IndirectList_H
#define IndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T> class IndirectList;

template<class T> Ostream& operator<<(Ostream&, const IndirectList<T>&);


/*---------------------------------------------------------------------------*\
                           Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
{
    // Private data

        const UList<T>& completeList_;
        List<label> addressing_;


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline IndirectList(const UList<T>&, const List<label>&);


    // Member Functions

        // Access

            inline label size() const;
            inline const UList<T>& completeList() const;
            inline const List<label>& addressing() const;

        // Edit

            //- Reset addressing
            void resetAddressing(const List<label>& addr);


        // Member Operators

            //- Return the addressed elements as a List
            inline List<T> operator()() const;

            //- Return const access to an element
            inline const T& operator[](const label) const;


    // Ostream operator

        // Write IndirectList to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const IndirectList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

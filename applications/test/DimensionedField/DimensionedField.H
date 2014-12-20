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
    DimensionedField

Description

SourceFiles
    DimensionedFieldI.H
    DimensionedField.C
    DimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedField_H
#define DimensionedField_H

#include "regIOobject.H"
#include "Field.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //
 
template<class Type, class GeoMesh> class DimensionedField;
template<class Type, class GeoMesh> Ostream& operator<<
(
    Ostream&,
    const DimensionedField<Type, GeoMesh>&
);


/*---------------------------------------------------------------------------*\
                           Class DimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedField
:
    public regIOobject,
    public Field<Type>
{
    // Private data

        //- Dimension set for this field
        dimensionSet dimensions_;


public:

    //- Runtime type information
    TypeName("DimensionedField");


    // Constructors

        //- Construct from components
        DimensionedField
        (
            const IOobject&,
            const dimensionSet&,
            const Field<Type>&
        );

        //- Construct from Istream
        DimensionedField
        (
            const IOobject&,
            const word& fieldDictEntry="fieldValues"
        );

        //- Construct as copy
        DimensionedField(const DimensionedField&);


    // Destructor

        ~DimensionedField();


    // Member Functions

        // Access

            //- Return dimensions
            inline dimensionSet& dimensions();

            //- Return dimensions
            inline const dimensionSet& dimensions() const;

            //- Calculate and return arithmetic average
            dimensioned<Type> average() const;


        // Check

        // Edit

        // Write

            bool writeData(Ostream&) const;


    // Member Operators

        void operator=(const DimensionedField<Type, GeoMesh>&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream&,
            const DimensionedField<Type, GeoMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldI.H"

#ifdef NoRepository
#   include "DimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

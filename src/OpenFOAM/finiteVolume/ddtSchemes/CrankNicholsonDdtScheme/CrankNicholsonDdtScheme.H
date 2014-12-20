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
    CrankNicholsonDdtScheme

Description
    Second-oder CrankNicholson implicit ddt using the current and
    two previous time-step as well as the previous time-step ddt.

SourceFiles
    CrankNicholsonDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicholsonDdtScheme_H
#define CrankNicholsonDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class CrankNicholsonDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CrankNicholsonDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        CrankNicholsonDdtScheme(const CrankNicholsonDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const CrankNicholsonDdtScheme&);

        template<class GeoMesh>
        GeometricField<Type, fvPatchField, GeoMesh>& ddt0_
        (
            const word& name,
            const dimensionSet& dims
        );

        template<class GeoMesh>
        bool evaluate
        (
            const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
        ) const;

        template<class GeoMesh>
        dimensionedScalar rDeltaT_
        (
            const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
        ) const;

        template<class GeoMesh>
        dimensionedScalar rDeltaT0_
        (
            const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
        ) const;


public:

    //- Runtime type information
    TypeName("CrankNicholson");


    // Constructors

        //- Construct from mesh
        CrankNicholsonDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        CrankNicholsonDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt0
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt0
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, surfaceMesh> > fvcDdt0
        (
            const GeometricField<Type, fvPatchField, surfaceMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt0
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CrankNicholsonDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

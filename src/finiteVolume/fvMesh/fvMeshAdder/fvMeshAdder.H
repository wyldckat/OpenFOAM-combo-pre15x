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
    fvMeshAdder

Description
    Adds two fvMeshes without using any polyMesh morphing.
    Uses fvMeshAdder.

SourceFiles
    fvMeshAdder.C
    fvMeshAdderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshAdder_H
#define fvMeshAdder_H

#include "polyMeshAdder.H"
#include "fvPatchFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class IOobject;
class faceCoupleInfo;
class IOobjectList;
class fvMesh;
class volMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMeshAdder Declaration
\*---------------------------------------------------------------------------*/

class fvMeshAdder
:
    public polyMeshAdder
{

private:
    // Private data


    // Private Member Functions

        //- Map from old to new according to map. Handles map = -1.
        template<class Type>
        static void map
        (
            const Field<Type>&,
            const labelList& map,
            Field<Type>&
        );

        //- Map from old to new fvPatchField. Handles map = -1.
        template<class Type>
        static void mapPatchField
        (
            const Field<Type>& oldFld,
            const label oldStart,
            const labelList& oldToNew,
            fvPatchField<Type>& fld
        );

        //- Update single volField.
        template<class Type>
        void MapVolField
        (
            const labelList& oldPatchStarts,
            GeometricField<Type, fvPatchField, volMesh>& fld,
            const GeometricField<Type, fvPatchField, volMesh>& fldToAdd
        ) const;

        template<class Type>
        void MapVolFields
        (
            const labelList& oldPatchStarts,
            const fvMesh& mesh,
            const fvMesh& meshToAdd
        ) const;

public:

    // Constructors

        //- Construct null
        fvMeshAdder();


    // Member Functions

        //- Inplace add mesh to fvMesh. Maps all stored fields.
        void add
        (
            fvMesh& mesh0,
            const fvMesh& mesh1,
            const faceCoupleInfo& coupleInfo
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMeshAdderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Description
    generic Geometric Type field class

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return reference to mesh
template<class Type, template<class> class PatchField, class GeoMesh>
inline const typename GeometricField<Type, PatchField, GeoMesh>::Mesh&
GeometricField<Type, PatchField, GeoMesh>::mesh() const
{
    return mesh_;
}


// Return dimensions
template<class Type, template<class> class PatchField, class GeoMesh>
inline dimensionSet&
GeometricField<Type, PatchField, GeoMesh>::dimensions()
{
    return dimensions_;
}

// Return dimensions
template<class Type, template<class> class PatchField, class GeoMesh>
inline const dimensionSet&
GeometricField<Type, PatchField, GeoMesh>::dimensions() const
{
    return dimensions_;
}

// Return internal field
template<class Type, template<class> class PatchField, class GeoMesh>
inline const Field<Type>&
GeometricField<Type, PatchField, GeoMesh>::internalField() const
{
    return *this;
}

// Return reference to GeometricBoundaryField for const field
template<class Type, template<class> class PatchField, class GeoMesh>
inline const typename GeometricField<Type, PatchField, GeoMesh>::
GeometricBoundaryField&
GeometricField<Type, PatchField, GeoMesh>::boundaryField() const
{
    return boundaryField_;
}


// Return reference to GeometricBoundaryField for const field
template<class Type, template<class> class PatchField, class GeoMesh>
inline label GeometricField<Type, PatchField, GeoMesh>::timeIndex() const
{
    return timeIndex_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

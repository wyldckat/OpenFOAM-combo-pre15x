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

\*---------------------------------------------------------------------------*/


#include "TractionPointPatchVectorField.H"
#include "constraints.H"
#include "tetFemMatrix.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch
>
void TractionPointPatchVectorField<PatchField, PointPatch>
::checkFieldSize() const
{
    if
    (
        this->size() != this->patchMesh().size()
     || traction_.size() != this->patchMesh().size()
     || pressure_.size() != this->patchMesh().size()
    )
    {
        FatalErrorIn
        (
            "void TractionPointPatchVectorField::checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << this->size()
            << " traction size: " << traction_.size()
            << " pressure size: " << pressure_.size()
            << " patch size: " << this->patchMesh().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const vectorField& iF
)
:
    PatchField<vector>(p, iF),
    traction_(p.size()),
    pressure_(p.size())
{}


template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const vectorField& iF,
    const vectorField& v,
    const scalarField& vf
)
:
    PatchField<vector>(p, iF),
    traction_(v),
    pressure_(vf)
{
    checkFieldSize();
}


template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const vectorField& iF,
    Istream& is
)
:
    PatchField<vector>(p, iF),
    traction_(is),
    pressure_(is)
{
    checkFieldSize();
}


template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    PatchField<vector>(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    this->updateBoundaryField();
}


template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const TractionPointPatchVectorField<PatchField, PointPatch>& ptf,
    const PointPatch& p,
    const vectorField& iF,
    const PointPatchFieldMapper& mapper
)
:
    PatchField<vector>(p, iF),
    traction_(ptf.traction_, (const FieldMapper&)mapper),
    pressure_(ptf.pressure_, (const FieldMapper&)mapper)
{}


template
<
    template<class> class PatchField,
    class PointPatch
>
TractionPointPatchVectorField<PatchField, PointPatch>::
TractionPointPatchVectorField
(
    const TractionPointPatchVectorField<PatchField, PointPatch>& ptf,
    const vectorField& iF
)
:
    PatchField<vector>(ptf, iF),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Grab the values using rmap
template
<
    template<class> class PatchField,
    class PointPatch
>
void TractionPointPatchVectorField<PatchField, PointPatch>::rmap
(
    const PointPatchField<PatchField, PointPatch, vector>& ptf,
    const labelList& addr
)
{
    const TractionPointPatchVectorField& mptf =
        refCast<const TractionPointPatchVectorField>(ptf);

    traction_.rmap(mptf.traction_, addr);
    pressure_.rmap(mptf.pressure_, addr);
}


// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class PointPatch
>
void TractionPointPatchVectorField<PatchField, PointPatch>::
addBoundarySourceDiag
(
    tetFemMatrix<vector>& matrix
) const
{
    const vectorField& points = this->patchMesh().localPoints();

    const labelList& meshPoints = this->patchMesh().meshPoints();
    const vectorField& pointNormals = this->patchMesh().pointNormals();
    const label nFaces = this->patchMesh().nFaces();

    vectorField& source = matrix.source();

    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        List<triFace> faceTriangles =
            this->patchMesh().faceTriangles(faceI);

        forAll (faceTriangles, triI)
        {
            const triFace& tri = faceTriangles[triI];

            source[meshPoints[tri[0]]] +=
                tri.mag(points)*
                (
                    traction_[tri[0]]/6.0
                  - pressure_[tri[0]]*pointNormals[tri[0]]/6.0
                  + traction_[tri[1]]/12.0
                  - pressure_[tri[1]]*pointNormals[tri[1]]/12.0
                  + traction_[tri[2]]/12.0
                  - pressure_[tri[2]]*pointNormals[tri[2]]/12.0
                );

            source[meshPoints[tri[1]]] +=
                tri.mag(points)*
                (
                    traction_[tri[1]]/6.0
                  - pressure_[tri[1]]*pointNormals[tri[1]]/6.0
                  + traction_[tri[2]]/12.0
                  - pressure_[tri[2]]*pointNormals[tri[2]]/12.0
                  + traction_[tri[0]]/12.0
                  - pressure_[tri[0]]*pointNormals[tri[0]]/12.0
                );

            source[meshPoints[tri[2]]] +=
                tri.mag(points)*
                (
                    traction_[tri[2]]/6.0
                  - pressure_[tri[2]]*pointNormals[tri[2]]/6.0
                  + traction_[tri[0]]/12.0
                  - pressure_[tri[0]]*pointNormals[tri[0]]/12.0
                  + traction_[tri[1]]/12.0
                  - pressure_[tri[1]]*pointNormals[tri[1]]/12.0
                );
        }
    }
}


// Write
template
<
    template<class> class PatchField,
    class PointPatch
>
void TractionPointPatchVectorField<PatchField, PointPatch>::
write(Ostream& os) const
{
    PatchField<vector>::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

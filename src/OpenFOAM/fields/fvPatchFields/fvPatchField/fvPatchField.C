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

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void fvPatchField<Type>::checkPatch() const
{
    if (this->size() != patch_.size())
    {
        FatalErrorIn("fvPatchField<Type>::checkPatch() const")
            << "number of fvPatchField entries (" << this->size()
            << ") != number of patch faces (" << patch_.size() << ')'
            << " for patch: " << patch_.name() << " index " << patch_.index()
            << abort(FatalError);
    }
}


template<class Type>
void fvPatchField<Type>::checkInternalField() const
{
    if (&internalField_ && internalField_.size())
    {
        label iFs = internalField_.size();

        const fvMesh& mesh = patch_.boundaryMesh().mesh();
        label nPoints = mesh.nPoints();
        label nInternalFaces = mesh.nInternalFaces();
        label nCells = mesh.nCells();

        if (iFs != nPoints && iFs != nInternalFaces && iFs != nCells)
        {
            FatalErrorIn("fvPatchField<Type>::checkInternalField() const")
                << "internal field size " << iFs << " does not match" << nl
                << "    the number of points " << nPoints << nl
                << "    the number of internal faces " << nInternalFaces << nl
                << "    or the number of cells " << nCells
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false)
{
    checkInternalField();
}


template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false)
{
    checkInternalField();
}


template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(ptf, mapper),
    patch_(p),
    internalField_(iF),
    updated_(false)
{
    checkInternalField();
}


template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(pTraits<Type>::zero);
    }
}


template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf
)
:
    lduCoupledInterface(),
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false)
{}


template<class Type>
fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false)
{
    checkInternalField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const objectRegistry& fvPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh();
}


template<class Type>
void fvPatchField<Type>::check(const fvPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorIn("PatchField<Type>::check(const fvPatchField<Type>&)")
            << "different patches for fvPatchField<Type>s"
            << abort(FatalError);
    }
}


// Does this patchField correspond to a volTypeField
template<class Type>
bool fvPatchField<Type>::isVolField() const
{
    return
        internalField().size()
     == patch().boundaryMesh().mesh().nCells();
}

// Does this patchField correspond to a volTypeField
template<class Type>
void fvPatchField<Type>::checkVolField() const
{
    if (!isVolField())
    {
        FatalErrorIn("fvPatchField<Type>::checkVolField() const")
            << "This " << type() << " patchField"
            << " is not part of a volTypeField which may cause "
               "undefined behaviour from the evaluate and other functions"
            << abort(FatalError);
    }
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > fvPatchField<Type>::snGrad() const
{
    return (*this - patchInternalField())*patch_.deltaCoeffs();
}


// Return internal field next to patch as patch field
template<class Type>
template<class Type2>
tmp<Field<Type2> > fvPatchField<Type>::patchInternalField
(
    const Field<Type2>& iField
) const
{
    tmp<Field<Type2> > tpif(new Field<Type2>(patch_.size()));
    Field<Type2>& pif = tpif();

    const labelList::subList FaceCells = patch_.faceCells();

    forAll(patch_, faceI)
    {
        pif[faceI] = iField[FaceCells[faceI]];
    }

    return tpif;
}


// Return internal field next to patch as patch field
template<class Type>
tmp<Field<Type> > fvPatchField<Type>::patchInternalField() const
{
    return patchInternalField(internalField_);
}


//- Return neighbour coupled given internal cell data
template<class Type>
tmp<Field<Type> > fvPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    return patchInternalField(iField);
}


//- Return neighbour coupled given internal cell data
template<class Type>
tmp<labelField> fvPatchField<Type>::nbrColour
(
    const labelField&
) const
{
    // Dummy return to avoid unnecessary oprations on an uncoupled interface
    // 
    return labelField(0);
}


//- Return patchField of the values on the patch or on the
//  opposite patch
template<class Type>
tmp<Field<Type> > fvPatchField<Type>::patchNeighbourField() const
{
    return *this;
}


template<class Type>
template<class GeometricField, class Type2>
const fvPatchField<Type2>& fvPatchField<Type>::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[patch_.index()];
}


template<class Type>
template<class GeometricField, class Type2>
const fvPatchField<Type2>& fvPatchField<Type>::lookupPatchField
(
    const word& name,
    const GeometricField*,
    const Type2*
) const
{
    return patchField<GeometricField, Type2>
    (
        db().objectRegistry::lookupObject<GeometricField>(name)
    );
}


// Map from self
template<class Type>
void fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (m.resizeOnly())
    {
        setSize(m.size());
    }
    else
    {
        Field<Type>::autoMap(m);
    }
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


// Evaluate the patch field, sets Updated to false
template<class Type>
void fvPatchField<Type>::evaluate()
{
    if (!updated_)
    {
        updateCoeffs();
    }
    
    updated_ = false;
}


// Write
template<class Type>
void fvPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator*=(const fvPatchField<scalar>& ptf)"
        )   << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator/=(const fvPatchField<scalar>& ptf)"
        )   << "    incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void fvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void fvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void fvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const fvPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fvPatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newFvPatchField.C"

// ************************************************************************* //

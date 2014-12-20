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

#include "timeVaryingUniformFixedValueFvPatchField.H"
#include "graph.H"
#include "IFstream.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataPtr_(NULL)
{}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper&
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL)
{}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand()),
    timeDataPtr_(NULL)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        updateCoeffs();
    }
}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void timeVaryingUniformFixedValueFvPatchField<Type>::checkTable()
{
    if (!timeDataPtr_.valid())
    {
        timeDataPtr_.reset
        (
            new graph("title", "x", "y", IFstream(timeDataFileName_)())
        );
    }

    if (this->db().time().value() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "timeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << this->db().time().value()
            << ") is less than the minimum in the data table ("
            << min(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

    if (this->db().time().value() > max(timeDataPtr_().x()))
    {
        WarningIn
        (
            "timeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << this->db().time().value()
            << ") is greater than the maximum in the data table ("
            << max(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}


// Write
template<class Type>
void timeVaryingUniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

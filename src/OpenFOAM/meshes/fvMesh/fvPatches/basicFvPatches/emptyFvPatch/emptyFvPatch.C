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

#include "emptyFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const vectorField emptyFvPatch::emptyNf_(0);

// Patch name
defineTypeNameAndDebug(emptyFvPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable(fvPatch, emptyFvPatch, polyPatch);

// Over-riding the face normals return from the underlying patch
// This is the only piece of info used out of the underlying primitivePatch
// I choose to store it there because it is used in primitive patch operations
// and it should not be duplicated as before.  However, to ensure everything
// in the empty patch is sized to zero, we shall here return a reference to
// a zero-sized field (it does not matter what the field is
// 
const vectorField& emptyFvPatch::nf() const
{
    return emptyNf_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

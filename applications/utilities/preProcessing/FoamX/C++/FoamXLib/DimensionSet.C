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

#include "DimensionSet.H"

// * * * * * * * * * * * * * * * Global Operators * * * * * * * * * * * * * //

void FoamX::operator==
(
    FoamXServer::DimensionSet& fxDs,
    const Foam::dimensionSet& ds
)
{
    fxDs.mass = ds[Foam::dimensionSet::MASS];
    fxDs.length = ds[Foam::dimensionSet::LENGTH];
    fxDs.time = ds[Foam::dimensionSet::TIME];
    fxDs.temperature = ds[Foam::dimensionSet::TEMPERATURE];
    fxDs.moles = ds[Foam::dimensionSet::MOLES];
    fxDs.current = ds[Foam::dimensionSet::CURRENT];
    fxDs.luminousIntensity = ds[Foam::dimensionSet::LUMINOUS_INTENSITY];
}


// ************************************************************************* //

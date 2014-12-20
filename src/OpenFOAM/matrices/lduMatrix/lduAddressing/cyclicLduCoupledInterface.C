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

#include "cyclicLduCoupledInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicLduCoupledInterface, 0);


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicLduCoupledInterface::~cyclicLduCoupledInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::cyclicLduCoupledInterface::transformCyclicCoupleField
(
    scalarField& pnf,
    const direction cmpt
) const
{
    if (doTransform())
    {
        label sizeby2 = pnf.size()/2;

        scalar forwardScale = 
            pow(diag(forwardT()[0]).component(cmpt), rank());

        scalar reverseScale =
            pow(diag(reverseT()[0]).component(cmpt), rank());

        for (label facei=0; facei<sizeby2; facei++)
        {
            pnf[facei] *= forwardScale;
            pnf[facei + sizeby2] *= reverseScale;
        }
    }
}

// ************************************************************************* //

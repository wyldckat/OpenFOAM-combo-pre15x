/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "channelIndex.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
channelIndex::channelIndex(const fvMesh& m)
:
    indexingDict
    (
        IOobject
        (
            "postChannelDict",
            m.time().constant(),
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Nx(readLabel(indexingDict.lookup("Nx"))),
    Ny(indexingDict.lookup("Ny")),
    Nz(readLabel(indexingDict.lookup("Nz"))),
    cumNy(Ny.size()),
    Nlayers(Ny[0])
{
    // initialise the layers
    cumNy[0] = Ny[0];

    for (label j = 1; j<Ny.size(); j++)
    {
        Nlayers += Ny[j];
        cumNy[j] = Ny[j]+cumNy[j-1];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

channelIndex::~channelIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField channelIndex::collapse
(
    const volScalarField& vsf
) const
{
    scalarField cs(nLayers(), 0.0);

    forAll(cs, j)
    {
        // sweep over all cells in this layer
        for (label i=0; i<nx(); i++)
        {
            for (label k=0; k<nz(); k++)
            {
                cs[j] += vsf[operator()(i,j,k)];
            }
        }

        // and divide by the number of cells in the layer
        cs[j] /= scalar(nx()*nz());
    }

    return cs;
}


scalarField channelIndex::y
(
    const volVectorField& cellCentres
) const
{
    scalarField Y(nLayers());

    label j;
    for (j=0; j<nLayers(); j++)
    {
        Y[j] = cellCentres[operator()(0,j,0)].y();
    }

    return Y;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

label channelIndex::operator()
(
    const label Jx,
    const label Jy,
    const label Jz
) const
{
    label index(0);

    // count up `full' layers in the mesh
    label j(0);
    label tmpJy(Jy);

    while (Jy>=cumNy[j])
    {
        index += Nx*Ny[j]*Nz;
        tmpJy -= Ny[j];
        j++;
    }

    index += Jx + Nx*tmpJy + Nx*Ny[j]*Jz;

    return index;
}


// ************************************************************************* //

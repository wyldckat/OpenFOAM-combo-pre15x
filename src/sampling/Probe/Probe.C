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

Description
    Probes a volTypeField on call to write() and writes the values
    into a file.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "Probe.H"
#include "fvMesh.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class volTypeField>
Probe<volTypeField>::Probe(const dictionary& dict, const volTypeField& field)
:
    field_(field),
    probeLocations_(dict.lookup("probeLocations")),
    cellList_(probeLocations_.size()),
    probes_(false),
    probeFilePtr_(NULL)
{
    const fvMesh& mesh = field_.mesh();

    label probesFoundi = 0;
    forAll(probeLocations_, probei)
    {
        label probeCelli = mesh.findCell(probeLocations_[probei]);

        if (probeCelli >= 0)
        {
            cellList_[probesFoundi] = probeCelli;
            probeLocations_[probesFoundi] = probeLocations_[probei];
            probesFoundi++;
        }
    }

    probeLocations_.setSize(probesFoundi);
    cellList_.setSize(probesFoundi);

    if (cellList_.size() > 0)
    {
        probes_ = true;
    }

    if (probes_)
    {
        mkDir(field_.time().path()/"probes"/field_.time().timeName());

        probeFilePtr_ = new OFstream
        (
            field_.time().path()
            /"probes"
            /field_.time().timeName()
            /field_.name()
        );

        *probeFilePtr_ << "Time" << tab;

        forAll(cellList_, probei)
        {
            *probeFilePtr_ << probeLocations_[probei] << tab;
        }
        *probeFilePtr_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class volTypeField>
void Probe<volTypeField>::write()
{
    if (probes_)
    {
        *probeFilePtr_ << field_.time().value() << tab;

        forAll(cellList_, probei)
        {
            *probeFilePtr_ << field_[cellList_[probei]] << tab;
        }
        *probeFilePtr_ << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

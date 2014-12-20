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
    write function for regIOobjects

\*---------------------------------------------------------------------------*/

#include "regIOobject.H"
#include "Time.H"
#include "OSspecific.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//  Write object to Ostream according to read/write options
bool regIOobject::write
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (!good())
    {
        SeriousErrorIn("regIOobject::write()")
            << "bad object " << name()
            << endl;

        return false;
    }

    if (!instance().size())
    {
        SeriousErrorIn("regIOobject::write()")
            << "instance undefined for object " << name()
            << endl;

        return false;
    }

    if
    (
        instance() != time().timeName()
     && instance() != time().constant()
     && instance() != "constant"
     && instance() != time().system()
    )
    {
        const_cast<regIOobject&>(*this).instance() = time().timeName();
    }

    mkDir(path());

    if (OFstream::debug)
    {
        Info<< "regIOobject::write() : "
            << "writing file " << objectPath();
    }

    // Try opening a Ostream on object
    OFstream os(objectPath(), fmt, ver, cmp);

    // If this has failed, return (leave error handling to  Ostream class)
    if (!os.good())
    {
        return false;
    }

    if (!writeHeader(os))
    {
        return false;
    }

    // Write the data to the Ostream
    if (!writeData(os))
    {
        return false;
    }

    os  << "\n\n// ************************************************************************* //"
        << endl;

    if (OFstream::debug)
    {
        Info<< " .... written" << endl;
    }

    lastModified_ = lastModified(objectPath());

    return os.good();
}


//  Write using setting from DB
bool regIOobject::write() const
{
    return write
    (
        time().writeFormat(),
        IOstream::currentVersion,
        time().writeCompression()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

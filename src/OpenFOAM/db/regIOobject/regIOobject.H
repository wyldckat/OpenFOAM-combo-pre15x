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

Class
    regIOobject

Description
    regIOobject is an abstract class derived from IOobject to handle
    automatic object registration with the objectRegistry.

SourceFiles
    regIOobject.C
    regIOobjectRead.C
    regIOobjectWrite.C

\*---------------------------------------------------------------------------*/

#ifndef regIOobject_H
#define regIOobject_H

#include "IOobject.H"
#include "typeInfo.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regIOobject Declaration
\*---------------------------------------------------------------------------*/

class regIOobject
:
    public IOobject
{

private:

    // Private data

        //- Is this object registered with the registry
        bool registered_;

        //- Is this object owned by the registry
        bool registries_;

        //- Time of last modification
        mutable time_t lastModified_;

        //- Has the objects file been updated
        bool updated_;

        //- Istream for reading
        Istream* isPtr_;


    // Private member functions

        //- Return Istream
        Istream& readStream();

        //- Dissallow assignment
        void operator=(const regIOobject&);


public:

    // Static data

        //- Runtime type information
        TypeName("regIOobject");

        static int fileModificationSkew;


    // Constructors

        //- Construct from IOobject
        regIOobject(const IOobject&);

        //- Construct as copy
        regIOobject(const regIOobject&);

        //- Construct as copy, and transfering registry registration to copy
        //  if registerCopy is true
        regIOobject(const regIOobject&, bool registerCopy);


    // Destructor

        virtual ~regIOobject();


    // Member functions

        // Registration

            //- Register object with registry
            void checkIn();

            //- Check-out object from registry
            void checkOut();

            //- Is this object owned by the registry
            bool registries() const
            {
                return registries_;
            }

            //- Transfer ownership of this object to it's registry
            void transferToRegistry()
            {
                registries_ = true;
            }


        // Update

            //- Reset the updated flag
            virtual void resetUpdate();

            //- Update object if it hasn't already been updated
            virtual bool update();

            //- Has the objects file been updated
            bool updated() const
            {
                return updated_;
            }


        // Edit

            //- Rename
            virtual void rename(const word& newName);


        // Reading

            //- Return Istream and check object type against that given
            Istream& readStream(const word&);

            //- Close Istream
            void close();

            //- Virtual readData function.
            //  Must be defined in derived types for which
            //  re-reading is required
            virtual bool readData(Istream&);

            //- Read object
            virtual bool read();

            //- Read object if modified
            virtual bool readIfModified();


        // Writing

            //- Pure virtual writaData function.
            //  Must be defined in derived types
            virtual bool writeData(Ostream&) const = 0;

            //- Write using given format, version and compression
            virtual bool write
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;

            //- Write using setting from DB
            bool write() const;


    // Member operators

        void operator=(const IOobject&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

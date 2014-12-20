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
    Time

Description
    Class to control time during FOAM simulations which is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    findInstance.C

\*---------------------------------------------------------------------------*/

#ifndef Time_H
#define Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "NamedEnum.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
    // Private data

        //- The controlDict
        IOdictionary controlDict_;


public:

        //- Write control options
        enum writeControls
        {
            wcTimeStep,
            wcRunTime,
            wcAdjustableRunTime,
            wcClockTime,
            wcCpuTime
        };

        //- Stop-run control options
        enum stopAtControls
        {
            saEndTime,
            saNoWriteNow,
            saWriteNow,
            saNextWrite
        };

        //- Suported time directory name formats
        enum fmtflags
        {
            general = 0,
            fixed = ios_base::fixed,
            scientific = ios_base::scientific
        };


protected:

    // Protected data

        scalar startTime_;
        scalar endTime_;

        static const NamedEnum<stopAtControls, 4> stopAtControlNames_;
        stopAtControls stopAt_;

        static const NamedEnum<writeControls, 5> writeControlNames_;
        writeControls writeControl_;

        scalar writeInterval_;
        label cycleWrite_;

        //- Time directory name format
        static fmtflags format_;

        //- Time directory name precision
        static int precision_;

        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- Default write option
        IOstream::streamFormat writeFormat_;

        //- Default output file format version number
        IOstream::versionNumber writeVersion_;

        //- Default output compression
        IOstream::compressionType writeCompression_;

        //- Default graph format
        word graphFormat_;

        //- Is runtim modification of dictionaries allowed
        Switch runTimeModifiable_;


public:

    TypeName("time");

    //- Return the default control dictionary name
    static word controlDictName;


    // Constructors

        //- Construct given name, rootPath and casePath
        Time
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


    // Destructor

        //- Virtual destructor
        virtual ~Time();


    // Member functions

        // Database functions

            //- Return root path
            const fileName& rootPath() const
            {
                return TimePaths::rootPath();
            }

            //- Return case name
            const fileName& caseName() const
            {
                return TimePaths::caseName();
            }

            //- Return path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            const dictionary& controlDict() const
            {
                return controlDict_;
            }

            virtual const fileName& dbDir() const
            {
                return fileName::null;
            }

            //- Return current time path
            fileName timePath() const
            {
                return path()/timeName();
            }

            //- Default write format
            IOstream::streamFormat writeFormat() const
            {
                return writeFormat_;
            }

            //- Default write version number
            IOstream::versionNumber writeVersion() const
            {
                return writeVersion_;
            }

            //- Default write compression
            IOstream::compressionType writeCompression() const
            {
                return writeCompression_;
            }

            //- Default graph format
            const word& graphFormat() const
            {
                return graphFormat_;
            }

            //- Read control dictionary, update controls and time
            virtual bool read();

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Return the location of "dir" containing the file "name".
            //  (Used in reading mesh data)
            word findInstance
            (
                const fileName& dir,
                const word& name,
                const IOobject::readOption rOpt = IOobject::MUST_READ
            ) const;

            //- Search tha case for valid time directories
            instantList times() const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const instant&) const;

            //- Search the case for the time closest to the given time
            instant findClosestTime(const scalar) const;

            //- Write all the objects
            virtual bool write() const;


        // Access

            //- Return time name of given scalar time
            static word timeName(const scalar);

            //- Return current time name
            virtual word timeName() const;

            //- Search a given directory for valid time directories
            static instantList findTimes(const fileName&);

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;


        // Check

            //- Return true if run should continue
            virtual bool run() const;

            //- Return true if end of run
            virtual bool end() const;


        // Edit

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time&);

            //- Reset the time and time-index
            virtual void setTime(const instant&, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime(const dimensionedScalar&,const label newIndex);

            //- Reset the time and time-index
            virtual void setTime(const scalar, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar&);

            //- Reset end time
            virtual void setEndTime(const scalar);

            //- Reset time step
            virtual void setDeltaT(const dimensionedScalar&);

            //- Reset time step
            virtual void setDeltaT(const scalar);

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Reset time after sub-cycling back to given TimeState
            virtual void endSubCycle(const TimeState&);


    // Member operators

        //- Set deltaT to that specified and increment time
        virtual Time& operator+=(const dimensionedScalar&);

        //- Set deltaT to that specified and increment time
        virtual Time& operator+=(const scalar);

        //- Prefix increment
        virtual Time& operator++();

        //- Postfix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

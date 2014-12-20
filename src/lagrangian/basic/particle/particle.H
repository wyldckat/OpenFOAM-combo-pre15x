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
    particle

Description

\*---------------------------------------------------------------------------*/

#ifndef particle_H
#define particle_H

#include "vector.H"
#include "IDLList.H"
#include "labelList.H"
#include "pointField.H"
#include "faceList.H"
#include "typeInfo.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class particle>
class Cloud;


// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class particleType>
class particle;

template<class particleType>
Ostream& operator<<
(
    Ostream&,
    const particle<particleType>&
);


/*---------------------------------------------------------------------------*\
                           Class particle Declaration
\*---------------------------------------------------------------------------*/

template<class particleType>
class particle
:
    public IDLList<particleType>::link
{
    // Private data

        //- Reference to the particle cloud
        const Cloud<particleType>& cloud_;

        //- Position of particle
        vector position_;

        //- Index of the cell it is
        label celli_;

        //- is particle on boundary/outside domain
        bool onBoundary_;

    // Private member functions

        //- Return the 'lambda' value for the position, p, on the face,
        // where, p = from + lamda*(to - from)
        // for non-static meshes
        inline scalar lambda
        (
            const vector& from,
            const vector& to,
            const label facei,
            const scalar fraction
        );

        //- Return the 'lambda' value for the position, p, on the face,
        // where, p = from + lamda*(to - from) 
        // for static meshes
        inline scalar lambda
        (
            const vector& from,
            const vector& to,
            const label facei
        );

        //- Return the faces between position and cell centre
        labelList findFaces
        (
            const vector& position
        );

        //- Return the faces between position and cell centre
        labelList findFaces
        (
            const vector& position,
            const label celli,
            const scalar fraction
        );

        //- Convert global addressing to the processor patch
        //  local equivalents
        void prepareForParallelTransfer
        (
            const label patchi,
            const label facei
        );

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        void correctAfterParallelTransfer
        (
            const label patchi,
            const label facei
        );


public:

    friend class Cloud<particleType>;


    //- Runtime type information
    TypeName("particle");


    // Constructors

        //- Construct from components
        particle
        (
            const Cloud<particleType>&,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        particle
        (
            const Cloud<particleType>&,
            Istream&
        );


    // Member Functions

        // Access

            //- Return true if particle is in cell
            inline bool inCell();

            //- Return true if position is in cell i
            inline bool inCell
            (
                const vector& position,
                const label celli,
                const scalar fraction
            );

            //- Return current particle position
            inline const vector& position() const;

            //- Return current cell particle is in
            inline label cell() const;

            //- Return reference to the particle cloud
            inline const Cloud<particleType>& cloud() const;


        // Check

            //- Is the particle on the boundary/( or outside the domain)?
            inline bool onBoundary() const;

            //- Which patch is particle on
            inline label patch(const label facei) const;

            //- Which face of this patch is this particle on
            inline label patchFace(const label patchi, const label facei) const;
    
            //- The nearest distance to a wall that 
            // the particle can be in the n direction
            virtual scalar wallImpactDistance(const vector& n) const;


        // Track

            //- Track particle to end of trajectory
            //  or until it hits the boundary,
            // on entry 'fraction' is the fraction of the integration step
            // at which the tracking starts (only used for moving meshes)
            // and on exit it contains the fraction of the traveled trajectory
            label track
            (
                const vector& endPosition,
                scalar& fraction
            );

            //- Track particle to a given position
            // If particle reaches that position without hitting a face
            // we return -1, otherwise track to face and return face index
            // on entry 'fraction' is the fraction of the integration step
            // at which the tracking starts (only used for moving meshes)
            // upon exit fraction contains the fraction of the actual
            // distance moved
            label trackToFace
            (
                const vector& endPosition,
                scalar& fraction
            );

        // Particle operations

            //- Transform the position the particle
            //  according to the given transformation tensor
            void transformPosition(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);


    // Ostream Operator

        friend Ostream& operator<< <particleType>
        (
            Ostream&,
            const particle<particleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleI.H"

#define defineParticleTypeNameAndDebug(Type, DebugSwitch)                     \
    template<>                                                                \
    const Foam::word particle<Type>::typeName(#Type);                         \
    template<>                                                                \
    int particle<Type>::debug(Foam::debug::debugSwitch(#Type, DebugSwitch));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "particle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    lduMatrix

Description
    lduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing arrays must
    be supplied for the upper and lower triangles.  Diagonal, ICCG and BCG
    solvers are supplied and selected automatically by the solve member
    function according the the particular matrix structure.

    It would be better if this class were organised as a hierachy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.  A job for the future.

SourceFiles
    lduMatrixATmul.C
    lduMatrix.C
    lduMatrixHoperations.C
    lduMatrixOperations.C
    lduMatrixSolver.C
    lduMatrixTests.C
    lduMatrixUpdateCoupledInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef lduMatrix_H
#define lduMatrix_H

#include "primitiveFields.H"
#include "FieldField.H"
#include "lduAddressing.H"
#include "processorTopology.H"
#include "lduCoupledInterfacePtrsList.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lduMatrix Declaration
\*---------------------------------------------------------------------------*/

class lduMatrix
{

public:

    //- Class returned by the solver
    //  containing performance statistics
    class solverPerformance
    {
        word   solverName_;
        word   fieldName_;
        scalar initialResidual_;
        scalar finalResidual_;
        label  noIterations_;
        bool   converged_;
        bool   singular_;


    public:

        // Constructors

            solverPerformance()
            :
                initialResidual_(0),
                finalResidual_(0),
                noIterations_(0),
                converged_(false),
                singular_(false)
            {}


            solverPerformance
            (
                const word&  solverName,
                const word&  fieldName,
                const scalar iRes = 0,
                const scalar fRes = 0,
                const label  nIter = 0,
                const bool   converged = false,
                const bool   singular = false
            )
            :
                solverName_(solverName),
                fieldName_(fieldName),
                initialResidual_(iRes),
                finalResidual_(fRes),
                noIterations_(nIter),
                converged_(converged),
                singular_(singular)
            {}


        // Member functions

            //- Return solver name
            const word& solverName() const
            {
                return solverName_;
            }

            //- Return initial residual
            scalar initialResidual() const
            {
                return initialResidual_;
            }

            //- Return initial residual
            scalar& initialResidual()
            {
                return initialResidual_;
            }


            //- Return final residual
            scalar finalResidual() const
            {
                return finalResidual_;
            }

            //- Return final residual
            scalar& finalResidual()
            {
                return finalResidual_;
            }


            //- Return number of iterations
            label nIterations() const
            {
                return noIterations_;
            }

            //- Return number of iterations
            label& nIterations()
            {
                return noIterations_;
            }


            //- Has the solver converged?
            bool converged() const
            {
                return converged_;
            }

            //- Is the matrix singular?
            bool singular() const
            {
                return singular_;
            }

            //- Convergence test
            bool checkConvergence
            (
                const scalar tolerance,
                const scalar relTolerance
            );

            //- Singularity test
            bool checkSingularity(const scalar residual);

            //- Print summary of solver performance
            void print() const;
    };


    class solver
    {

    protected:

        // Protected data

            word fieldName_;
            scalarField& psi_;
            const lduMatrix& matrix_;
            const scalarField& source_;
            const FieldField<Field, scalar>& coupleBouCoeffs_;
            const FieldField<Field, scalar>& coupleIntCoeffs_;
            const lduCoupledInterfacePtrsList& interfaces_;
            const direction cmpt_;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    scalarField& psi,
                    lduMatrix& matrix,
                    const scalarField& source,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduCoupledInterfacePtrsList& interfaces,
                    const direction cmpt,
                    Istream& solverData
                ),
                (
                    fieldName,
                    psi,
                    matrix,
                    source,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    cmpt,
                    solverData
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                asymMatrix,
                (
                    const word& fieldName,
                    scalarField& psi,
                    lduMatrix& matrix,
                    const scalarField& source,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduCoupledInterfacePtrsList& interfaces,
                    const direction cmpt,
                    Istream& solverData
                ),
                (
                    fieldName,
                    psi,
                    matrix,
                    source,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    cmpt,
                    solverData
                )
            );


        // Constructors

            solver
            (
                const word& fieldName,
                scalarField& psi,
                const lduMatrix& matrix,
                const scalarField& source,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const direction cmpt
            )
            :
                fieldName_(fieldName),
                psi_(psi),
                matrix_(matrix),
                source_(source),
                coupleBouCoeffs_(coupleBouCoeffs),
                coupleIntCoeffs_(coupleIntCoeffs),
                interfaces_(interfaces),
                cmpt_(cmpt)
            {}


        // Selectors

            //- Return a new solver
            static autoPtr<solver> New
            (
                const word& fieldName,
                scalarField& psi,
                lduMatrix& matrix,
                const scalarField& source,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const direction cmpt,
                Istream& solverData
            );


        // Destructor

            virtual ~solver()
            {}


        // Member functions

            virtual solverPerformance solve() = 0;
    };


private:

    // private data

        // L-U addressing
        const lduAddressing& lduAddr_;

        // lduCoupledInterface evaluation schedule
        const patchScheduleList& lduCoupledInterfaceSchedule_;

        //- Face matrix elements
        scalarField *lowerPtr_, *diagPtr_, *upperPtr_;


public:

    // Static data

        // Declare name of the class and it's debug switch
        ClassName("lduMatrix");

        //- Large scalar for the use in solvers
        static const scalar great_;

        //- Small scalar for the use in solvers
        static const scalar small_;


    // Constructors

        //- Construct given a field to solve, L-U matrix addressing
        //  and tolerence
        lduMatrix
        (
            const lduAddressing& ldu,
            const patchScheduleList& lduCoupledInterfaceSchedule
        );

        //- Construct as copy
        lduMatrix(const lduMatrix&);

        //- Construct as copy or re-use as specified.
        lduMatrix(lduMatrix&, bool reUse);

        //- Construct from Istream
        lduMatrix
        (
            const lduAddressing& ldu,
            const patchScheduleList& lduCoupledInterfaceSchedule,
            Istream&
        );


    // Destructor

        ~lduMatrix();


    // Member functions

        // Access

            scalarField& lower();
            scalarField& diag();
            scalarField& upper();

            const scalarField& lower() const;
            const scalarField& diag() const;
            const scalarField& upper() const;

            const lduAddressing& lduAddr() const
            {
                return lduAddr_;
            }

            const patchScheduleList& lduCoupledInterfaceSchedule() const
            {
                return lduCoupledInterfaceSchedule_;
            }

            bool hasDiag() const
            {
                return (diagPtr_);
            }

            bool hasUpper() const
            {
                return (upperPtr_);
            }

            bool hasLower() const
            {
                return (lowerPtr_);
            }

            bool diagonal() const
            {
                return (diagPtr_ && !lowerPtr_ && !upperPtr_);
            }

            bool symmetric() const
            {
                return (diagPtr_ && (!lowerPtr_ && upperPtr_));
            }

            bool asymmetric() const
            {
                return (diagPtr_ && lowerPtr_ && upperPtr_);
            }


        // operations

            void sumDiag();
            void negSumDiag();

            //- Relax matrix (for steady-state solution)
            //  alpha=1 : diagonally equal
            //  alpha<1 :    ,,      dominant
            void relax
            (
                const FieldField<Field, scalar>& intCoeffsCmptAvg,
                const FieldField<Field, scalar>& magCoupleBouCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const scalar alpha=1.0
            );


            //- Matrix multiplication with updated coupled interfaces.
            void Amul
            (
                scalarField&,
                const tmp<scalarField>&,
                const FieldField<Field, scalar>&,
                const lduCoupledInterfacePtrsList&,
                const direction cmpt
            ) const;
    
            //- Matrix transpose multiplication with updated coupled interfaces.
            void Tmul
            (
                scalarField&,
                const tmp<scalarField>&,
                const FieldField<Field, scalar>&,
                const lduCoupledInterfacePtrsList&,
                const direction cmpt
            ) const;


            //- Initialise the update of coupled interfaces
            //  for matrix operations
            void initMatrixInterfaces
            (
                const FieldField<Field, scalar>& coupleCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const scalarField& psiif,
                scalarField& result,
                const direction cmpt
            ) const;
    
            //- Update coupled interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const FieldField<Field, scalar>& coupleCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const scalarField& psiif,
                scalarField& result,
                const direction cmpt
            ) const;
    

            tmp<scalarField> residual
            (
                const scalarField& psi,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const direction cmpt
            ) const;

            tmp<scalarField> residual
            (
                const scalarField& psi,
                const scalarField& source,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const lduCoupledInterfacePtrsList& interfaces,
                const direction cmpt
            ) const;


            template<class Type>
            tmp<Field<Type> > H(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > H(const tmp<Field<Type> >&) const;

            template<class Type>
            tmp<Field<Type> > faceH(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > faceH(const tmp<Field<Type> >&) const;


    // Member operators

        void operator=(const lduMatrix&);

        void negate();

        void operator+=(const lduMatrix&);
        void operator-=(const lduMatrix&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const lduMatrix&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lduMatrixHoperations.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

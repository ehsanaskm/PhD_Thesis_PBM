/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    coalescenceLaw

Description
    Bubbles coalescence selector

SourceFiles
    coalescenceLaw.C
    newcoalescenceLaw.C

Author
    brahim selma, PhD 

\*---------------------------------------------------------------------------*/

#ifndef coalescenceLaw_H
#define coalescenceLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coalescenceLaw Declaration
\*---------------------------------------------------------------------------*/

class coalescenceLaw
{
    // Private data

        //- Name
        const word name_;

		  //- Reference to volume fracion of dispersed phase 
		  const volScalarField& alpha_;


    // Private Member Functions

        //- Disallow copy construct
        coalescenceLaw(const coalescenceLaw&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return alpha_.mesh();
        }

		  //- Return reference to alpha
		  const volScalarField& alpha() const
        {
            return alpha_;
        }
        
       
        
public:

    //- Runtime type information
    TypeName("coalescenceLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceLaw,
            dictionary,
            (
                const word name,
                const volScalarField& alpha,
                const dictionary& dict
            ),
            (name, alpha, dict)
        );


    // Selectors

        //- Return a reference to the selected breakup model
        static autoPtr<coalescenceLaw> New
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& dict
        );

    // Constructors

        //- Construct from dictionary
        coalescenceLaw
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& dict
        );

    // Destructor

        virtual ~coalescenceLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Calculate and return coalescence rate 
        virtual dimensionedScalar coalescenceRate
        (
				scalar&, scalar&, volScalarField&
			) const = 0;
       
       
       // references to liquid and gas density
        virtual dimensionedScalar rhoa() const = 0;
        virtual dimensionedScalar rhob() const = 0;
       
      
       
        //- If the breakup coalescence models needs population model
        virtual bool populationModelNeeded() const
        {
            return false;
        };

        //- Correct the bubble breakup model
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

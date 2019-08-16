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
    populationBalanceReturn

Description
    Method to calcluate the Sauter mean diameter, returning the
    d32 values

SourceFiles
    populationBalanceReturn.C
    newpopulationBalanceReturn.C

Author
    Brahim SELMA 
    USh 2015

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceReturn_H
#define populationBalanceReturn_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class populationBalanceModel;

/*---------------------------------------------------------------------------*\
                         Class populationBalanceReturn Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceReturn
{
    // Private data

        //- Name
        const word name_;


    // Private Member Functions

        //- Disallow copy construct
        populationBalanceReturn(const populationBalanceReturn&);

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceReturn&);


  //protected:

public:

    //- Runtime type information
    TypeName("populationBalanceReturn");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            populationBalanceReturn,
            dictionary,
            (
					const word name,
					populationBalanceModel& populationBalanceModel
            ),
            (name, populationBalanceModel)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<populationBalanceReturn> New
        (
				const word& name,
				populationBalanceModel& populationBalanceModel
        );


    // Constructors

        //- Construct from dictionary
        populationBalanceReturn
        (
				const word& name,
				populationBalanceModel& populationBalanceModel
        );


    // Destructor

        virtual ~populationBalanceReturn()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

		  //- If population is active
        virtual bool populationActive() const
        {
			    return true;
		  }

        //- Update the yield stress
        //- called at end of time-step
        virtual void updateSauterDiameter() = 0;

       
        //- Correct model
        virtual void correct() = 0;

            

        //- Write model dictionary
        //virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

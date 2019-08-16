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
    populationBalanceModel

Description
    Class of population balance methods 

SourceFiles
    populationBalanceModel.C

Author
    brahim selma
    Université de Sherbrooke

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "IOdictionary.H"
#include "populationBalanceReturn.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "breakupLaw.H"
#include "coalescenceLaw.H"
#include "Switch.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
:
  public IOdictionary
{
    // Private data

		  //- Reference to volume fraction of dispersed phases 
		  const volScalarField& alpha_;
		  
        //- breakup law ptr
        autoPtr<breakupLaw> breakupLawPtr_;

		  //- coalescence law ptr
		  autoPtr<coalescenceLaw> coalescenceLawPtr_;

        //- Run-time selectable method to calculate f and source term of PBE
        autoPtr<populationBalanceReturn> populationBalanceReturnPtr_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        populationBalanceModel(const populationBalanceModel&);

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceModel&);


public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Constructors

        //- Construct from components
        populationBalanceModel
        (
				PtrList<volScalarField>& S,
				PtrList<volScalarField>& f,
				volScalarField& alpha
        );


    // Destructor

        virtual ~populationBalanceModel();


    // Member Functions

        // Access

		  //- If population balance equation is active
        bool populationActive() const
        {
				if(populationBalanceReturnPtr_.valid())
				{
					return populationBalanceReturnPtr_->populationActive();
				}
				return false;
        }

        //- Return reference to volume fraction 
        const volScalarField& alpha() const
        {
            return alpha_;
        }
       
        //- Return breakup model
        const breakupLaw& breakupModel() const
        {
            return breakupLawPtr_();
        }
        
        //- Return coalescence model
        const coalescenceLaw& coalescenceModel() const
        {
            return coalescenceLawPtr_();
        }
        
        //- Return breakup rate   
        tmp<volScalarField> breakupRate(scalar&, scalar&, scalar&) const; 
        
        //- Return breakup frequency
        tmp<volScalarField> breakupFrequency(scalar&, scalar&, scalar&) const; 
        
        //- Return coalescence rate 
        tmp<volScalarField> coalescenceRate(scalar&, scalar&, scalar&) const; 
        
        //--- reference to gas density using pointer 
        dimensionedScalar rhoa() const
        {
            return breakupLawPtr_->rhoa();
        }
        
        //--- reference to liquid density using pointer 
        dimensionedScalar rhob() const
        {
            return breakupLawPtr_->rhob();
        }
        
        //- Correct d32
        virtual void correct();

        //- Update yield stress
        void updateSauterDiameter();

        //- Read plasticityProperties dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    classMethod

SourceFiles
    classMethod.C
    
Description
	 solve population balance equation  	
    
References

    [1] 
    [2] 
     
Author
    Ehsan Askari
	 Sep. 2016

\*---------------------------------------------------------------------------*/

#ifndef classMethod_H
#define classMethod_H

#include "populationBalanceReturn.H"
#include "breakupLaw.H"
//Added by Ehsan
#include "diameterModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class populationBalanceModel;

/*---------------------------------------------------------------------------*\                                                                                         
                         Class classMethod Declaration
\*---------------------------------------------------------------------------*/

class classMethod
:
    public populationBalanceReturn
{
    // Private data

			
			// Reference to populationBalance model
			populationBalanceModel& populationBalanceModel_;

			// volume fraction 
			volScalarField alpha_;

			//- Menétrey-willam dictionary to define all constant
			dictionary dict_;
			
			// Sauter diameter d32
			volScalarField d32_;


            
                       volScalarField dsauter;
                        
		  
			//- mean diameter
			dimensionedScalar dm_;


                       // correction parameters
                       dimensionedScalar ResidualAlphaForDsauter_;
                       dimensionedScalar ResidualAlphaForAdjust_;
                       dimensionedScalar ResidualAlphaForCorrect_;

			//- number of classes  (is a label not scalar!!)
			label Nc_;
    
			//- coalescence: initial thickness
			dimensionedScalar h0_;
   
			//- coalescence: rupture thickness  	
			dimensionedScalar hf_;
    
			//-
			dimensionedScalar beta1_;
   
			//-
			dimensionedScalar beta2_;
   
			//-
			dimensionedScalar betaPB_;
   
			//-
			dimensionedScalar k1_;
     
			//- surface tension 
			dimensionedScalar sigma_;
      
			//- breakup coefficient
			dimensionedScalar breakCoeff_;
  
			//- coalescence coefficient
			dimensionedScalar coalCoeff_;
     
			//- Ramkrishna model paramter: equivalent to V(i+1) = S*V(i) for geometric evolution 
			dimensionedScalar S_;
			
			//-
			dimensionedScalar dMin_;
			dimensionedScalar dMax_;
			
			//- Max of iterations
			label maxIters_; 
			
			//- tolerance
			scalar loopTolerance_;

            //           PtrList<volScalarField> f;


                  
                      			
			//- class method selector
			//word populationBalanceMethodName_; 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        classMethod(const classMethod&);

        //- Disallow default bitwise assignment
        void operator=(const classMethod&);

public:
   

     
    //- Runtime type information
    TypeName("classMethod");

    // Static data members


    // Constructors

        //- Construct from dictionary
        classMethod
        (
            const word& name,
            populationBalanceModel& populationBalanceModel
        );


    // Destructor

        virtual ~classMethod();


     


    // Member Functions
/////////////////////////////////////////////////////////////////////
    //    virtual tmp<volScalarField> dSauter() const;
//////////////////////////////////////////////////////////////////////        

        //- Correct the rheological model
        virtual void correct();

		  scalar funcTest(scalar&, scalar&, scalar&) const;
	
		  scalar kronecker(label&, label&) const;   
	
		  //------------------------------ bubbles breakup kernel -----------------------//
		  virtual void breakupKernel
		  ( 
				PtrList<volScalarField>&, 
				PtrList<volScalarField>&,
				const volScalarField&
		  ) const;  
		  
		  //------------------------------ bubbles coalescence kernel -------------------//
		  virtual void coalescenceKernel
		  ( 
				PtrList<volScalarField>&, 
				PtrList<volScalarField>&,
				const volScalarField&,
				const volScalarField&
		  ) const;
		  	
	     //------------------------------ Sauter diameter function ---------------------//
	     tmp<volScalarField> SauterDiameter(PtrList<volScalarField>&, const volScalarField&) const ;
           
                 
	    
	     //--- 
	     void adjust
			(
				PtrList<volScalarField>&, 
				PtrList<volScalarField>&, 
				const volScalarField&
			)
			const; 

	    
        //- Return slave friction value fraction
        virtual void updateSauterDiameter();

        
        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

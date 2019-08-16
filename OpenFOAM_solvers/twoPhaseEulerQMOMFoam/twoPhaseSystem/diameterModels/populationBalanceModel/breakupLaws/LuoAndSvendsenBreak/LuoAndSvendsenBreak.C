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

\*---------------------------------------------------------------------------*/

#include "LuoAndSvendsenBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "incompleteGammaFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LuoAndSvendsenBreak, 0);
    addToRunTimeSelectionTable(breakupLaw, LuoAndSvendsenBreak, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::LuoAndSvendsenBreak::LuoAndSvendsenBreak
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& dict
)
:
    breakupLaw(name, alpha, dict),
    dictLuoAndSvendsenBreakup_
    (
        dict.subDict("LuoAndSvendsenCoeffs")
    ),
    
    // Luo and Svendsen model parameters
    rhoa_(dictLuoAndSvendsenBreakup_.lookup("rhoa")),
    rhob_(dictLuoAndSvendsenBreakup_.lookup("rhob")),
    sigma_(dictLuoAndSvendsenBreakup_.lookup("sigma")),
    beta_(readScalar(dictLuoAndSvendsenBreakup_.lookup("beta"))),
    epsMax_(dictLuoAndSvendsenBreakup_.lookup("epsilonMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LuoAndSvendsenBreak::~LuoAndSvendsenBreak()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::LuoAndSvendsenBreak::rhoa() const
{
    return rhoa_; 
}

Foam::dimensionedScalar Foam::LuoAndSvendsenBreak::rhob() const
{
    return rhob_;
}

Foam::dimensionedScalar Foam::LuoAndSvendsenBreak::breakupRate(volScalarField& vf, scalar& d1, scalar& d2) const

{


  const volScalarField& epsilon = 
	mesh().objectRegistry::lookupObject<volScalarField>("epsilon.water");   			
	
  dimensionedScalar result("SMALL",dimMass/pow3(dimLength)/dimTime, SMALL);

 // scalar fbv =
//					pow3(min(d1,d2)/(max(d1,d2)+SMALL));

   scalar fbv=0.5;
  scalar Cf = 
               pow(fbv,2./3.) + pow((1.-fbv), 2./3.) - 1.;
  
  scalar a = 
			k_.value()*beta_*max(vf).value()*Foam::pow(min(average(epsilon).value(),epsMax_.value())/sqr(d1), 1./3.);
                        //0.2727*k_.value()*beta_*Foam::pow(average(epsilon).value()/(sqr(d1)+SMALL), 1./3.); 
                        
			
			
  scalar b = 
			12.0*Cf*sigma_.value()/
			(  
			         beta_*rhob_.value()*Foam::pow(min(average(epsilon).value(),epsMax_.value()), 2./3.)
                                 //    beta_*rhob_.value()*Foam::pow(average(epsilon).value(), 2./3.)
                                 
			                          *Foam::pow(d1, 5./3.) 
			     + SMALL
			);			
			
  scalar breakRate = a*
			 (
					  (3./11.)*(1./Foam::pow(b, 8./11.))*incompleteGammaFunction().gammQ(8./11., b)
					+ (6./11.)*(1./Foam::pow(b, 5./11.))*incompleteGammaFunction().gammQ(5./11., b)
					+ (3./11.)*(1./Foam::pow(b, 2./11.))*incompleteGammaFunction().gammQ(2./11., b)
			 );
  
  if(breakRate < 0.0)
  {
	  result = dimensionedScalar
				  (
						"zero",
						dimMass/pow3(dimLength)/dimTime,
						0.0
				  );
  } 
  else
  {
	  result = dimensionedScalar
				  (
						"value",
						dimMass/pow3(dimLength)/dimTime,
						breakRate
				  );
  }
    
  return result ;
}

Foam::dimensionedScalar Foam::LuoAndSvendsenBreak::breakupFrequency(volScalarField& vf, scalar& d, scalar& fbv) const
{
 // const volScalarField& epsilon = 
			//mesh().objectRegistry::lookupObject<volScalarField>("epsilon");  

  const volScalarField& epsilon = 
			mesh().objectRegistry::lookupObject<volScalarField>("epsilon.water");			
  dimensionedScalar result("SMALL",pow(dimTime,-1), SMALL);
  
  scalar Cf = pow(fbv, 2./3.) + pow((1.-fbv),2./3.) - 1.;
  
  scalar b  = 
		 12.0*Cf*sigma_.value()/(
			                          k_.value()*rhob_.value()
			                         *Foam::pow(min(average(epsilon).value(),epsMax_.value()), 2./3.)
			                         *Foam::pow(d,5./3.) 
			                         + SMALL
			                     );


  scalar a = 
            k_.value()*beta_*max(vf).value()*pow(min(average(epsilon).value(),epsMax_.value())/sqr(d), 1./3.);
  
  
  scalar breakFrequency = a*
			(
					  (3./11.)*(1./Foam::pow(b, 8./11.))*incompleteGammaFunction().gammQ(8./11., b)
					+ (6./11.)*(1./Foam::pow(b, 5./11.))*incompleteGammaFunction().gammQ(5./11., b)
					+ (3./11.)*(1./Foam::pow(b, 2./11.))*incompleteGammaFunction().gammQ(2./11., b)

			);
  
  if(breakFrequency < 0.0)
  {
	  result = dimensionedScalar
				  (
						"zero",
						pow(dimTime,-1),
						0.0
				  );
  } 
  else
  {
	  result = dimensionedScalar
				  (
						"value",
						pow(dimTime,-1),
						breakFrequency
				  );
  }
    
  return result ;
 
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LuoSvendsen.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        LuoSvendsen,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen
::LuoSvendsen
(
    const dictionary& dict
)
:
    breakupKernel(dict),
   Cb_(dict.lookup("Cb")),
    epsilonExp_(readScalar(dict.lookup("epsilonExp"))),
   nuExp_(readScalar(dict.lookup("nuExp"))),
    k_(dict.lookup("k")),
    beta_(readScalar(dict.lookup("beta"))),
    rhoa_(dict.lookup("rhoa")),
    rhob_(dict.lookup("rhob")),
    sigma_(dict.lookup("sigma")),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::~LuoSvendsen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::Kb
(
    const volScalarField& abscissa
) const
{   
/*    const compressible::turbulenceModel& fluidTurb =
        abscissa.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
*/
  const volScalarField& epsilon = 
	abscissa.mesh().lookupObject<volScalarField>("epsilon.water"); 

  const volScalarField& alpha = 
	abscissa.mesh().lookupObject<volScalarField>("alpha.air");

  const fvMesh& mesh = abscissa.mesh();


 // volScalarField result("SMALL",pow(dimTime,-1), SMALL);

  volScalarField result
  (
			IOobject
			(
				"result",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			  alpha*
			dimensionedScalar
			(
			   "one", 
			   pow(dimTime,-1), 
			   1.
		        )
                         

  );



  const  volScalarField& vf = 1-alpha;



  scalar fbv=0.5;
  
  scalar Cf = pow(fbv, 2./3.) + pow((1.-fbv),2./3.) - 1.; 


  volScalarField a = 
			k_*beta_*vf*Foam::pow(epsilon/sqr(abscissa), 1./3.);


  volScalarField b =   12.0*Cf*sigma_/
			(  
			         beta_*rhob_*Foam::pow(epsilon, 2./3.)
                                 
			                          *Foam::pow(abscissa, 5./3.) 
			     + dimensionedScalar("SMALL", dimMass/pow(dimTime,2), SMALL)
			);

  volScalarField breakRate=a;

  forAll(breakRate,II)
  {
     
    breakRate[II]= breakRate[II]*
			 (
					  (3./11.)*(1./Foam::pow(b[II], 8./11.))*incompleteGammaFunction().gammQ(8./11., b[II])
					+ (6./11.)*(1./Foam::pow(b[II], 5./11.))*incompleteGammaFunction().gammQ(5./11., b[II])
					+ (3./11.)*(1./Foam::pow(b[II], 2./11.))*incompleteGammaFunction().gammQ(2./11., b[II])
					- (1.0)*incompleteGammaFunction().gammQ(8./11., b[II])
					- (2.0)*(Foam::pow(b[II], 3./11.))*incompleteGammaFunction().gammQ(5./11., b[II])
                                        - (1.0)*(Foam::pow(b[II], 6./11.))*incompleteGammaFunction().gammQ(2./11., b[II])
			 ); 
   
  }

                      			
	forAll(breakRate,ii)
	{
		  if(breakRate[ii] < 0.0)
		  {
			  result[ii] = 1e-15;
		  } 

		  else
		  {
			  result[ii] = breakRate[ii]+SMALL;
						 
		  }
         }
			
  
  Info<<max(result)<<endl;
  Info<<min(result)<<endl;
  Info<<average(result)<<endl;

    dimensionedScalar minAbs
    (
        "minAbs", 
        abscissa.dimensions(), 
        minAbscissa_.value()
    );

    return Cb_*pos(abscissa - minAbs);


   // return Cb_*pow(epsilon, epsilonExp_)
     //   *pow(epsilon, nuExp_)*pow(abscissa, sizeExp_);

  //  return result;
}

// ************************************************************************* //

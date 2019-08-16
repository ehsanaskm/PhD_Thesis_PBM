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

#include "luoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(luoModel, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        luoModel,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::luoModel
::luoModel
(
    const dictionary& dict
)
:
    aggregationKernel(dict),
    rhoa_(dict.lookup("rhoa")),
    rhob_(dict.lookup("rhob")),
    sigma_(dict.lookup("sigma")),
    epsMax_(dict_.lookup("epsilonMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::luoModel
::~luoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::luoModel::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{  
	 const fvMesh& mesh = abscissa1.mesh();

	// const volScalarField& epsilon = 
		//abscissa1.mesh().lookupObject<volScalarField>("epsilon.water"); 

	  const volScalarField& alpha = 
		abscissa1.mesh().lookupObject<volScalarField>("alpha");



	scalar omega_coal, xi, p_coal, weber, num, den, u_i, u_j, u_i_j,coalRate,d_i,d_j;

    //dimensionedScalar result("coalRate",dimensionSet(0,3,-1,0,0,0,0),0.0); 

      dimensionedScalar result("coalRate",dimMass/pow3(dimLength)/dimTime,SMALL); 


  volScalarField coalRate2
  (
			IOobject
			(
				"coalRate2",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			  alpha*
			dimensionedScalar
			(
			   "one", 
			   dimensionSet(0,3,-1,0,0,0,0), 
			   1.
		        )
                         

  );



	xi=max(abscissa1).value()/1000000/(max(abscissa2).value()+1e-15)/1000000;

        

        d_i=average(abscissa1).value()/1000000;
        d_j=average(abscissa2).value()/1000000;

        

	//u_i = 1.43*pow(average(epsilon).value()*d_i, 1./3);

	//u_j = 1.43*pow(average(epsilon).value()*d_j/1000, 1./3);
	u_i = d_i;

	u_j =  d_j;

	u_i_j = pow((u_i+u_j),3);

	weber = rhob_.value()*d_i*pow(u_i_j,2.0)/sigma_.value();

	num = pow((0.75*(1.+pow(xi,2.))*(1.+pow(xi,3.))), 0.5);
        den = pow((0.5+(rhoa_/rhob_).value()),0.5)*pow((1.+xi),3.);
        p_coal = exp(-(num/den)*pow(weber, 0.5));
 
        omega_coal = (4/3)*u_i_j*534.7;

        coalRate = omega_coal*p_coal;



       forAll(alpha,cellI)
       {

	   if(coalRate < 0.0) 
	   {
		  coalRate2[cellI]=0.0;
	   } 

	   else
	   {
		  coalRate2[cellI]=coalRate;
	   }

       }

  
 

     

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "constantAggregationK",
                    abscissa1.mesh().time().timeName(),
                    abscissa1.mesh()
                ),
                   coalRate2
                
            )
        );

}

// ************************************************************************* //

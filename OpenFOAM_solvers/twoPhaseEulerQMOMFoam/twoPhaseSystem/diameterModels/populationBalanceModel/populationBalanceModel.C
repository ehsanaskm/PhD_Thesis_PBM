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

#include "populationBalanceModel.H"
#include "volFields.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(populationBalanceModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

populationBalanceModel::populationBalanceModel
(
	PtrList<volScalarField>& S,
	PtrList<volScalarField>& f,
	volScalarField& alpha
)
:
    IOdictionary
    (
        IOobject
        (
            "populationBalanceProperties",
            alpha.time().constant(),
            alpha.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    alpha_(alpha),
    breakupLawPtr_
    (
		breakupLaw::New
		(
			"law", 
			alpha_, 
			subDict("breakupModels")
		)
	 ),
	 coalescenceLawPtr_
    (
		coalescenceLaw::New
		(
			"law", 
			alpha_, 
			subDict("coalescenceModels")
		)
	 ),
    populationBalanceReturnPtr_(NULL)


{
	 Info << "Creating population balance methods" << endl;
  
	 if (breakupLawPtr_->populationModelNeeded())
	 {
	 	populationBalanceReturnPtr_ =
		populationBalanceReturn::New(lookup("populationModel"), *this);
	 }
	 else if (found("populationModel"))
	 {
		Warning << "populationModel is specified in constant/populationBalanceProperties"
		<< " but population is not active in the current model"
		<< endl;
	 }
}


populationBalanceModel::~populationBalanceModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void populationBalanceModel::correct()
{
  if (populationBalanceReturnPtr_.valid())
    {
			populationBalanceReturnPtr_->correct();
    }
  else
    {
			breakupLawPtr_->correct();
    }
}

void populationBalanceModel::updateSauterDiameter()
{
return (populationBalanceReturnPtr_->updateSauterDiameter());

  if (populationBalanceReturnPtr_.valid())
    {
		  populationBalanceReturnPtr_->updateSauterDiameter();
    }

}

bool populationBalanceModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


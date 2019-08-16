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
	 Selector of population balance methods. Available methods are:
	 1. class method
	 2. sqmom: standard method of moments
	 3. qmom: quadrature method of moments 
	 4. dqmom: direct quadrature method of moments 
    
Author
	 Brahim Selma, PhD
	 Université de Sherbrooke

\*---------------------------------------------------------------------------*/

#include "populationBalanceReturn.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<populationBalanceReturn> populationBalanceReturn::New
(
	const word& name,
	populationBalanceModel& populationBalanceModel
)
{
    Info<< "\tPopulation balance return method: " << name << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
			FatalIOErrorIn
			(
				" populationBalanceReturn::New(\n"
				" const word& name\n"
				" populationBalanceModel& populationBalanceModel\n"
				")",
				"populationBalanceModelCoeffs"
		)
		<< "Unknown populationBalanceReturn type "
		<< name << endl << endl
		<< "Valid  populationBalanceReturns methods are : " << endl
		<< dictionaryConstructorTablePtr_->sortedToc()
		<< exit(FatalIOError);
    }

    return autoPtr<populationBalanceReturn>
    (
			cstrIter()
			(
				name,
				populationBalanceModel
			)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

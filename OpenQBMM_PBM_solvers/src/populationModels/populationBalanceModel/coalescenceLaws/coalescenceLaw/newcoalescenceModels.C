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
    coalescence Law

Description
	 Selecting coalescence models

\*---------------------------------------------------------------------------*/

#include "coalescenceLaw.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<coalescenceLaw> coalescenceLaw::New
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& dict
)
{
    word coalescenceModelName = dict.lookup("modelName");

    Info<< "Selecting coalescence model name " << coalescenceModelName << endl;

    dictionaryConstructorTable::iterator 
    cstrIter = dictionaryConstructorTablePtr_->find(coalescenceModelName);
	
    if (
				cstrIter == dictionaryConstructorTablePtr_->end()
		 )
    {
        FatalIOErrorIn
        (
            "coalescenceLaw::New(\n"
            "    const word& name,\n"
            "    const volScalarField& alpha,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown coalescence model name "
            << coalescenceModelName << endl << endl
            << "Valid coalescence models are: " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coalescenceLaw>
    (
			cstrIter()
			(
				name, 
				alpha, 
				dict
			)
	 );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

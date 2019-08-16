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

#include "ConstantBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "incompleteGammaFunction.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ConstantBreak, 0);
    addToRunTimeSelectionTable(breakupLaw, ConstantBreak, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::ConstantBreak::ConstantBreak
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& dict
)
:
    breakupLaw(name, alpha, dict),
    dictConstantBreak_
    (
        dict.subDict("ConstantBreakCoeffs")
    ),

    // Luo and Svendsen model parameters
    rhoa_(dictConstantBreak_.lookup("rhoa")),
    rhob_(dictConstantBreak_.lookup("rhob")),
    sigma_(dictConstantBreak_.lookup("sigma")),
    beta_(readScalar(dictConstantBreak_.lookup("beta"))),
    k_(dictConstantBreak_.lookup("k"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ConstantBreak::~ConstantBreak()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::ConstantBreak::rhoa() const
{
    return rhoa_;
}

Foam::dimensionedScalar Foam::ConstantBreak::rhob() const
{
    return rhob_;
}

Foam::dimensionedScalar Foam::ConstantBreak::breakupRate(volScalarField& vf, scalar& d1, scalar& d2) const


{


  dimensionedScalar result("SMALL",dimMass/pow3(dimLength)/dimTime, SMALL); //Rate


	  result = dimensionedScalar
				  (
						"value",
						dimMass/pow3(dimLength)/dimTime,
						0.02
				  );

  return result ;


}

Foam::dimensionedScalar Foam::ConstantBreak::breakupFrequency(volScalarField& vf, scalar& d, scalar& fbv) const
{

  dimensionedScalar result("SMALL",pow(dimTime,-1), SMALL); //Frequency



  return result ;

}



// ************************************************************************* //

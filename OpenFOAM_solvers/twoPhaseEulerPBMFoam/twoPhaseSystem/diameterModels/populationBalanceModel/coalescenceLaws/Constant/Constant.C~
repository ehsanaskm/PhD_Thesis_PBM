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

#include "Constant.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "incompleteGammaFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(coalescenceLaw, Constant, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::Constant::Constant
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& dict
)
:
    coalescenceLaw(name, alpha, dict),
    dictConstant_
    (
        dict.subDict("ConstantCoeffs")
    ),
    
    // Luo and Svendsen model parameters
    rhoa_(dictConstant_.lookup("rhoa")),
    rhob_(dictConstant_.lookup("rhob")),
    sigma_(dictConstant_.lookup("sigma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Constant::~Constant()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::Constant::rhoa() const
{
    return rhoa_; 
}

Foam::dimensionedScalar Foam::Constant::rhob() const
{
    return rhob_;
}

Foam::dimensionedScalar Foam::Constant::coalescenceRate(scalar& di, scalar& dj, volScalarField& epsf) const
{


       dimensionedScalar result("coalRate",dimMass/pow3(dimLength)/dimTime,0.0); 
	  result = dimensionedScalar
				  (
						"value",
						dimMass/pow3(dimLength)/dimTime,
						1.0
				  );

     
   return result ;

}

// ************************************************************************* //

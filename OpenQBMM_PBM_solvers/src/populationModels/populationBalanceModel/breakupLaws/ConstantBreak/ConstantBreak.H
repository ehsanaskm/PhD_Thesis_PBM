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
    ConstantBreak

Description
    Constant breakage kernel
    [Ref] Daniele L. Marchisio et al. "Quadrature method of moents for aggregation-breakage processes",
			                         Journal of Colloid and Interface Science, Vol. 258, pp 322-334 (2003)

SourceFiles
    ConstantBreak.C

Author
     Ehsan Askari,  ehsan.askari@usherbrooke.ca  
\*---------------------------------------------------------------------------*/

#ifndef ConstantBreak_H
#define ConstantBreak_H

#include "breakupLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bubbles breakup Declaration
\*---------------------------------------------------------------------------*/

class ConstantBreak
:
    public breakupLaw
{
    // Private data

			//- subDictionary Luo and Svendsen's model
		  dictionary dictConstantBreak_;


  
         //- Gas density
        dimensionedScalar rhoa_;

        //- Liquid density
        dimensionedScalar rhob_;

         //- Surface tension liquid-gas
        dimensionedScalar sigma_;

         //- Luo and Svendsen's model parameters
		  scalar beta_;
		  dimensionedScalar k_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ConstantBreak(const ConstantBreak&);

        //- Disallow default bitwise assignment
        void operator=(const ConstantBreak&);


public:

    //- Runtime type information
    TypeName("ConstantBreak");

    // Static data members


    // Constructors

        //- Construct from dictionary
        ConstantBreak
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& dict
        );


    // Destructor

        virtual ~ConstantBreak();


    // Member Functions

        //- Return bubbles breakup rate
        
        virtual dimensionedScalar breakupRate
        (
				volScalarField&,
				scalar&,
				scalar&

		  ) const;


		  //- Return bubbles breakup frequency
		  virtual dimensionedScalar breakupFrequency
                  (
				volScalarField&,
				scalar&,
				scalar&
		  ) const;

		  //- Return gas density rhoa as scalar
        virtual dimensionedScalar rhoa() const;

		  //- Return liquid density rhob as scalar
        virtual dimensionedScalar rhob() const;



        //- If the breakup models needs population models
        virtual bool populationModelNeeded() const
        {
            return true;
        };

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

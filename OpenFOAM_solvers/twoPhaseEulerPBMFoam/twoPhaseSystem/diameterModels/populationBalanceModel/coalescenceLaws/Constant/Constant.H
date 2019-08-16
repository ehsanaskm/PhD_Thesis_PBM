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
    Constant

Description
    Constant aggregation kernel
    [Ref] Daniele L. Marchisio et al. "Quadrature method of moents for aggregation-breakage processes",
			                         Journal of Colloid and Interface Science, Vol. 258, pp 322-334 (2003)

SourceFiles
    Constant.C

Author
    Ehsan Askari,  ehsan.askari@usherbrooke.ca  
	 

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "coalescenceLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bubbles breakup Declaration
\*---------------------------------------------------------------------------*/

class Constant
:
    public coalescenceLaw
{
    // Private data

			//- subDictionary Luo and Svendsen's model
		  dictionary dictConstant_;
    
         //- Gas density
        dimensionedScalar rhoa_;
        
        //- Liquid density
        dimensionedScalar rhob_;
        
        //- Surface tension liquid-gas
        dimensionedScalar sigma_;
        
  


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Constant(const Constant&);

        //- Disallow default bitwise assignment
        void operator=(const Constant&);


public:

    //- Runtime type information
    TypeName("Constant");

    // Static data members


    // Constructors

        //- Construct from dictionary
        Constant
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& dict
        );


    // Destructor

        virtual ~Constant();


    // Member Functions

        //- Return bubbles coalescence rate 
        virtual dimensionedScalar coalescenceRate
        (
				scalar&, 
				scalar&,
				volScalarField&
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

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
    LuoAndSvendsenCoal

Description
    Luo and Svendsen bubbles coalescence model
    [Ref] H. Luo and F. Svendsen "Theoretical model for drop and bubble
			                         breakup in turbulent dispersions",
			                         AIChE Journal, Vol. 42, pp 1225-1233 (1996)

SourceFiles
    LuoAndSvendsenCoal.C

Author
    brahim selma, consultant 
	 

\*---------------------------------------------------------------------------*/

#ifndef LuoAndSvendsenCoal_H
#define LuoAndSvendsenCoal_H

#include "coalescenceLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bubbles breakup Declaration
\*---------------------------------------------------------------------------*/

class LuoAndSvendsenCoal
:
    public coalescenceLaw
{
    // Private data

			//- subDictionary Luo and Svendsen's model
		  dictionary dictLuoAndSvendsenCoal_;
    
         //- Gas density
        dimensionedScalar rhoa_;
        
        //- Liquid density
        dimensionedScalar rhob_;
        
        //- Surface tension liquid-gas
        dimensionedScalar sigma_;
        
  


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LuoAndSvendsenCoal(const LuoAndSvendsenCoal&);

        //- Disallow default bitwise assignment
        void operator=(const LuoAndSvendsenCoal&);


public:

    //- Runtime type information
    TypeName("LuoAndSvendsen");

    // Static data members


    // Constructors

        //- Construct from dictionary
        LuoAndSvendsenCoal
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& dict
        );


    // Destructor

        virtual ~LuoAndSvendsenCoal();


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

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
    classMethod
     
sourceFiles
    classMethod.H

Author
    Ehsan Askari, M.Sc
    ehsan.askari@usherbrooke.ca

\*---------------------------------------------------------------------------*/

#include "classMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "populationBalanceModel.H"
#include "bound.H"
#include "fvm.H"
//Added by Ehsan
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(classMethod, 0);
    addToRunTimeSelectionTable(populationBalanceReturn, classMethod, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
classMethod::classMethod
(
	const word& name,
	populationBalanceModel& populationBalanceModel
)
:
	populationBalanceReturn(name, populationBalanceModel),
	populationBalanceModel_(populationBalanceModel),
	alpha_(populationBalanceModel.alpha()),
	
	dict_
    (
        populationBalanceModel_.subDict("classMethodCoeffs")
    ),

	dsauter
	(
		IOobject
		(
			"dsauter",
			alpha_.time().timeName(),
			alpha_.db(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		alpha_.mesh()
	),

	// Sauter diameter d32
	d32_
	(
		IOobject
		(
			"d32",
			alpha_.time().timeName(),
			alpha_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
         	alpha_.mesh(),
		dimensionedScalar("SMALL", dimLength, SMALL)
	),

       


	ResidualAlphaForDsauter_(dict_.lookup("ResidualAlphaForDsauter")), 
        ResidualAlphaForAdjust_(dict_.lookup("ResidualAlphaForAdjust")), 
        ResidualAlphaForCorrect_(dict_.lookup("ResidualAlphaForCorrect")), 
        dm_(dict_.lookup("dm")),    
	Nc_(readLabel(dict_.lookup("Nc"))),
	h0_(dict_.lookup("h0")),
	hf_(dict_.lookup("hf")),
	beta1_(dict_.lookup("beta1")),
	beta2_(dict_.lookup("beta2")),
	betaPB_(dict_.lookup("betaPB")),
	k1_(dict_.lookup("k1")),
	sigma_(dict_.lookup("sigma")),
	breakCoeff_(dict_.lookup("breakCoeff")),
	coalCoeff_(dict_.lookup("coalCoeff")),
	S_(dict_.lookup("S")),
	dMin_(dict_.lookup("dMin")),
	dMax_(dict_.lookup("dMax")),
	maxIters_(readLabel(dict_.lookup("maxIters"))),
	loopTolerance_(readScalar(dict_.lookup("loopTolerance")))


	
{
	Info << "Selecting method for population balance model" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

classMethod::~classMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar classMethod::funcTest(scalar& x, scalar& y,	scalar& z) const
{
	 if ( y >= x && y <= z)
    {   
      return 1.0;
    }
	 else
    {
      return 0.0;
    }
}

scalar  classMethod::kronecker(int& i, int& j) const
{
  if (i == j) return 1.0;
   
  else return 0.0;
} 

void classMethod::breakupKernel(PtrList<volScalarField>& S, PtrList<volScalarField>& f, const volScalarField& alpha) const
{
	scalar d_i, x_i, d_i1, x_i1;
   scalar  d_j, x_j, d_i_1, x_i_1;
   scalar *diam = new scalar[Nc_+ 1];
   scalar fbv25 = 0.25, fbv50 = 0.50, fbv75 = 0.75; 
  diam[0] = dm_.value()/10.0;
   
   for(int II=1;II<=Nc_;II++)
    {
      diam[II]= Foam::pow(S_.value(), scalar(II-Nc_/2)/3)*dm_.value(); 
    } 

   volScalarField vf = 1.-alpha;                      
   
	for (int i = 1; i<= Nc_-1; i++)
	{
			d_i = diam[i];
			d_i_1 = diam[i-1];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;
			x_i_1 =  M_PI*pow(d_i_1, 3.0)/6.0;
	              
			for(int j = i+1 ; j <= Nc_ ; j++)
			{	      
				d_j = diam[j];
				x_j =  M_PI*pow(d_j, 3.0)/6.0;
				         
				S[i] +=   populationBalanceModel_.rhoa().value()
									* populationBalanceModel_.breakupModel().breakupRate(vf, d_j, d_i)
                           * f[j]*f[j]*(x_i-x_i_1)*sqr(alpha)
                           * x_i/(x_j+SMALL)/(x_j+SMALL)/2.; 	      
			}
	}
	
	for (int i = 0; i<= Nc_-1; i++)
	{
			d_i = diam[i];
			d_i1 = diam[i+1];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;
			x_i1 =  M_PI*pow(d_i1, 3.0)/6.0;
	              
			for(int j = i+1 ; j <= Nc_; j++)
			{
				d_j = diam[j];
				x_j = M_PI*pow(d_j, 3.0)/6.0;	            
				S[i] += populationBalanceModel_.rhoa().value()
				                  * populationBalanceModel_.breakupModel().breakupRate(vf, d_j, d_i)
                              * f[j]*f[j]*(x_i1-x_i)*sqr(alpha)
                              * x_i/(x_j+SMALL)/(x_j+SMALL)/2.;       
	      }
	}
	
	for (int i = 1; i<= Nc_; i++)
	{
			d_i = diam[i];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;	              
			S[i] -= populationBalanceModel_.rhoa()
			       *( 
			            populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv25)
                   + populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv50)
                   + populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv75)
                 )
                 * f[i]*f[i]*sqr(alpha)
                 * x_i*x_i/(x_i+SMALL)/(x_i+SMALL)/4.; 	 
	}
	
	return;
}

void classMethod::coalescenceKernel
(
	PtrList<volScalarField>& S, 
	PtrList<volScalarField>& f, 
	const volScalarField& alpha,
	const volScalarField& epsilon
) 
const
{
	scalar *diam=new scalar[Nc_+1];
	scalar d_i, d_i1, d_i_1, x_i, x_i1;
	scalar  d_j, d_k, x_k,  x_j, x_i_1, v;
	scalar d_0, d_1, d_n, d_n_1, x_0, x_1, x_n, x_n_1;
	
	diam[0] = dm_.value()/10.0; 

	for(int II=1;II<=Nc_;II++)
	{
      diam[II]= Foam::pow(S_.value(), scalar(II-Nc_/2)/3)*dm_.value(); 
	}

	
      d_0 = diam[0];
      d_1 = diam[1];
      d_n = diam[Nc_];
      d_n_1 = diam[Nc_-1];
      x_0 =  M_PI*::pow(d_0, 3.0)/6.0;
      x_1 =  M_PI*::pow(d_1, 3.0)/6.0;
      x_n =  M_PI*::pow(d_n, 3.0)/6.0;
      x_n_1 =  M_PI*::pow(d_n_1, 3.0)/6.0;
		volScalarField epsf = epsilon; 
	
	for (int ki = 0; ki<=Nc_; ki++)
	{ 
			d_k = diam[ki];
			x_k =  M_PI*::pow(d_k, 3.0)/6.0;
      
			for(int  j = ki ; j <= Nc_ ; j++)
			{
				d_j = diam[j];
				x_j = M_PI*Foam::pow(d_j, 3.0)/6.0;
				v = x_j+x_k;	     
				
				S[0] += populationBalanceModel_.rhoa().value()
				                  * (funcTest(x_0,v,x_1)*(x_1-v)/(x_1-x_0 + SMALL))
                              * (1.0-kronecker(ki,j)/2.)
                              * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                              * sqr(alpha)*x_0/(x_j+SMALL)/(x_k+SMALL)*f[ki]*f[j] ;
	      
				S[Nc_] += populationBalanceModel_.rhoa().value()
				                  * (funcTest(x_n_1,v,x_n)*(v-x_n_1)/(x_n-x_n_1+SMALL))
                              * (1.0-kronecker(ki,j)/2.)
                              * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                              * sqr(alpha)*x_n/(x_j+SMALL)/(x_k+SMALL)*f[ki]*f[j] ;
			}
      
			S[0] -= populationBalanceModel_.rhoa().value()
			             * populationBalanceModel_.coalescenceModel().coalescenceRate(d_0,d_k,epsf)
			             * f[ki]*f[0]
                      * sqr(alpha)*x_0/(x_k+SMALL)/(x_0+SMALL);
	}
		
	for (int i = 1; i<= Nc_-1; i++)
	{
			d_i = diam[i];
			d_i_1 = diam[i-1];
			d_i1 = diam[i+1];
			x_i =  M_PI*::pow(d_i, 3.0)/6.0;
			x_i1 =  M_PI*::pow(d_i1, 3.0)/6.0;
			x_i_1 =  M_PI*::pow(d_i_1, 3.0)/6.0;
	  
			for(int  k = 0 ; k <= Nc_ ; k++)
			{
				d_k = diam[k];
				x_k =  M_PI*::pow(d_k, 3.0)/6.0;
				for(int j = k ; j <= Nc_ ; j++)
				{  
					d_j = diam[j];
					x_j =  M_PI*::pow(d_j, 3.0)/6.0;
					v = x_j+x_k;
	 
					S[i] += populationBalanceModel_.rhoa().value()
					                  * ((funcTest(x_i_1, v, x_i)*(v-x_i_1)/(x_i-x_i_1+SMALL))
                                 + (funcTest(x_i,v,x_i1)*(x_i1-v)/(x_i1-x_i+SMALL)))
                                 * (1.0-kronecker(k,j)/2.)
                                 * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                                 * f[j]*f[k]*sqr(alpha)*x_i/(x_j+SMALL)/(x_k+SMALL); 
				}

				   S[i] -= populationBalanceModel_.rhoa().value()
				                     * populationBalanceModel_.coalescenceModel().coalescenceRate(d_i,d_k,epsf)
				                     * f[k]*f[i]*sqr(alpha)*x_i/(x_k+SMALL)/(x_i+SMALL); 
			}
	} 	
		
	return ; 	
}

tmp<volScalarField> classMethod::SauterDiameter(PtrList<volScalarField>& f, const volScalarField& alpha) const
{ 
  const fvMesh& mesh = alpha_.mesh();
	
  label N=lrint(Nc_) - 1;  
  scalar *diameter=new scalar[N+1];



  volScalarField sum_v
  (
			IOobject
			(
				"sum_v",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
                         alpha
			//dsauter
  );
  
  volScalarField sum_s
  (
			IOobject
			(
				"sum_s",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
                         alpha
			//dsauter
  );
  
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
			  alpha
                         // dsauter
                            
  );

        
  for(label II=0; II<=N; II++)
    {
      diameter[II] = Foam::pow(S_.value(), scalar(II-N/2)/3)*dm_.value(); 
   //   Info<<"Diameter\n"<<diameter[II]<<"\n"<<endl;
    }

					    
		
  		for(label i=0;i<=N;i++)
  		{


	   		forAll(f[i], cI)
				{    
 		   			if(f[i][cI] < 0.0) f[i][cI] = 0.0;      
		   			if(f[i][cI] > 1.0) f[i][cI] = 1.0;
				}
	   	sum_v += f[i];
	   	sum_s += (f[i]/diameter[i]);
   		} 


  // partial result of Sauter diameter d32
    result = sum_v/(sum_s + SMALL);

   


     
  forAll(result, cellI)
  {
	  if(alpha[cellI] > ResidualAlphaForDsauter_.value())
	  {
			//if (result[cellI] > 1.e-6 )
			//{
				result[cellI] = max(dMin_.value(), min(result[cellI], dMax_.value()));		    
			//}
	  }
	  else 
	  {
			  // result[cellI] = diameter[0];
                           result[cellI] = dMin_.value();	
	  }
  } 
	


  return tmp<volScalarField>
  (
		new volScalarField
		(
			IOobject
			(
				"d32",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			result//*
			//dimensionedScalar
			//(
			  // "one", 
			   //dimLength, 
			   //1.
			//)
		)
  );



} 



void classMethod::adjust(PtrList<volScalarField>& S, PtrList<volScalarField>& f, const volScalarField& alpha) const
{

	for(int j=1; j<=Nc_; j++)
   {	  
      forAll(S[j], II)
		{ 	
			S[j][II]=0; 
		}
   }


	for(int i=1; i<=Nc_; i++)
   {
      forAll(f[i], I)
		{
			if(alpha[I] > ResidualAlphaForAdjust_.value())    
			{
				if(f[i][I]<SMALL)
				
				{
				
					f[i][I]=SMALL;
				}
				else if(f[i][I]>1)
				{
					f[i][I]=1.;
				}
			}    
			else 
                        {
                         f[i][I]=1e-15;
                         f[0][I]=1;
                        }
         
		 }	  
   }

}

 
// *******************************************************************************************************************************

// calculacte and correct variables 
void classMethod::correct()
{
  const fvMesh& mesh = alpha_.mesh();

 PtrList<volScalarField> source(2*Nc_); 
  PtrList<volScalarField> f(2*Nc_); 

  
  volScalarField source_ini =
    mesh.objectRegistry::lookupObject<volScalarField>("source_ini");  

  volScalarField f_ini = 
    mesh.objectRegistry::lookupObject<volScalarField>("f_ini");  



  for(label j=0;j<=Nc_+1;j++)
	 {
		 word sourceName = "Sb_" + Foam::name(j); 
		 word fName = "f_" + Foam::name(j); 
                 
		 
		 source.set
		 ( j,
	      volScalarField 
	      ( 
	         IOobject 
	         ( 
				    sourceName, 
					 mesh.time().timeName(),
					 mesh, 
					 IOobject::NO_READ, 
				         IOobject::AUTO_WRITE 
			   ),
	         source_ini
	      ) 
	    );   
	 
		 f.set
		 ( j,
	      volScalarField 
	      ( 
	         IOobject 
	         ( 
				    fName, 
					 mesh.time().timeName(),
					 mesh, 
					 IOobject::NO_READ, 
				         IOobject::AUTO_WRITE 
			   ),
                f_ini
	      ) 
	    );


	 }



//////////////////////////////////////////////////////////////////////////////////////////////////////////


  volScalarField alpha = 
    mesh.objectRegistry::lookupObject<volScalarField>("alpha.air");
  
  volScalarField epsilon = 
   mesh.objectRegistry::lookupObject<volScalarField>("epsilon.water");


  surfaceScalarField phia = 
    mesh.objectRegistry::lookupObject<surfaceScalarField>("phi.air");  



		
		classMethod::adjust(source,f,alpha);
	      
                     

            
 		for(label i=0;i<Nc_;i++)
		{
			// update bubbles breakage/coalescence

			classMethod::coalescenceKernel(source,f,alpha,epsilon);
			classMethod::breakupKernel(source,f,alpha) ;
			volScalarField Sb = source[i]/
			                           populationBalanceModel_.breakupModel().rhoa(); 


			
			if(i != (Nc_-1)/2)
                        
			{
				word fScheme("div(phia,fi)");
				fvScalarMatrix fiEqn
				(
						fvm::ddt(alpha, f[i])	
					 + fvm::div(fvc::flux(phia, alpha, fScheme), f[i], fScheme)	 
					 - fvm::Sp(fvc::div(phia), f[i])  
				);
				solve(fiEqn == Sb);
                          

			}
		}





		classMethod::adjust(source,f,alpha);
   
   // bounding f 
   volScalarField alphai=f[0];	
    
	for(int i=1; i<Nc_; i++)
	{
		if(i!=(Nc_-1)/2) alphai = alphai + f[i];
              		
	}

	f[(Nc_-1)/2] = 1. - alphai;        

  
	f[(Nc_-1)/2] = max(f[(Nc_-1)/2], scalar(0));        
	f[(Nc_-1)/2] = min(f[(Nc_-1)/2], scalar(1));       

	alphai = alphai + f[(Nc_-1)/2];
        

	
	for(int i=0; i<Nc_; i++)
	{		  
		f[i]=f[i]/(alphai + SMALL);
	}


	forAll(f[(Nc_-1)/2], cellI)
	{
		if(alpha[cellI] < ResidualAlphaForCorrect_.value())
		{
			f[(Nc_-1)/2][cellI]=0;
		}
	}


                           
                                

	classMethod::adjust(source,f,alpha);

       
        
      
	// update Sauter diameter d32, relax and correct BC
        dsauter = SauterDiameter(f,alpha);




	d32_ = SauterDiameter(f,alpha)*dimensionedScalar
			(
			   "one", 
			   dimLength, 
			   1.
			); 
                
        dsauter.relax();
        dsauter.correctBoundaryConditions();


        d32_.relax();
        d32_.correctBoundaryConditions();



     
}

// update final calculation 






void classMethod::updateSauterDiameter()
{
   // mesh declaration
 

  

   Info << nl << " *********************************************************************** "<< endl;
   Info << nl << "           Updating Sauter diameter d32 by using " << populationBalanceReturn::name() << endl;
   Info << nl << "           min/avg/max d32: " << min(dsauter).value()*1000
              << "/" << average(dsauter).value()*1000<< "/" << max(dsauter).value()*1000<< "mm." << endl;
   Info << nl << " *********************************************************************** "<< endl;
	Info << " " << endl; 
}


} // end of namespace
// **************************************************************************************************************************** //

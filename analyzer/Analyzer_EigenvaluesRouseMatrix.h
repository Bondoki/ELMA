/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Ron Dockhorn)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef Analyzer_EigenvaluesRouseMatrix_H
#define Analyzer_EigenvaluesRouseMatrix_H

#include <vector>
#include <string>
#include <utility>      // std::pair
#include <map>
#include <vector>
#include <algorithm>    // std::find_if
#include <iostream>
#include <functional>
#include <queue>
#include <cmath> //sqrt

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;


#include <LeMonADE/utility/Vector3D.h>

#include "StatisticMoment.h"

template<class IngredientsType>
class Analyzer_EigenvaluesRouseMatrix:public AbstractAnalyzer
{
public:
	Analyzer_EigenvaluesRouseMatrix(const IngredientsType& ing, std::string dstDir_);


	virtual ~Analyzer_EigenvaluesRouseMatrix(){

	};

	//typedef typename IngredientsType::molecules_type molecules_type;
	const typename IngredientsType::molecules_type& molecules ;

	const IngredientsType& getIngredients() const {return ingredients;}

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();


private:

	const IngredientsType& ingredients;


	uint32_t counterFrames;

	std::string filename;
	std::string dstdir;

	// RNG
	RandomNumberGenerators rng;

	StatisticMoment Statistic_Rg2;
};




/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
Analyzer_EigenvaluesRouseMatrix<IngredientsType>::Analyzer_EigenvaluesRouseMatrix(const IngredientsType& ing, std::string dstDir_)
:ingredients(ing), molecules(ing.getMolecules()),  dstdir(dstDir_)
 {
	counterFrames = 1;
	Statistic_Rg2.clear();

 }

template<class IngredientsType>
void Analyzer_EigenvaluesRouseMatrix<IngredientsType>::initialize()
{
}

template<class IngredientsType>
bool Analyzer_EigenvaluesRouseMatrix<IngredientsType>::execute()
{


	 MatrixXd A(molecules.size(),molecules.size());

		  for(int row = 0; row < molecules.size(); row++)
			  for(int col = 0; col < molecules.size(); col++)
			  {
				  if(molecules.areConnected(row, col) == true)
				  {
					  A(row, col) = -1;
				  }
				  else
				  {
					  A(row, col) = 0;
				  }

				  if(row == col)
				  {
					  A(row, col) = molecules.getNumLinks(row);
				  }
			  }

		 // std::cout << "Here is the matrix m:\n" << A << std::endl;

		   SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
		   if (eigensolver.info() != Success) abort();
		   //std::cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
		   //std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
		   //     << "corresponding to these eigenvalues:\n"
		   //     << eigensolver.eigenvectors() << std::endl;

		   double RG2 = 0.0;

		   //first ommited as this is zero
		   for(int i = 1; i < eigensolver.eigenvalues().rows(); i++)
		   {
			   double lambda = eigensolver.eigenvalues()[i];
			   //std::cout << "Consider the " << (i+1) << " eigenvalue, lambda = " << lambda << std::endl;
			   RG2 += 1.0/lambda;
		   }

		   RG2 /= eigensolver.eigenvalues().rows();

		   Statistic_Rg2.AddValue(RG2);

	counterFrames++;
	std::cout << " run : " << counterFrames << std::endl;

}



struct PathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
void Analyzer_EigenvaluesRouseMatrix<IngredientsType>::cleanup()
{

	// print results into a file
	// get the filename and path
	// find the filename without path and extensions
			std::string filenameGeneral= std::string( std::find_if( ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator() ).base(), ingredients.getName().end() );

			std::string::size_type const p(filenameGeneral.find_last_of('.'));
			filenameGeneral = filenameGeneral.substr(0, p);


	// construct a list
	std::vector < std::vector<double> > tmpResultsRG;

	// we have 3 columns and row
	uint32_t columns = 6;
	uint32_t rows = 1;



	// we have columns
	tmpResultsRG.resize(columns);

	// we have rows
	for(int i = 0; i < columns; i++)
		tmpResultsRG[i].resize(rows);

	tmpResultsRG[0][0] = molecules.size();
	tmpResultsRG[1][0] = Statistic_Rg2.ReturnM1();
	tmpResultsRG[2][0] = Statistic_Rg2.ReturnM2();
	tmpResultsRG[3][0] = Statistic_Rg2.ReturnVar();
	tmpResultsRG[4][0] = 2.0*Statistic_Rg2.ReturnSigma()/ std::sqrt(1.0* Statistic_Rg2.ReturnN() );
	tmpResultsRG[5][0] = Statistic_Rg2.ReturnN();


	std::stringstream comment;
	comment << " File produced by analyzer Analyzer_EigenvaluesRouseMatrix" << std::endl
			<< " RG² = 1/N sum _{i=2} ^{N} 1/lambda_i " << std::endl
			<< " with N = " << ingredients.getMolecules().size() << " monomers" << std::endl
			<< std::endl
			<< "N <RG²> <RG²>² Varianz(RG²) Error SampleSize";



	//new filename
	std::string filenameRG = filenameGeneral + "_EigenvaluesRouseMatrix.dat";

	//ResultFormattingTools::writeResultFile(dstdir+"/"+filenameRG, this->ingredients, tmpResultsRG, comment.str());
	ResultFormattingTools::appendToResultFile(dstdir+"/"+filenameRG, tmpResultsRG);

}

#endif /*Analyzer_EigenvaluesRouseMatrix_H*/

/* **************************************************************
 * A simple analyzer example for static properties
 *
 * In this case it calculates the average RgSquared of the molecules in
 * the system. Plain and simple...well...at least plain. Have fun.
 * *************************************************************/

#ifndef ANALYZER_CHAIN_WALKING_RG2_H
#define ANALYZER_CHAIN_WALKING_RG2_H

#include <iostream>
#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>



/*****************************************************************************
 * CLASS DEFINITION (implementation of methods below)
 * **************************************************************************/

template<class IngredientsType>
class Analyzer_ChainWalking_RG2:public AbstractAnalyzer
{
public:

  //constuctor
  Analyzer_ChainWalking_RG2(const IngredientsType& ing, long evalulation_time_);

  //initializes the groups. called explicitly or by Taskmanager::init()
  virtual void initialize();

  //does all the calculations esp in every !mcs
  virtual bool execute();

  //write out your results, in this case to standard output
  virtual void cleanup();

private:



 //holds a reference of the complete system
 const IngredientsType& ingredients;


 //for calculating the average RgSquared and the components
 double sumRg2;
 double sumRg2_x;
 double sumRg2_y;
 double sumRg2_z;

 uint32_t nValues;

 double sumBondLength2;
 uint32_t nValuesBondLength2;

 //only used to make sure you initialize your groups before you do things
 bool initialized;

 long evalulation_time;
};


/*****************************************************************************
 * IMPLEMENTATION OF METHODS
 * **************************************************************************/



/* ****************************************************************************
 * constructor. only initializes some variables
 * ***************************************************************************/
template<class IngredientsType>
Analyzer_ChainWalking_RG2<IngredientsType>::Analyzer_ChainWalking_RG2(const IngredientsType& ing, long evalulation_time_)
 :ingredients(ing),initialized(false),sumRg2(0.0),nValues(0),evalulation_time(evalulation_time_)
{
}

/* **********************************************************************
 * initialize()
 *
 * this is called in the beginning only once - need for setting up your analyzer
 * **********************************************************************/
template<class IngredientsType>
void Analyzer_ChainWalking_RG2<IngredientsType>::initialize()
{
	sumRg2 = 0.0;
	sumRg2_x = 0.0;
	sumRg2_y = 0.0;
	sumRg2_z = 0.0;

	nValues = 0;

	sumBondLength2 = 0.0;
	nValuesBondLength2 = 0;

	//set the initialized tag to true
	initialized=true;

}

/* ***********************************************************************
 * execute()
 * this is where the calculation happens
 * ***********************************************************************/
template<class IngredientsType>
bool Analyzer_ChainWalking_RG2<IngredientsType>::execute()
{
	//check if analyzer have been initialized. if not, exit and explain
	if(initialized==false)
	{
		std::stringstream errormessage;
		errormessage<<"Analyzer_ChainWalking_RG2::execute()...Analyzer_ChainWalking_RG2 not initialized\n"
			    <<"Use Analyzer_ChainWalking_RG2::initialize() or Taskmanager::init()\n";

		throw std::runtime_error(errormessage.str());
	}


	//doing the calculation for Rg2 of the whole molecule in Ingredients

	//wait relaxation time
	if(ingredients.getMolecules().getAge() > evalulation_time)
	{
		std::cout << "SimpleAnalyzer_Rg2.execute() at MCS:" << ingredients.getMolecules().getAge() << std::endl;

		int monomerCounter = 0;


		// The radius of gyration is defined as
		// Rg2 = 1/N * SUM_i=1 (r_i - r_COM)^2 = 1/N^2 * SUM_i=1 SUM_j=i (r_i - r_j)^2
		// The components in the same manner

		double Rg2 = 0.0;
		double Rg2_x = 0.0;
		double Rg2_y = 0.0;
		double Rg2_z = 0.0;

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l= k; l < ingredients.getMolecules().size(); l++)
                 if((ingredients.getMolecules()[l].getAttributeTag()==1) && (ingredients.getMolecules()[k].getAttributeTag()==1))
			{
				Rg2_x += (ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX())*(ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX());
				Rg2_y += (ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY())*(ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY());
				Rg2_z += (ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ())*(ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ());


			}
			
			// only count the monomers which contribute to RG2
			if((ingredients.getMolecules()[l].getAttributeTag()==1))
				monomerCounter++;
		}
		/*if(monomerCounter != ingredients.getMolecules().size())
		{
			throw std::runtime_error("invalid number of monomers in Rg2-calculation");
		}

		Rg2_x /= 1.0*(ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_y /= 1.0*(ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_z /= 1.0*(ingredients.getMolecules().size()*ingredients.getMolecules().size());
        */

        Rg2_x /= 1.0*(monomerCounter*monomerCounter);
		Rg2_y /= 1.0*(monomerCounter*monomerCounter);
		Rg2_z /= 1.0*(monomerCounter*monomerCounter);

		Rg2 = Rg2_x+Rg2_y+Rg2_z;


		// add value to the average
		sumRg2   += Rg2;
		sumRg2_x += Rg2_x;
		sumRg2_y += Rg2_y;
		sumRg2_z += Rg2_z;

		// increase the counter for average calculation
		nValues++;


		// calculate the bond length of the connected structure
		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l = 0; l < ingredients.getMolecules().getNumLinks(k); l++)
                	{
				if( k < ingredients.getMolecules().getNeighborIdx(k, l))
					if((ingredients.getMolecules()[k].getAttributeTag()==1) && (ingredients.getMolecules()[ingredients.getMolecules().getNeighborIdx(k, l)].getAttributeTag()==1))
					{
					double BondLength = (ingredients.getMolecules()[k]-(ingredients.getMolecules()[ingredients.getMolecules().getNeighborIdx(k, l)])).getLength();
					sumBondLength2 += BondLength*BondLength;
					nValuesBondLength2++;
					}
			}
		}
	}

}

/****************************************************************************
 * cleanup
 * this is where you get the final value and write it to file, std output, or
 * whatever
 * this is called in the end only once - need for write output to a file or stdout
 * *************************************************************************/
template<class IngredientsType>
void Analyzer_ChainWalking_RG2<IngredientsType>::cleanup()
{
	// calculate the real average
	sumRg2   = sumRg2  /(double (nValues));
	sumRg2_x = sumRg2_x/(double (nValues));
	sumRg2_y = sumRg2_y/(double (nValues));
	sumRg2_z = sumRg2_z/(double (nValues));

	// print results to stdout
	std::cout<<"Average Rg2   : <Rg2  > = " << sumRg2   << std::endl;
	std::cout<<"Average Rg2_x : <Rg2_x> = " << sumRg2_x << std::endl;
	std::cout<<"Average Rg2_y : <Rg2_y> = " << sumRg2_y << std::endl;
	std::cout<<"Average Rg2_z : <Rg2_z> = " << sumRg2_z << std::endl;

	// calculate the real average of bond length
	sumBondLength2   = sumBondLength2  /(double (nValuesBondLength2));
	std::cout<<"Average (bond length)^2   : < b^2 > = " << sumBondLength2   << std::endl;
	// print results into a file

	// get the filename and path
	std::string filenameGeneral=ingredients.getName();
	// delete the .bfm in the name
	filenameGeneral.erase (ingredients.getName().length()-4, ingredients.getName().length());

	// construct a list
	std::vector < std::vector<double> > tmpResultsRg2;

	// we have 4 columns and 1 row
	uint32_t columns = 9;
	uint32_t rows = 1;

	// we have columns
	tmpResultsRg2.resize(columns);

	// we have rows
	for(int i = 0; i < columns; i++)
		tmpResultsRg2[i].resize(rows);

	// fill the list

	tmpResultsRg2[0][0]=sumRg2;
	tmpResultsRg2[1][0]=sumRg2_x;
	tmpResultsRg2[2][0]=sumRg2_y;
	tmpResultsRg2[3][0]=sumRg2_z;
	tmpResultsRg2[4][0]=sumRg2/sumBondLength2;
	tmpResultsRg2[5][0]=sumRg2_x/sumBondLength2;
	tmpResultsRg2[6][0]=sumRg2_y/sumBondLength2;
	tmpResultsRg2[7][0]=sumRg2_z/sumBondLength2;
	tmpResultsRg2[8][0]=sumBondLength2;


	std::stringstream comment;
	comment << "File produced by analyzer Analyzer_ChainWalking_RG2" << std::endl
			<< "Radius of Gyration Rg2 with " << ingredients.getMolecules().size() << " monomers" << std::endl
			<< "Average squared bond length <b^2>=" << sumBondLength2  << std::endl
			<< std::endl
			<< "<Rg2> <Rg2_x> <Rg2_y> <Rg2_z> <Rg2>/<b^2> <Rg2_x>/<b^2> <Rg2_y>/<b^2> <Rg2_z>/<b^2> <b^2>";



	//new filename
	std::string filenameRg2 = filenameGeneral + "_Rg2.dat";

	ResultFormattingTools::writeResultFile(filenameRg2, this->ingredients, tmpResultsRg2, comment.str());

}



#endif /*ANALYZER_CREATOR_SLOW_GROWTH_RG2_H*/

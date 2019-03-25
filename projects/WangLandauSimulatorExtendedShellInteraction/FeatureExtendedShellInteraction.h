/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
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


#ifndef FEATURE_EXTENDED_SHELL_INTERACTION_H
#define FEATURE_EXTENDED_SHELL_INTERACTION_H

/**
 * @file
 * @date 2018/09/21
 * @author Ron Dockhorn
 * @brief Definition and implementation of class template FeatureExtendedShellInteraction
**/

#include <algorithm>
#include <vector>
#include <iterator>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>


#include "FeatureExtendedShellInteractionReadWrite.h"
#include "FeatureExcludedVolumeSingleNode.h"
#include "Histogram1D.h"
#include "HistogramGeneralStatistik1D.h"

/**
 * @class FeatureExtendedShellInteraction
 * @brief Provides interaction of monomers on distances d<=sqrt(6) for standard BFM
 *
 * @tparam FeatureLatticeType Underlying lattice feature, e.g. FeatureLattice or 
 * FeatureLatticePowerOfTwo (template template parameter)
 *
 * @details
 * The interaction energy can be set for pairs of monomer-types A,B, where
 * the type can be any integer between 1 and 255. 
 * The feature automatically adds FeatureExcludedVolumeSc<FeatureLatticeType<uint8_t>
 * to the system. Given an energy E in kT between  two types, the interaction potential 
 * as a function of the distance d is:
 * - inf d<2 (implicitly through excluded volume)
 * - 4*E d=2
 * - 2*E d=sqrt(5)
 * - 1*E d=sqrt(6)
 * - 0   d>sqrt(6)
 * .
 * Usage: In the feature list defining Ingredients use this feature as
 * FeatureExtendedShellInteraction<FeatureLattice> (arbitrary lattices), or as
 * FeatureExtendedShellInteraction<FeatureLatticePowerOfTwo> (2**n lattices)
 * The feature adds the bfm-file command !nn_interaction A B E
 * for monomers of types A B with interaction energy of E in kT.
**/

template<template<typename> class FeatureLatticeType>
class FeatureExtendedShellInteraction:public Feature
{

private:
  //! Type for the underlying lattice, used as template parameter for FeatureLatticeType<...>
  typedef uint8_t lattice_value_type;

  //! Interaction energies between monomer types. Max. type=255 given by max(uint8_t)=255 
  double interactionTable[256][256];
  
  //! Lookup table for exp(-interactionTable[a][b])
  double probabilityLookup[256][256];

  //! Returns this feature's factor for the acceptance probability for the given Monte Carlo move
  template<class IngredientsType>
  double calculateAcceptanceProbability(const IngredientsType& ingredients, 
					const MoveLocalSc& move) const;

  //! Occupies the lattice with the attribute tags of all monomers
  template<class IngredientsType>
  void fillLattice(IngredientsType& ingredients);

  //! Access to array probabilityLookup with extra checks in Debug mode
  double getProbabilityFactor(int32_t typeA,int32_t typeB) const;


  //! Extended Shell for Interaction
  std::vector<VectorInt3> ExtendedShell;

  //! Extended Shell for Interaction
  std::vector<VectorInt3> ExtendedShellMinusXDirNew;
  std::vector<VectorInt3> ExtendedShellMinusXDirOld;

  std::vector<VectorInt3> ExtendedShellMinusYDirNew;
  std::vector<VectorInt3> ExtendedShellMinusYDirOld;

  std::vector<VectorInt3> ExtendedShellMinusZDirNew;
  std::vector<VectorInt3> ExtendedShellMinusZDirOld;

  std::vector<VectorInt3> ExtendedShellPlusXDirNew;
  std::vector<VectorInt3> ExtendedShellPlusXDirOld;

  std::vector<VectorInt3> ExtendedShellPlusYDirNew;
  std::vector<VectorInt3> ExtendedShellPlusYDirOld;

  std::vector<VectorInt3> ExtendedShellPlusZDirNew;
  std::vector<VectorInt3> ExtendedShellPlusZDirOld;



  //! Extended Shell radius
  double ExtendedShellRadius;

  Histogram1D HG_VisitsEnergyStates;
  HistogramGeneralStatistik1D HG_LnDOS;
  double Energy;
  double diffEnergy;
  double modificationFactor;

  Histogram1D HG_TotalVisitsEnergyStates;

  //! isEnergyInWindow
  bool windowingState;
  double minWin;
  double maxWin;

public:

  FeatureExtendedShellInteraction();
  ~FeatureExtendedShellInteraction(){}

    
  //This feature adds interaction energies, so it requires FeatureBoltzmann
  typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

  //FeatureExcludedVolumeSc needs to be in front, because FeatureExtendedShellInteraction
  //re-initializes the lattice and overwrites what FeatureExcludedVolumeSc has written.
  //FeatureAttributes needs to be in front, because when a monomer is added to the system
  //by a MoveAddMonomerSc, its attribute has to be set before it is written to the lattice.
  typedef LOKI_TYPELIST_2(
      FeatureAttributes,
      FeatureExcludedVolumeSingleNode<FeatureLatticeType<lattice_value_type> >)
    required_features_front;


  //! check for all Monte Carlo moves without special check functions (always true)
  template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;

  //! check for standard sc-BFM local move
  template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move);

  //! check move for bcc-BFM local move. always throws std::runtime_error
  template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveLocalBcc& move) const;

  //! apply function for all Monte Carlo moves without special apply functions (does nothing)
  template<class IngredientsType>
    void applyMove(const IngredientsType& ing, const MoveBase& move){}

  //! apply function for bcc-BFM local move (always throws std::runtime_error)
  template<class IngredientsType>
    void applyMove(const IngredientsType& ing, const MoveLocalBcc& move);

  //for moves of type MoveLocalSc
  template<class IngredientsType>
  void applyMove(IngredientsType& ing, const MoveLocalSc& move);

  //! apply function for adding a monomer in sc-BFM
  template<class IngredientsType>
  void applyMove(IngredientsType& ing, const MoveAddMonomerSc& move);


  template<class IngredientsType>
    	void rejectMove(IngredientsType& ingredients);

  //note: apply function for sc-BFM local move is not necessary, because 
  //job of moving lattice entries is done by the underlying FeatureLatticeType

  //! guarantees that the lattice is properly occupied with monomer attributes
  template<class IngredientsType>
  void synchronize(IngredientsType& ingredients);

  //!adds interaction energy between two types of monomers
  void setExtendedShellInteraction(int32_t typeA,int32_t typeB,double energy);

  //!returns the interaction energy between two types of monomers
  double getExtendedShellInteraction(int32_t typeA,int32_t typeB) const;

  //!sets radius of interaction energy calculation between two types of monomers
  void setExtendedShellInteractionRadius(double radius);

  //!returns the radius of interaction energy calculation between two types of monomers
  double getExtendedShellInteractionRadius() const;

  //!export bfm-file read command !nn_interaction
  template <class IngredientsType>
  void exportRead(FileImport <IngredientsType>& fileReader);

  //!export bfm-file write command !nn_interaction
  template <class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

  HistogramGeneralStatistik1D& modifyHGLnDOS() {return HG_LnDOS;}
  const HistogramGeneralStatistik1D& getHGLnDOS() const {return HG_LnDOS;}

  Histogram1D& modifyVisitsEnergyStates() {return HG_VisitsEnergyStates;}
  const Histogram1D& getVisitsEnergyStates() const {return HG_VisitsEnergyStates;}

  Histogram1D& modifyTotalVisitsEnergyStates() {return HG_TotalVisitsEnergyStates;}
  const Histogram1D& getTotalVisitsEnergyStates() const {return HG_TotalVisitsEnergyStates;}

  double getModificationFactor() const {return modificationFactor;}

  void setModificationFactor(double modificationFactor) {

	  std::cout << "Modification Factor f=" << modificationFactor << std::endl;
	  this->modificationFactor = modificationFactor;

  }

  bool isEnergyInWindow() const {
	  return windowingState;
  }

  void setWindowState(bool isEnergyInWindow, double minWin, double maxWin) {
	  this->windowingState = isEnergyInWindow;
	  this->minWin = minWin;
	  this->maxWin = maxWin;
  }

  double getMaxWin() const {
	  return maxWin;
  }

  double getMinWin() const {
	  return minWin;
  }

  void setMinMaxWin(double minWin, double maxWin) {
	  this->minWin = minWin;
	  this->maxWin = maxWin;
  }

  //! Returns this feature's factor for the energy difference for the given Monte Carlo move
  template<class IngredientsType>
  double calculateInteractionDifference(const IngredientsType& ingredients,
		  const MoveLocalSc& move) const;

  template<class IngredientsType>
  double getInternalEnergyCurrentConfiguration(const IngredientsType& ingredients) const;

};


//////////////////  IMPLEMENTATION OF MEMBERS //////////////////////////////////

/**
 * @brief Constructor
 **/
template<template<typename> class LatticeClassType>
FeatureExtendedShellInteraction<LatticeClassType>::FeatureExtendedShellInteraction()
{
	Energy=0.0;
	diffEnergy=0.0;
	modificationFactor=1.01;
	windowingState=false;

  //initialize the energy and probability lookups with default values
  for(size_t n=0;n<256;n++)
    {
      for(size_t m=0;m<256;m++)
        {
	  interactionTable[m][n]=0.0;
	  probabilityLookup[m][n]=1.0;
        }
    }

  ExtendedShellRadius = 4.0;
  ExtendedShell.clear();

  ExtendedShellMinusXDirNew.clear();
  ExtendedShellMinusYDirNew.clear();
  ExtendedShellMinusZDirNew.clear();
  ExtendedShellPlusXDirNew.clear();
  ExtendedShellPlusYDirNew.clear();
  ExtendedShellPlusZDirNew.clear();

  ExtendedShellMinusXDirOld.clear();
  ExtendedShellMinusYDirOld.clear();
  ExtendedShellMinusZDirOld.clear();
  ExtendedShellPlusXDirOld.clear();
  ExtendedShellPlusYDirOld.clear();
  ExtendedShellPlusZDirOld.clear();
}

/**
 * @details Returns true for all moves other than the ones that have specialized versions of this function.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move other than moves with specialized functions.
 * @return true (always)
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureExtendedShellInteraction<LatticeClassType>::checkMove(const IngredientsType& ingredients,
							 const MoveBase& move) const
{
    return true;
}

/**
 * @details calculates the factor for the acceptance probability of the move
 * arising from the contact interactions and adds it to the move.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalSc
 * @return true (always)
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureExtendedShellInteraction<LatticeClassType>::checkMove(const IngredientsType& ingredients,
							 MoveLocalSc& move)
{
  //add the probability factor coming from this feature, then return true,
  //because the total probability is evaluated by FeatureBoltzmann at the end
  //double prob=calculateAcceptanceProbability(ingredients,move);
  //move.multiplyProbability(prob);
  //return true;

  if(windowingState == false)
  	{
  	diffEnergy = calculateInteractionDifference(ingredients,move);
  	//std::cout << "EnergyOld: " << Energy <<"\t EnergyNew: " << (Energy+diffEnergy)  << " \t dEnergy: " << (diffEnergy) << std::endl;

  	double p=1.0;

  	// p = min(1,g(E_old)/g(E_new))
  	if(HG_LnDOS.getCountAt(Energy) < HG_LnDOS.getCountAt(Energy+diffEnergy))
  		p=std::exp(HG_LnDOS.getCountAt(Energy)-HG_LnDOS.getCountAt(Energy+diffEnergy));

  	//add move probability according to current potentialOfMeanForce
  	move.multiplyProbability(p);
  	}
  	else
  	{
  		diffEnergy = calculateInteractionDifference(ingredients,move);

  		//if( (((Energy+diffEnergy) > minWin-2.0*HG_LnDOS.getBinwidth()) && ((Energy+diffEnergy) < maxWin+2.0*HG_LnDOS.getBinwidth())) && ( (Energy > minWin-2.0*HG_LnDOS.getBinwidth()) && (Energy < maxWin+2.0*HG_LnDOS.getBinwidth()) ))
  		if( ( ((Energy+diffEnergy) > minWin) && ((Energy+diffEnergy) < maxWin)) )
  		{
  			//std::cout << "EnergyOld: " << Energy <<"\t EnergyNew: " << (Energy+diffEnergy)  << " \t dEnergy: " << (diffEnergy) << std::endl;


  			double p=1.0;

  			// p = min(1,g(E_old)/g(E_new))

  			if(HG_LnDOS.getCountAt(Energy) < HG_LnDOS.getCountAt(Energy+diffEnergy))
  				p=std::exp(HG_LnDOS.getCountAt(Energy)-HG_LnDOS.getCountAt(Energy+diffEnergy));


  			move.multiplyProbability(p);
  		}
  		else
  		{
  			move.multiplyProbability(0);
  			return false;
  		}

  	}
  	return true;
}

/**
 * @details Because moves of type MoveLocalBcc must not be used with this
 * feature, this function always throws an exception when called. The function
 * is only implemented for savety purposes.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalBcc
 * @throw std::runtime_error
 * @return false always throws exception before returning
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureExtendedShellInteraction<LatticeClassType>::checkMove(const IngredientsType& ingredients,
							 const MoveLocalBcc& move) const
{
  //throw exception in case someone accidentaly uses a bcc-BFM move with this feature
  std::stringstream errormessage;
  errormessage<<"FeatureExtendedShellInteraction::checkMove(...):\n";
  errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
  throw std::runtime_error(errormessage.str());

  return false;
}


/**
 * @details When a new monomer is added to the system, the lattice sites occupied
 * by this monomer must be filled with the attribute tag. This function takes care
 * of this.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveAddMonomerSc
 * @throw std::runtime_error if attribute tag is not in range [1,255]
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::applyMove(IngredientsType& ing, const MoveLocalSc& move)
{
	//update the book keeping variables in case the move was accepted
	//EnergyOld=EnergyNew;
	Energy += diffEnergy;

	//HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
	//HG_VisitsEnergyStates.addValue(Energy, 1.0);
	//HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);

	if(windowingState == false)
	{
		HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
		HG_VisitsEnergyStates.addValue(Energy, 1.0);
		HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);
		;
	}
	else
	{

		HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
		HG_VisitsEnergyStates.addValue(Energy, 1.0);
		HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);

	}

/*
    //get the position and attribute tag of the monomer to be inserted
    VectorInt3 pos=move.getPosition();
    VectorInt3 dx(1,0,0);
    VectorInt3 dy(0,1,0);
    VectorInt3 dz(0,0,1);
    lattice_value_type type=lattice_value_type(move.getTag());

    //the feature is based on a uint8_t lattice, thus the max type must not
    //exceed the max value of uint8_t (255)
    if(int32_t(type) != move.getTag() || int32_t(type)==0)
      {
	std::stringstream errormessage;
	errormessage<<"FeatureExtendedShellInteraction::applyMove(MoveAddMonomerSc)\n";
	errormessage<<"Trying to add monomer with type "<<int32_t(type)<<">maxType=255\n";
	throw std::runtime_error(errormessage.str());
      }
    
    //update lattice
    ing.setLatticeEntry(pos,type);
    */
}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::rejectMove(IngredientsType& ingredients)//, const MoveLocalSc& move)
{
	//HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
	//HG_VisitsEnergyStates.addValue(Energy, 1.0);
	//HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);


	if(windowingState == false)
	{
			HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
			HG_VisitsEnergyStates.addValue(Energy, 1.0);
			HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);
		;
	}
	else
	{
		HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
		HG_VisitsEnergyStates.addValue(Energy, 1.0);
		HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);

	}
}

/**
 * @details Because moves of type MoveLocalBcc must not be used with this
 * feature, this function always throws an exception when called. The function
 * is only implemented for savety purposes.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalBcc
 * @throw std::runtime_error
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::applyMove(const IngredientsType& ing,
							 const MoveLocalBcc& move)
{
  //throw exception in case someone accidentaly uses a bcc-BFM move with this feature
  std::stringstream errormessage;
  errormessage<<"FeatureExtendedShellInteraction::applyMove(...):\n";
  errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
  throw std::runtime_error(errormessage.str());

}

/**
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::synchronize(IngredientsType& ingredients)
{

	//refill the lattice with attribute tags
	//caution: this overwrites, what is currently written on the lattice
	std::cout << "FeatureExtendedShellInteraction::synchronizing lattice occupation...";
	fillLattice(ingredients);
	std::cout << "done\n";

	//calculate the Extended Shell between r_min = 2 <= r <= r_max (4)

	ExtendedShell.clear();

	ExtendedShellMinusXDirNew.clear();
	ExtendedShellMinusYDirNew.clear();
	ExtendedShellMinusZDirNew.clear();
	ExtendedShellMinusXDirOld.clear();
	ExtendedShellMinusYDirOld.clear();
	ExtendedShellMinusZDirOld.clear();

	ExtendedShellPlusXDirNew.clear();
	ExtendedShellPlusYDirNew.clear();
	ExtendedShellPlusZDirNew.clear();
	ExtendedShellPlusXDirOld.clear();
	ExtendedShellPlusYDirOld.clear();
	ExtendedShellPlusZDirOld.clear();



	if(ExtendedShellRadius > 6.0){
		std::stringstream errormessage;
		errormessage<<"***FeatureExtendedShellInteraction::synchronize***\n";
		errormessage<<"Seriously, you want to use a interaction shell with "<<ExtendedShellRadius<<"?!" <<std::endl;
		errormessage<<"Not with this approach. Use your own implementation :P";
		throw std::runtime_error(errormessage.str());
	}

	uint32_t counter = 0;

	std::vector<VectorInt3> tmpMX;
	std::vector<VectorInt3> tmpPX;
	std::vector<VectorInt3> tmpMY;
	std::vector<VectorInt3> tmpPY;
	std::vector<VectorInt3> tmpMZ;
	std::vector<VectorInt3> tmpPZ;

	tmpMX.clear();
	tmpPX.clear();
	tmpMY.clear();
	tmpPY.clear();
	tmpMZ.clear();
	tmpPZ.clear();


	for(int dx = -10; dx < 11 ; dx++)
		for(int dy = -10; dy < 11 ; dy++)
			for(int dz = -10; dz < 11 ; dz++)
			{
				double length = std::sqrt(dx*dx+dy*dy+dz*dz);

				// due to ex vol the minimal distance can only be rmin=2
				if((dx*dx+dy*dy+dz*dz) >= 4 && (dx*dx+dy*dy+dz*dz) <= ExtendedShellRadius*ExtendedShellRadius)
				{
					counter ++;
					std::cout << counter << " with ( " << dx << " , " << dy << " , " << dz << " ) -> " << length << std::endl;

					VectorInt3 r(dx, dy, dz);
					ExtendedShell.push_back(r);

					// vectors smaller than 3^0.5 can be neglegted due to ex vol
					if( (r+VectorInt3(-1,0,0)).getLength() > 1.74)
						tmpMX.push_back(r+VectorInt3(-1,0,0));

					if( (r+VectorInt3(1,0,0)).getLength() > 1.74)
						tmpPX.push_back(r+VectorInt3(1,0,0));

					if( (r+VectorInt3(0,-1,0)).getLength() > 1.74)
						tmpMY.push_back(r+VectorInt3(0,-1,0));

					if( (r+VectorInt3(0,1,0)).getLength() > 1.74)
						tmpPY.push_back(r+VectorInt3(0,1,0));

					if( (r+VectorInt3(0,0,-1)).getLength() > 1.74)
						tmpMZ.push_back(r+VectorInt3(0,0,-1));

					if( (r+VectorInt3(0,0,1)).getLength() > 1.74)
						tmpPZ.push_back(r+VectorInt3(0,0,1));

				}
			}

	// calculate symmetric difference M = (A\B) U (B\A) = (A U B) \ (A ^ B)
	std::cout << "Calulate symmetric difference set." << std::endl;

	// (A\B)
	for( std::vector<VectorInt3>::iterator it = ExtendedShell.begin(); it != ExtendedShell.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iterMX = std::find( tmpMX.begin(), tmpMX.end(), *it );
		std::vector<VectorInt3>::iterator iterPX = std::find( tmpPX.begin(), tmpPX.end(), *it );
		std::vector<VectorInt3>::iterator iterMY = std::find( tmpMY.begin(), tmpMY.end(), *it );
		std::vector<VectorInt3>::iterator iterPY = std::find( tmpPY.begin(), tmpPY.end(), *it );
		std::vector<VectorInt3>::iterator iterMZ = std::find( tmpMZ.begin(), tmpMZ.end(), *it );
		std::vector<VectorInt3>::iterator iterPZ = std::find( tmpPZ.begin(), tmpPZ.end(), *it );

		// if not found element is only in A  == find returns tmpMX.end()
	    if( iterMX == tmpMX.end() ) {	ExtendedShellMinusXDirOld.push_back(*it); }
	    if( iterPX == tmpPX.end() ) {	ExtendedShellPlusXDirOld.push_back(*it); }

	    if( iterMY == tmpMY.end() ) {	ExtendedShellMinusYDirOld.push_back(*it); }
	    if( iterPY == tmpPY.end() ) {	ExtendedShellPlusYDirOld.push_back(*it); }

	    if( iterMZ == tmpMZ.end() ) {	ExtendedShellMinusZDirOld.push_back(*it); }
	    if( iterPZ == tmpPZ.end() ) {	ExtendedShellPlusZDirOld.push_back(*it); }
	}

	// (B\A)
	for( std::vector<VectorInt3>::iterator it =  tmpMX.begin(); it != tmpMX.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellMinusXDirNew.push_back(*it); }
	}

	for( std::vector<VectorInt3>::iterator it =  tmpPX.begin(); it != tmpPX.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellPlusXDirNew.push_back(*it); }
	}

	for( std::vector<VectorInt3>::iterator it =  tmpMY.begin(); it != tmpMY.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellMinusYDirNew.push_back(*it); }
	}

	for( std::vector<VectorInt3>::iterator it =  tmpPY.begin(); it != tmpPY.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellPlusYDirNew.push_back(*it); }
	}

	for( std::vector<VectorInt3>::iterator it =  tmpMZ.begin(); it != tmpMZ.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellMinusZDirNew.push_back(*it); }
	}

	for( std::vector<VectorInt3>::iterator it =  tmpPZ.begin(); it != tmpPZ.end(); ++it )
	{
		std::vector<VectorInt3>::iterator iter = std::find( ExtendedShell.begin(), ExtendedShell.end(), *it );
		// if not found element is only in B == find returns ExtendedShell.end()
		if( iter == ExtendedShell.end() ) { ExtendedShellPlusZDirNew.push_back(*it); }
	}

	/*
	std::sort(ExtendedShell.begin(), ExtendedShell.end());
	std::sort(tmpMX.begin(), tmpMX.end());
	// this approach does not work
	//set_symmetric_difference(ExtendedShell.begin(), ExtendedShell.end(), tmpMX.begin(), tmpMX.end(), std::inserter(ExtendedShellMinusXDir, ExtendedShellMinusXDir.begin()) );
	//set_symmetric_difference(ExtendedShell.begin(), ExtendedShell.end(), tmpMX.begin(), tmpMX.end(), std::back_inserter(ExtendedShellMinusXDir));
	*/

	for (std::vector<VectorInt3>::iterator i = ExtendedShellMinusXDirOld.begin(); i != ExtendedShellMinusXDirOld.end(); i++) {
		std::cout << *i << " " << std::endl;
	    }

	std::cout << "Added MXNew" << ExtendedShellMinusXDirNew.size() << " position query on -X to the extended interaction shell. " << std::endl;
	std::cout << "Added MXOld" << ExtendedShellMinusXDirOld.size() << " position query on -X to the extended interaction shell. " << std::endl;

	std::cout << "Added PXNew" << ExtendedShellPlusXDirNew.size() << " position query on +X to the extended interaction shell. " << std::endl;
	std::cout << "Added PXOld" << ExtendedShellPlusXDirOld.size() << " position query on +X to the extended interaction shell. " << std::endl;

	std::cout << "Added MYNew" << ExtendedShellMinusYDirNew.size() << " position query on -Y to the extended interaction shell. " << std::endl;
	std::cout << "Added MYOld" << ExtendedShellMinusYDirOld.size() << " position query on -Y to the extended interaction shell. " << std::endl;

    std::cout << "Added PYNew" << ExtendedShellPlusYDirNew.size() << " position query on +Y to the extended interaction shell. " << std::endl;
    std::cout << "Added PYOld" << ExtendedShellPlusYDirOld.size() << " position query on +Y to the extended interaction shell. " << std::endl;

	std::cout << "Added MZNew" << ExtendedShellMinusZDirNew.size() << " position query on -Z to the extended interaction shell. " << std::endl;
	std::cout << "Added MZOld" << ExtendedShellMinusZDirOld.size() << " position query on -Z to the extended interaction shell. " << std::endl;

	std::cout << "Added PZNew" << ExtendedShellPlusZDirNew.size() << " position query on +Z to the extended interaction shell. " << std::endl;
	std::cout << "Added PZOld" << ExtendedShellPlusZDirOld.size() << " position query on +Z to the extended interaction shell. " << std::endl;

	std::cout << "Added " << ExtendedShell.size() << " position to the extended interaction shell. " << std::endl;



	Energy = getInternalEnergyCurrentConfiguration(ingredients);

	std::cout << "Synchronize FeatureWangLandauExtendedShellInteraction -> E = " << Energy << std::endl;
}

/**
 * @details If not compiled with DEBUG flag this function only returns the content
 * of the lookup table probabilityLookup. If compiled with DEBUG flag it checks 
 * that the attribute tags typeA, typeB are within the allowed range.
 * @param typeA monomer attribute type in range [1,255]
 * @param typeB monomer attribute type in range [1,255]
 * @throw std::runtime_error In debug mode, if types are not in range [1,255]
 **/
template<template<typename> class LatticeClassType>
inline double FeatureExtendedShellInteraction<LatticeClassType>::getProbabilityFactor(int32_t typeA,
									     int32_t typeB) const
{
#ifdef DEBUG
  //extra checks only in debug mode, because this is very frequently called
  //and this costs performance
  if(typeA<0 || typeA>255 || typeB<0 || typeB>255){
    std::stringstream errormessage;
    errormessage<<"***FeatureNaNInteractionSc::getInteraction(typeA,typeB)***\n";
    errormessage<<"probability undefined between types "<<typeA<<" and "<<typeB<<std::endl;
    errormessage<<"types are out of the allowed range";
    throw std::runtime_error(errormessage.str());
  }
#endif /*DEBUG*/

  return probabilityLookup[typeA][typeB];

}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * - !contactInteraction
 * .
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileReader File importer for the bfm-file
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::exportRead(FileImport< IngredientsType >& fileReader)
{
  typedef FeatureExtendedShellInteraction<LatticeClassType> my_type;
  fileReader.registerRead("!nn_interaction",new ReadExtendedShellInteraction<my_type>(*this));
  fileReader.registerRead("!nn_interaction_shell_radius",new ReadExtendedShellInteractionRadius<my_type>(*this));

}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * - !contact_interaction
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileWriter File writer for the bfm-file.
 */
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  typedef FeatureExtendedShellInteraction<LatticeClassType> my_type;
  fileWriter.registerWrite("!nn_interaction",new WriteExtendedShellInteraction<my_type>(*this));
  fileWriter.registerWrite("!nn_interaction_shell_radius",new WriteExtendedShellInteractionRadius<my_type>(*this));
  fileWriter.registerWrite("#!energy",new WriteExtendedShellInteractionEnergy<IngredientsType>(fileWriter.getIngredients_()));

}

/**
 * @details occupies the lattice with the attribute tags of the monomers
 * as this is required to determine the contact interactions in this feature.
 * An additional check is performed asserting that the tags are in the range [1,255]
 *
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @throw std::runtime_error In case a monomer has attribute tag not in [1,255]
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureExtendedShellInteraction<LatticeClassType>::fillLattice(IngredientsType& ingredients)
{
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

    for(size_t n=0;n<molecules.size();n++)
    {
    	// iterates over the monomer cube.

    	// check for multiple occupation
    	/*for ( int x = -1 ; x < 2; ++x)
    		for ( int y = -1 ; y < 2; ++y)
    			for ( int z = -1 ; z < 2; ++z)
    			{
    				VectorInt3 pos = molecules[n]+VectorInt3(x,y,z);

    				if (ingredients.getLatticeEntry(pos))
    				{
    					std::ostringstream errorMessage;
    					errorMessage << "FeatureExcludedVolume<>::fillLattice(): lattice already occupied when trying to write monomer ";
    					errorMessage << n << " at " << molecules[n] << ", cube corner at " << pos << ".\n";
    					throw std::runtime_error(errorMessage.str());
    				}
    				//here we simply set a one on every occupied lattice
    				//site. this assumes that 1 can be cast to LatticeType,
    				//even though in principle LatticeType could be anything.
    				//note that this may be just a preliminiary initialization,
    				//as other features may assign more specific values to the
    				//lattice site.
    			}*/

    	lattice_value_type attribute=lattice_value_type(molecules[n].getAttributeTag());

    	if(int32_t(attribute)!=molecules[n].getAttributeTag()){
    		std::stringstream errormessage;
    		errormessage<<"***FeatureExtendedShellInteraction::fillLattice()***\n";
    		errormessage<<"type "<<attribute<<" is out of the allowed range";

    		throw std::runtime_error(errormessage.str());
    	}

    	ingredients.setLatticeEntry(molecules[n],attribute);


    }

}


/**
 * @details The function calculates the factor for the acceptance probability
 * for the local move given as argument. The calculation is based on the lattice
 * entries in the vicinity of the monomer to be moved. If the move is accepted,
 * 12 new contacts can potentially be made, and 12 contacts are lost. Thus a number
 * of 24 lattice positions around the monomer have to be checked.
 *
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move reference to the local move for which the calculation is performed
 * @return acceptance probability factor for the move arising from nearest neighbor contacts
 **/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureExtendedShellInteraction<LatticeClassType>::calculateAcceptanceProbability(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{
    
    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();



    //double prob=1.0;
    int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

    double E_old = 0.0;
    double E_new = 0.0;

/*
    VectorInt3 newPos=oldPos+direction;
    for(size_t shell=0;shell<ExtendedShell.size();shell++)
    {
    	int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShell[shell]));
    	int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(newPos+ExtendedShell[shell]));

    	E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
    	E_new += getExtendedShellInteraction(monoType,latticeEntryNew);

    }


    std::cout << "Diff " << (E_new-E_old) ;
    E_old = 0.0;
    E_new = 0.0;
*/

    switch( (direction.getX() & 3) + ((direction.getY() &3) << 2) + ((direction.getZ() &3) << 4) ){

    case 1: //std::cout << " +X ";

				for(size_t shell=0;shell<ExtendedShellPlusXDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusXDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellPlusXDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusXDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 3: //std::cout << " -X ";

				for(size_t shell=0;shell<ExtendedShellMinusXDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusXDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusXDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusXDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 4: //std::cout << " +Y ";

				for(size_t shell=0;shell<ExtendedShellPlusYDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusYDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellPlusYDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusYDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 12: //std::cout << " -Y ";

				for(size_t shell=0;shell<ExtendedShellMinusYDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusYDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusYDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusYDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 16: //std::cout << " +Z ";
    			for(size_t shell=0;shell<ExtendedShellPlusZDirNew.size();shell++)
    	    	{
    	    		int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusZDirNew[shell]));
    	    		E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
    	    	}
    	    	for(size_t shell=0;shell<ExtendedShellPlusZDirOld.size();shell++)
    	    	{
    	    		int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusZDirOld[shell]));
    	    		E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
    	    	}

    	    //	std::cout <<  " -> " << (E_new-E_old) << std::endl;
    		break;

    case 48: //std::cout << " -Z ";

				for(size_t shell=0;shell<ExtendedShellMinusZDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusZDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusZDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusZDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    default: std::stringstream errormessage;
    			errormessage<<"FeatureExtendedShellInteraction::calculateAcceptanceProbability.\n";
    			errormessage<<"Mismatch in direction calculation\n";
    			throw std::runtime_error(errormessage.str());
    		break;

    }

    double prob = std::exp(-(E_new-E_old));
    return prob;


}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureExtendedShellInteraction<LatticeClassType>::calculateInteractionDifference(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{
    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();



    //double prob=1.0;
    int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

    double E_old = 0.0;
    double E_new = 0.0;

/*
    VectorInt3 newPos=oldPos+direction;
    for(size_t shell=0;shell<ExtendedShell.size();shell++)
    {
    	int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShell[shell]));
    	int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(newPos+ExtendedShell[shell]));

    	E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
    	E_new += getExtendedShellInteraction(monoType,latticeEntryNew);

    }


    std::cout << "Diff " << (E_new-E_old) ;
    E_old = 0.0;
    E_new = 0.0;
*/

    switch( (direction.getX() & 3) + ((direction.getY() &3) << 2) + ((direction.getZ() &3) << 4) ){

    case 1: //std::cout << " +X ";

				for(size_t shell=0;shell<ExtendedShellPlusXDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusXDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellPlusXDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusXDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 3: //std::cout << " -X ";

				for(size_t shell=0;shell<ExtendedShellMinusXDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusXDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusXDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusXDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 4: //std::cout << " +Y ";

				for(size_t shell=0;shell<ExtendedShellPlusYDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusYDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellPlusYDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusYDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 12: //std::cout << " -Y ";

				for(size_t shell=0;shell<ExtendedShellMinusYDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusYDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusYDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusYDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    case 16: //std::cout << " +Z ";
    			for(size_t shell=0;shell<ExtendedShellPlusZDirNew.size();shell++)
    	    	{
    	    		int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusZDirNew[shell]));
    	    		E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
    	    	}
    	    	for(size_t shell=0;shell<ExtendedShellPlusZDirOld.size();shell++)
    	    	{
    	    		int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellPlusZDirOld[shell]));
    	    		E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
    	    	}

    	    //	std::cout <<  " -> " << (E_new-E_old) << std::endl;
    		break;

    case 48: //std::cout << " -Z ";

				for(size_t shell=0;shell<ExtendedShellMinusZDirNew.size();shell++)
				{
					int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusZDirNew[shell]));
					E_new += getExtendedShellInteraction(monoType,latticeEntryNew);
				}
				for(size_t shell=0;shell<ExtendedShellMinusZDirOld.size();shell++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(oldPos+ExtendedShellMinusZDirOld[shell]));
					E_old += getExtendedShellInteraction(monoType,latticeEntryOld);
				}

			//	std::cout <<  " -> " << (E_new-E_old) << std::endl;
			break;

    default: std::stringstream errormessage;
    			errormessage<<"FeatureExtendedShellInteraction::calculateAcceptanceProbability.\n";
    			errormessage<<"Mismatch in direction calculation\n";
    			throw std::runtime_error(errormessage.str());
    		break;

    }

   // double prob = std::exp(-(E_new-E_old));
   // return prob;

    return (E_new-E_old);
}

/**
 * @param typeA monomer attribute tag in range [1,255]
 * @param typeB monomer attribute tag in range [1,255]
 * @param interaction energy between typeA and typeB
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 **/
template<template<typename> class LatticeClassType>
void FeatureExtendedShellInteraction<LatticeClassType>::setExtendedShellInteraction(int32_t typeA,
								     int32_t typeB,
								     double energy)
{
    if(0<typeA && typeA<=255 && 0<typeB && typeB<=255)
      {
        interactionTable[typeA][typeB]=energy;
        interactionTable[typeB][typeA]=energy;
        probabilityLookup[typeA][typeB]=exp(-energy);
        probabilityLookup[typeB][typeA]=exp(-energy);
        std::cout<<"set interation between types ";
	std::cout<<typeA<<" and "<<typeB<<" to "<<energy<<"kT\n";
      }
    else
      {
	std::stringstream errormessage;
	errormessage<<"FeatureExtendedShellInteraction::setNNInteraction(typeA,typeB,energy).\n";
	errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
	throw std::runtime_error(errormessage.str());
      }
}

/**
 * @param typeA monomer attribute tag in range [1,255]
 * @param typeB monomer attribute tag in range [1,255]
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 * @return interaction energy per nearest neighbor contact for typeA,typeB
 **/
template<template<typename> class LatticeClassType>
double FeatureExtendedShellInteraction<LatticeClassType>::getExtendedShellInteraction(int32_t typeA,
								       int32_t typeB) const
{

#ifdef DEBUG
  //extra checks only in debug mode, because this is very frequently called
  //and this costs performance
  if(typeA<0 || typeA>255 || typeB<0 || typeB>255){
    std::stringstream errormessage;
    errormessage<<"***FeatureExtendedShellInteraction::getExtendedShellInteraction((typeA,typeB)***\n";
    errormessage<<"Energy undefined between types "<<typeA<<" and "<<typeB<<std::endl;
    errormessage<<"types are out of the allowed range";
    throw std::runtime_error(errormessage.str());
  }
#endif /*DEBUG*/

  return interactionTable[typeA][typeB];
}


//calculate the center of mass of the affected monomers in z-direction
template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureExtendedShellInteraction<LatticeClassType>::getInternalEnergyCurrentConfiguration(const IngredientsType& ingredients) const
{

	//prepare the contact shell, in which the contacts are counted
	// contact shell is calulated in synchronize() -> ExtendedShell


	VectorInt3 pos;
	double EnergyActual=0.0;

	//loop through all polymer
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();
		pos=ingredients.getMolecules()[n];


		for(size_t contactNo=0;contactNo<ExtendedShell.size();contactNo++)
		{
			int32_t latticeEntry=int32_t(ingredients.getLatticeEntry(pos+ExtendedShell[contactNo]));

			EnergyActual += getExtendedShellInteraction(monoType,latticeEntry);
		}
	}

	// factor 0.5 due to the double sum in hamiltonian
	return 0.5*EnergyActual;
}


/**
 * @param radius maximum radius of interaction shell
 * @throw std::runtime_error In case radius out of range 0 < r <= 6.0
 **/
template<template<typename> class LatticeClassType>
void FeatureExtendedShellInteraction<LatticeClassType>::setExtendedShellInteractionRadius(double radius)
{
    if(2.0<=radius && radius<=6.0)
      {
    	ExtendedShellRadius = radius;
        std::cout<<"set interation radius to "<<ExtendedShellRadius<<"\n";
      }
    else
      {
	std::stringstream errormessage;
	errormessage<<"FeatureExtendedShellInteractionRadius::setExtendedShellInteractionRadius(radius).\n";
	errormessage<<"Interaction Radius out of range 2 <= r <= 6.0\n";
	throw std::runtime_error(errormessage.str());
      }
}

/**
 * @param typeA monomer attribute tag in range [1,255]
 * @param typeB monomer attribute tag in range [1,255]
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 * @return interaction energy per nearest neighbor contact for typeA,typeB
 **/
template<template<typename> class LatticeClassType>
double FeatureExtendedShellInteraction<LatticeClassType>::getExtendedShellInteractionRadius() const
{
  return ExtendedShellRadius;
}


#endif /*FEATURE_EXTENDED_SHELL_INTERACTION_H*/

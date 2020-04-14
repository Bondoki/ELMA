#ifndef FEATURE_WANG_LANDAU_NEXT_NEIGHBOR__H
#define FEATURE_WANG_LANDAU_NEXT_NEIGHBOR__H

#include <vector>
#include <iostream>
#include <algorithm>    // for min_element

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/utility/DistanceCalculation.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureNNInteractionReadWrite.h> //read/write
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>

#include "Histogram1D.h"
#include "HistogramGeneralStatistik1D.h"


template<template<typename> class FeatureLatticeType>
class FeatureWangLandauNextNeighbor:public Feature
{
public:
	typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

	//! Type for the underlying lattice, used as template parameter for FeatureLatticeType<...>
	typedef uint8_t lattice_value_type;

	typedef LOKI_TYPELIST_2(
	      FeatureAttributes,
	      FeatureExcludedVolumeSc<FeatureLatticeType<lattice_value_type> >)
	    required_features_front;
		
	FeatureWangLandauNextNeighbor()
	:allMonomersAffected(true)
	,Energy(0.0)
	,diffEnergy(0.0)
	,modificationFactor(1.01)//std::exp(1))
	,windowingState(false)
	{
		//initialize the energy and probability lookups with default values
		  for(size_t n=0;n<256;n++)
		    {
		      for(size_t m=0;m<256;m++)
		        {
			  interactionTable[m][n]=0.0;
			  probabilityLookup[m][n]=1.0;
		        }
		    }

		  lnDOSmin = 1.01;
		  numHistoVisits = 1.0;
	}
	
	virtual ~FeatureWangLandauNextNeighbor(){}
    
        //set the monomers affected
	void setAffectedMonomers(std::set<size_t> m){allMonomersAffected=false; affectedMonomers=m;}
	void setAffectedMonomersAll(bool flag){allMonomersAffected=flag;}
	
	//for all unknown moves: does nothing
	template < class IngredientsType> 
	bool checkMove( const IngredientsType& ingredients, const MoveBase& move ) const;

	//overload for MoveLocalSc
	template < class IngredientsType> 
	bool checkMove( const IngredientsType& ingredients, MoveLocalSc& move );

	//! check move for bcc-BFM local move. always throws std::runtime_error
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,const MoveLocalBcc& move) const;

	//for unknown moves(does nothing)
	template<class IngredientsType> 
	void applyMove(IngredientsType& ing, const MoveBase& move){}
	
	//for moves of type MoveLocalSc
	template<class IngredientsType> 
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);

	//! apply function for bcc-BFM local move (always throws std::runtime_error)
	template<class IngredientsType>
	void applyMove(const IngredientsType& ing, const MoveLocalBcc& move);


	
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	
	
	

	template<class IngredientsType>
	void rejectMove(IngredientsType& ingredients);



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


	template<class IngredientsType>
	double getInternalEnergy(const IngredientsType& ingredients,size_t index=0,VectorInt3 direction=VectorInt3(0,0,0)) const;


	template<class IngredientsType>
	double getInternalEnergyDifference(const IngredientsType& ingredients,size_t index,VectorInt3 direction) const;


	template<class IngredientsType>
	double getInternalEnergyCurrentConfiguration(const IngredientsType& ingredients) const;

	double getEnergy() const
		{
			return Energy;
		}

	//!adds interaction energy between two types of monomers
	  void setNNInteraction(int32_t typeA,int32_t typeB,double energy);

	  //!returns the interaction energy between two types of monomers
	  double getNNInteraction(int32_t typeA,int32_t typeB) const;

	bool isEnergyInWindow() const {
		return windowingState;
	}

	void setWindowState(bool isEnergyInWindow, double minWin, double maxWin) {
		this->windowingState = isEnergyInWindow;
		this->minWin = minWin;
		this->maxWin = maxWin;
	}

	  //!export bfm-file read command !nn_interaction
	    template <class IngredientsType>
	    void exportRead(FileImport <IngredientsType>& fileReader);

	    //!export bfm-file write command !nn_interaction
	    template <class IngredientsType>
	    void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

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

private:
	
	
	//contains the indices of the monomers with type affectedMonomerType
	bool allMonomersAffected;
	std::set<size_t> affectedMonomers;
	
	Histogram1D HG_VisitsEnergyStates;
	HistogramGeneralStatistik1D HG_LnDOS;
	double Energy;
	double diffEnergy;
	double modificationFactor;

	Histogram1D HG_TotalVisitsEnergyStates;
	
	//Histogram1D HG_VisitsEnergyStates;


  //! Interaction energies between monomer types. Max. type=255 given by max(uint8_t)=255
  double interactionTable[256][256];

  //! Lookup table for exp(-interactionTable[a][b])
  double probabilityLookup[256][256];

  //! isEnergyInWindow
  bool windowingState;
  double minWin;
  double maxWin;

  //! Returns this feature's factor for the acceptance probability for the given Monte Carlo move
  template<class IngredientsType>
  double calculateAcceptanceProbability(const IngredientsType& ingredients,
					const MoveLocalSc& move) const;


  //! Returns this feature's factor for the energy difference for the given Monte Carlo move
   template<class IngredientsType>
   double calculateInteractionDifference(const IngredientsType& ingredients,
			const MoveLocalSc& move) const;

  //! Occupies the lattice with the attribute tags of all monomers
  template<class IngredientsType>
  void fillLattice(IngredientsType& ingredients);

  //! Access to array probabilityLookup with extra checks in Debug mode
  double getProbabilityFactor(int32_t typeA,int32_t typeB) const;


  double lnDOSmin;
  double numHistoVisits;

};



///////////////////////////////////////////////////////////////////////////////
////////////////////Implementation of methods////////////////////////////

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
bool FeatureWangLandauNextNeighbor<LatticeClassType>::checkMove(const IngredientsType& ingredients,
							 const MoveLocalBcc& move) const
{
  //throw exception in case someone accidentaly uses a bcc-BFM move with this feature
  std::stringstream errormessage;
  errormessage<<"FeatureWangLandauNextNeighbor::checkMove(...):\n";
  errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
  throw std::runtime_error(errormessage.str());

  return false;
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
void FeatureWangLandauNextNeighbor<LatticeClassType>::applyMove(const IngredientsType& ing,
							 const MoveLocalBcc& move)
{
  //throw exception in case someone accidentaly uses a bcc-BFM move with this feature
  std::stringstream errormessage;
  errormessage<<"FeatureNNInteractionSc::applyMove(...):\n";
  errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
  throw std::runtime_error(errormessage.str());

}


template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureWangLandauNextNeighbor<LatticeClassType>::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;//nothing to do for unknown moves
}


template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureWangLandauNextNeighbor<LatticeClassType>::checkMove(const IngredientsType& ingredients, MoveLocalSc& move)
{
	
	/*if(!allMonomersAffected){
		if(affectedMonomers.count(move.getIndex())==0) return true;
	}
	*/
	
	// uncomment for testing and comparing different implementations
	//double EOld = getInternalEnergyCurrentConfiguration(ingredients);
	//double ENew = getInternalEnergy(ingredients,move.getIndex(),move.getDir());
	//double dE = getInternalEnergyDifference(ingredients,move.getIndex(),move.getDir());
	//std::cout << "EOld: " << EOld <<"\t ENew: " << (ENew)  << " \t (ENew-EOld): " << (ENew-EOld)  << " \t dE: " << (dE) << std::endl;

	if(windowingState == false)
	{
	diffEnergy = calculateInteractionDifference(ingredients,move);
	//std::cout << "EnergyOld: " << Energy <<"\t EnergyNew: " << (Energy+diffEnergy)  << " \t dEnergy: " << (diffEnergy) << std::endl;

	//double prob=calculateAcceptanceProbability(ingredients,move);
	//std::cout << "from probability dEnergy: " << (-std::log(prob)) << std::endl;

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
		if( ( ((Energy+diffEnergy) >= minWin) && ((Energy+diffEnergy) <= maxWin)) )
		{
			//std::cout << "EnergyOld: " << Energy <<"\t EnergyNew: " << (Energy+diffEnergy)  << " \t dEnergy: " << (diffEnergy) << std::endl;

			//double prob=calculateAcceptanceProbability(ingredients,move);
			//std::cout << "from probability dEnergy: " << (-std::log(prob)) << std::endl;

			double p=1.0;

			// p = min(1,g(E_old)/g(E_new))

			if(HG_LnDOS.getCountAt(Energy) < HG_LnDOS.getCountAt(Energy+diffEnergy))
				p=std::exp(HG_LnDOS.getCountAt(Energy)-HG_LnDOS.getCountAt(Energy+diffEnergy));

			/*double m = (HG_LnDOS.getCountAt(maxWin)-HG_LnDOS.getCountAt(minWin))/(maxWin-minWin);
			double n = HG_LnDOS.getCountAt(minWin)-m*minWin;
			double bOld = m*Energy+n;
			double bNew = m*(Energy+diffEnergy)+n;

			//std::cout << "m: " << m <<"\t n: " << (n)  << " \t b: " << (b) << std::endl;
			if((bOld != 0.0) && (bNew != 0))
				p=p*(std::exp(bOld)/std::exp(bNew));
			*/

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

template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureWangLandauNextNeighbor<LatticeClassType>::applyMove(IngredientsType& ingredients, const MoveLocalSc& move)
{
	//update the book keeping variables in case the move was accepted
	//EnergyOld=EnergyNew;
	Energy += diffEnergy;

	//HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
	//HG_VisitsEnergyStates.addValue(Energy, 1.0);
	//HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);

	/*if(windowingState == false)
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

	}*/

	/*if(	HG_LnDOS.getNumCountAt(Energy) != 0)
	{
		if(windowingState == false)
		{
			HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
			HG_VisitsEnergyStates.addValue(Energy, 1.0);
			HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);
		}
		else
		{
			if(  (Energy > minWin-4.0*HG_LnDOS.getBinwidth()) && (Energy < maxWin+4.0*HG_LnDOS.getBinwidth()) )
			{
				HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
				HG_VisitsEnergyStates.addValue(Energy, 1.0);
				HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);
			}
		}
		//if(lnDOSmin < HG_LnDOS.getCountAt(Energy))
		//	lnDOSmin = HG_LnDOS.getCountAt(Energy);
	}
	else
	{
		//find minimum

		double eln = HG_LnDOS.getFirstMomentInBin(0);
		double visits = 1.0;
		for (size_t n=1; n < HG_LnDOS.getNBins(); n++) {
			if ((HG_LnDOS.getFirstMomentInBin(n) < eln && HG_LnDOS.getFirstMomentInBin(n) != 0 ) || eln == 0)
			{
				eln = HG_LnDOS.getFirstMomentInBin(n);
				visits = HG_VisitsEnergyStates.getCountInBin(n);
			}
		}

		if(visits != 0)
		{
			numHistoVisits = visits;
			lnDOSmin = eln;
		}

		//std::cout << "lnDOSmin " << lnDOSmin  << " visits: " << numHistoVisits <<  std::endl;
		//std::cout << "e " << e   << " visits: " << visits <<  std::endl;




		//if(lnDOSmin < HG_LnDOS.getCountAt(Energy))
		//				lnDOSmin = HG_LnDOS.getCountAt(Energy);
		if(windowingState == false)
		{
			HG_LnDOS.resetValue(Energy, lnDOSmin );
			HG_VisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
			HG_TotalVisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
		}
		else
		{
			if(  (Energy > minWin-4.0*HG_LnDOS.getBinwidth()) && (Energy < maxWin+4.0*HG_LnDOS.getBinwidth()) )
			{
				HG_LnDOS.resetValue(Energy, lnDOSmin );
				HG_VisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
				HG_TotalVisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
			}
		}
	}
	*/

}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureWangLandauNextNeighbor<LatticeClassType>::rejectMove(IngredientsType& ingredients)//, const MoveLocalSc& move)
{
	//HG_LnDOS.resetValue(Energy, HG_LnDOS.getCountAt(Energy)+std::log(modificationFactor));
	//HG_VisitsEnergyStates.addValue(Energy, 1.0);
	//HG_TotalVisitsEnergyStates.addValue(Energy, 1.0);

	/*
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
	*/

	/*	//if(lnDOSmin < HG_LnDOS.getCountAt(Energy))
		//	lnDOSmin = HG_LnDOS.getCountAt(Energy);
	}
	else
	{
		//find minimum

		double eln = HG_LnDOS.getFirstMomentInBin(0);
		double visits = 1.0;
		for (size_t n=1; n < HG_LnDOS.getNBins(); n++) {
			if ((HG_LnDOS.getFirstMomentInBin(n) < eln && HG_LnDOS.getFirstMomentInBin(n) != 0 ) || eln == 0)
			{
				eln = HG_LnDOS.getFirstMomentInBin(n);
				visits = HG_VisitsEnergyStates.getCountInBin(n);
			}
		}

		if(visits != 0)
		{
			numHistoVisits = visits;
			lnDOSmin = eln;
		}

		//std::cout << "lnDOSmin " << lnDOSmin  << " visits: " << numHistoVisits <<  std::endl;
		//std::cout << "e " << e   << " visits: " << visits <<  std::endl;




		//if(lnDOSmin < HG_LnDOS.getCountAt(Energy))
		//				lnDOSmin = HG_LnDOS.getCountAt(Energy);
		if(windowingState == false)
		{
			HG_LnDOS.resetValue(Energy, lnDOSmin );
			HG_VisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
			HG_TotalVisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
		}
		else
		{
			if(  (Energy > minWin-4.0*HG_LnDOS.getBinwidth()) && (Energy < maxWin+4.0*HG_LnDOS.getBinwidth()) )
			{
				HG_LnDOS.resetValue(Energy, lnDOSmin );
				HG_VisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
				HG_TotalVisitsEnergyStates.addValue(Energy, numHistoVisits);//lnDOSmin/std::log(modificationFactor));//1.0);
			}
		}
	}
	*/
}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureWangLandauNextNeighbor<LatticeClassType>::synchronize(IngredientsType& ingredients)
{
	//refill the lattice with attribute tags
	//caution: this overwrites, what is currently written on the lattice
	fillLattice(ingredients);

	Energy = getInternalEnergyCurrentConfiguration(ingredients);
	

	std::cout << "Synchronize FeatureWangLandauNextNeighbor -> E = " << Energy << std::endl;
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
void FeatureWangLandauNextNeighbor<LatticeClassType>::fillLattice(IngredientsType& ingredients)
{
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

    for(size_t n=0;n<molecules.size();n++)
    {
        VectorInt3 pos=molecules[n];
	lattice_value_type attribute=lattice_value_type(molecules[n].getAttributeTag());

	if(int32_t(attribute)!=molecules[n].getAttributeTag()){
	  std::stringstream errormessage;
	  errormessage<<"***FeatureNNInteractionSc::fillLattice()***\n";
	  errormessage<<"type "<<attribute<<" is out of the allowed range";

	  throw std::runtime_error(errormessage.str());
	}

        ingredients.setLatticeEntry(pos,attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,0,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,1),attribute);

    }

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
inline double FeatureWangLandauNextNeighbor<LatticeClassType>::getProbabilityFactor(int32_t typeA,
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
double FeatureWangLandauNextNeighbor<LatticeClassType>::calculateAcceptanceProbability(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{

    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();

    double prob=1.0;
    int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

    /*get two directions perpendicular to vector directon of the move*/
    VectorInt3 perp1,perp2;
    /* first perpendicular direction is either (0 1 0) or (1 0 0)*/
    int32_t x1=((direction.getX()==0) ? 1 : 0);
    int32_t y1=((direction.getX()!=0) ? 1 : 0);
    perp1.setX(x1);
    perp1.setY(y1);
    perp1.setZ(0);

    /* second perpendicular direction is either (0 0 1) or (0 1 0)*/
    int32_t y2=((direction.getZ()==0) ? 0 : 1);
    int32_t z2=((direction.getZ()!=0) ? 0 : 1);
    perp2.setX(0);
    perp2.setY(y2);
    perp2.setZ(z2);

    //the probability is calculated by going through all possible lattice sites
    //at which the contacts may have changed. At every site the type of the
    //monomer sitting there is retrieved from the lattice. the additional
    //factor for the probability (exp(-deltaE/kT)) is retrieved from the
    //lookup using getProbabilityFactor. For new contacts this factor is multiplied
    //with the probability, for contacts taken away the probability is devided.
    VectorInt3 actual=oldPos;

    //first check front,i.e newly acquired contacts
    if(direction.getX()>0 || direction.getY()>0 || direction.getZ()>0) actual+=direction;
    actual+=direction;

    actual-=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+direction;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    //now check back side (contacts taken away)
    double prob_div=1.0;
    actual=oldPos;
    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) actual-=direction;
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2-direction;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    prob/=prob_div;
    return prob;

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
double FeatureWangLandauNextNeighbor<LatticeClassType>::calculateInteractionDifference(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{

    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();

    double deltaE=0.0;
    int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

    /*get two directions perpendicular to vector directon of the move*/
    VectorInt3 perp1,perp2;
    /* first perpendicular direction is either (0 1 0) or (1 0 0)*/
    int32_t x1=((direction.getX()==0) ? 1 : 0);
    int32_t y1=((direction.getX()!=0) ? 1 : 0);
    perp1.setX(x1);
    perp1.setY(y1);
    perp1.setZ(0);

    /* second perpendicular direction is either (0 0 1) or (0 1 0)*/
    int32_t y2=((direction.getZ()==0) ? 0 : 1);
    int32_t z2=((direction.getZ()!=0) ? 0 : 1);
    perp2.setX(0);
    perp2.setY(y2);
    perp2.setZ(z2);

    //the probability is calculated by going through all possible lattice sites
    //at which the contacts may have changed. At every site the type of the
    //monomer sitting there is retrieved from the lattice. the additional
    //factor for the probability (exp(-deltaE/kT)) is retrieved from the
    //lookup using getProbabilityFactor. For new contacts this factor is multiplied
    //with the probability, for contacts taken away the probability is devided.
    VectorInt3 actual=oldPos;

    //first check front,i.e newly acquired contacts
    if(direction.getX()>0 || direction.getY()>0 || direction.getZ()>0) actual+=direction;
    actual+=direction;

    actual-=perp1;
    deltaE += getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+direction;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    deltaE +=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    //now check back side (contacts taken away)
   // double prob_div=1.0;
    actual=oldPos;
    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) actual-=direction;
    actual-=perp1;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2-direction;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    deltaE -=getNNInteraction(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    //prob/=prob_div;
    //return prob;
    return deltaE;

}


/**
 * @param typeA monomer attribute tag in range [1,255]
 * @param typeB monomer attribute tag in range [1,255]
 * @param interaction energy between typeA and typeB
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 **/
template<template<typename> class LatticeClassType>
void FeatureWangLandauNextNeighbor<LatticeClassType>::setNNInteraction(int32_t typeA,
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
	errormessage<<"FeatureWangLandauNextNeighborFeatureNNInteractionSc::setNNInteraction(typeA,typeB,energy).\n";
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
double FeatureWangLandauNextNeighbor<LatticeClassType>::getNNInteraction(int32_t typeA,
								       int32_t typeB) const
{

#ifdef DEBUG
  //extra checks only in debug mode, because this is very frequently called
  //and this costs performance
  if(typeA<0 || typeA>255 || typeB<0 || typeB>255){
    std::stringstream errormessage;
    errormessage<<"***FeatureWangLandauNextNeighbor::setNNInteraction(typeA,typeB)***\n";
    errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
    errormessage<<"types are out of the allowed range";
    throw std::runtime_error(errormessage.str());
  }
#endif /*DEBUG*/

  return interactionTable[typeA][typeB];

}

//calculate the center of mass of the affected monomers in z-direction
template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureWangLandauNextNeighbor<LatticeClassType>::getInternalEnergy(const IngredientsType& ingredients,size_t index, VectorInt3 direction) const
{
	std::vector<VectorInt3> contactSites;

		//prepare the contact shell, in which the contacts are counted

		contactSites.push_back(VectorInt3(2,0,0));
		contactSites.push_back(VectorInt3(2,1,0));
		contactSites.push_back(VectorInt3(2,0,1));
		contactSites.push_back(VectorInt3(2,1,1));
		contactSites.push_back(VectorInt3(0,2,0));
		contactSites.push_back(VectorInt3(1,2,0));
		contactSites.push_back(VectorInt3(0,2,1));
		contactSites.push_back(VectorInt3(1,2,1));
		contactSites.push_back(VectorInt3(0,0,2));
		contactSites.push_back(VectorInt3(1,0,2));
		contactSites.push_back(VectorInt3(0,1,2));
		contactSites.push_back(VectorInt3(1,1,2));

		contactSites.push_back(VectorInt3(-1,0,0));
		contactSites.push_back(VectorInt3(-1,1,0));
		contactSites.push_back(VectorInt3(-1,0,1));
		contactSites.push_back(VectorInt3(-1,1,1));
		contactSites.push_back(VectorInt3(0,-1,0));
		contactSites.push_back(VectorInt3(1,-1,0));
		contactSites.push_back(VectorInt3(0,-1,1));
		contactSites.push_back(VectorInt3(1,-1,1));
		contactSites.push_back(VectorInt3(0,0,-1));
		contactSites.push_back(VectorInt3(1,0,-1));
		contactSites.push_back(VectorInt3(0,1,-1));
		contactSites.push_back(VectorInt3(1,1,-1));

		std::vector<VectorInt3> vicinitySites;
		vicinitySites.push_back(VectorInt3(0,0,0));
		vicinitySites.push_back(VectorInt3(1,0,0));
		vicinitySites.push_back(VectorInt3(0,1,0));
		vicinitySites.push_back(VectorInt3(1,1,0));
		vicinitySites.push_back(VectorInt3(0,0,1));
		vicinitySites.push_back(VectorInt3(1,0,1));
		vicinitySites.push_back(VectorInt3(0,1,1));
		vicinitySites.push_back(VectorInt3(1,1,1));



		VectorInt3 pos;
		double Energy=0.0;

		// new contacts
		VectorInt3 posIndex=ingredients.getMolecules()[index]+direction;
		VectorInt3 posIndexOld=ingredients.getMolecules()[index];

		//loop through all polymer
		for(size_t n=0;n<ingredients.getMolecules().size();n++)
		{
			if(n != index)
			{
				int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();
				pos=ingredients.getMolecules()[n];

				for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
				{
					//exclude the moveing monomer
					bool possibleCount = true;
					for(size_t vicinityNo=0; vicinityNo<vicinitySites.size();vicinityNo++)
					{
						//compare with old position due to not-uptdated lattice in checkmove
						if((pos+contactSites[contactNo]) == (posIndexOld+vicinitySites[vicinityNo]))
							possibleCount = false;
					}

					if(possibleCount == true)
					{
						int32_t latticeEntry=int32_t(ingredients.getLatticeEntry(pos+contactSites[contactNo]));
						Energy += 0.5*getNNInteraction(monoType,latticeEntry);
					}
				}
			}
		}

		// counting the moving monomer


		for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
		{
			//exclude the old position
			bool possibleCount = true;
			for(size_t vicinityNo=0; vicinityNo<vicinitySites.size();vicinityNo++)
			{
				if((posIndex+contactSites[contactNo]) == (posIndexOld+vicinitySites[vicinityNo]))
					possibleCount = false;
			}

			if(possibleCount == true)
			{
				int32_t monoTypeIndex=ingredients.getMolecules()[index].getAttributeTag();
				int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(posIndex+contactSites[contactNo]));
				Energy += getNNInteraction(monoTypeIndex,latticeEntryNew);
			}
		}


		// factor 0.5 due to the double sum in hamiltonian
		return Energy;
}

//calculate the center of mass of the affected monomers in z-direction
template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureWangLandauNextNeighbor<LatticeClassType>::getInternalEnergyCurrentConfiguration(const IngredientsType& ingredients) const
{
	std::vector<VectorInt3> contactSites;

	//prepare the contact shell, in which the contacts are counted

	contactSites.push_back(VectorInt3(2,0,0));
	contactSites.push_back(VectorInt3(2,1,0));
	contactSites.push_back(VectorInt3(2,0,1));
	contactSites.push_back(VectorInt3(2,1,1));
	contactSites.push_back(VectorInt3(0,2,0));
	contactSites.push_back(VectorInt3(1,2,0));
	contactSites.push_back(VectorInt3(0,2,1));
	contactSites.push_back(VectorInt3(1,2,1));
	contactSites.push_back(VectorInt3(0,0,2));
	contactSites.push_back(VectorInt3(1,0,2));
	contactSites.push_back(VectorInt3(0,1,2));
	contactSites.push_back(VectorInt3(1,1,2));

	contactSites.push_back(VectorInt3(-1,0,0));
	contactSites.push_back(VectorInt3(-1,1,0));
	contactSites.push_back(VectorInt3(-1,0,1));
	contactSites.push_back(VectorInt3(-1,1,1));
	contactSites.push_back(VectorInt3(0,-1,0));
	contactSites.push_back(VectorInt3(1,-1,0));
	contactSites.push_back(VectorInt3(0,-1,1));
	contactSites.push_back(VectorInt3(1,-1,1));
	contactSites.push_back(VectorInt3(0,0,-1));
	contactSites.push_back(VectorInt3(1,0,-1));
	contactSites.push_back(VectorInt3(0,1,-1));
	contactSites.push_back(VectorInt3(1,1,-1));

	VectorInt3 pos;
	double Energy=0.0;

	//loop through all polymer
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();
		pos=ingredients.getMolecules()[n];


		for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
		{
			int32_t latticeEntry=int32_t(ingredients.getLatticeEntry(pos+contactSites[contactNo]));

			Energy += getNNInteraction(monoType,latticeEntry);
		}
	}

	// factor 0.5 due to the double sum in hamiltonian
	return 0.5*Energy;
}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
double  FeatureWangLandauNextNeighbor<LatticeClassType>::getInternalEnergyDifference(const IngredientsType& ingredients,size_t index,VectorInt3 direction) const
{
	std::vector<VectorInt3> contactSites;

		//prepare the contact shell, in which the contacts are counted

		contactSites.push_back(VectorInt3(2,0,0));
		contactSites.push_back(VectorInt3(2,1,0));
		contactSites.push_back(VectorInt3(2,0,1));
		contactSites.push_back(VectorInt3(2,1,1));
		contactSites.push_back(VectorInt3(0,2,0));
		contactSites.push_back(VectorInt3(1,2,0));
		contactSites.push_back(VectorInt3(0,2,1));
		contactSites.push_back(VectorInt3(1,2,1));
		contactSites.push_back(VectorInt3(0,0,2));
		contactSites.push_back(VectorInt3(1,0,2));
		contactSites.push_back(VectorInt3(0,1,2));
		contactSites.push_back(VectorInt3(1,1,2));

		contactSites.push_back(VectorInt3(-1,0,0));
		contactSites.push_back(VectorInt3(-1,1,0));
		contactSites.push_back(VectorInt3(-1,0,1));
		contactSites.push_back(VectorInt3(-1,1,1));
		contactSites.push_back(VectorInt3(0,-1,0));
		contactSites.push_back(VectorInt3(1,-1,0));
		contactSites.push_back(VectorInt3(0,-1,1));
		contactSites.push_back(VectorInt3(1,-1,1));
		contactSites.push_back(VectorInt3(0,0,-1));
		contactSites.push_back(VectorInt3(1,0,-1));
		contactSites.push_back(VectorInt3(0,1,-1));
		contactSites.push_back(VectorInt3(1,1,-1));

		std::vector<VectorInt3> vicinitySites;
		vicinitySites.push_back(VectorInt3(0,0,0));
		vicinitySites.push_back(VectorInt3(1,0,0));
		vicinitySites.push_back(VectorInt3(0,1,0));
		vicinitySites.push_back(VectorInt3(1,1,0));
		vicinitySites.push_back(VectorInt3(0,0,1));
		vicinitySites.push_back(VectorInt3(1,0,1));
		vicinitySites.push_back(VectorInt3(0,1,1));
		vicinitySites.push_back(VectorInt3(1,1,1));


		VectorInt3 posOld;
		VectorInt3 posNew;
		double dEnergy=0.0;

				int32_t monoType=ingredients.getMolecules()[index].getAttributeTag();

				// old contacts
				posOld=ingredients.getMolecules()[index];

				for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
				{
					int32_t latticeEntryOld=int32_t(ingredients.getLatticeEntry(posOld+contactSites[contactNo]));

					dEnergy -= getNNInteraction(monoType,latticeEntryOld);
				}

				// new contacts
				posNew=ingredients.getMolecules()[index]+direction;

				for(size_t contactNo=0;contactNo<contactSites.size();contactNo++)
				{
					//exclude the old position
					bool possibleCount = true;
					for(size_t vicinityNo=0; vicinityNo<vicinitySites.size();vicinityNo++)
					{
						if((posNew+contactSites[contactNo]) == (posOld+vicinitySites[vicinityNo]))
							possibleCount = false;
					}

					if(possibleCount == true)
					{
						int32_t latticeEntryNew=int32_t(ingredients.getLatticeEntry(posNew+contactSites[contactNo]));
						dEnergy += getNNInteraction(monoType,latticeEntryNew);
					}
				}



		return dEnergy;



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
void FeatureWangLandauNextNeighbor<LatticeClassType>::exportRead(FileImport< IngredientsType >& fileReader)
{
  typedef FeatureWangLandauNextNeighbor<LatticeClassType> my_type;
  fileReader.registerRead("!nn_interaction",new ReadNNInteraction<my_type>(*this));
}


/**
 * @class WriteNNInteractionType
 * @brief Handles BFM-file write command !nn_interaction_shell_type
 *
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteNNInteractionType:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
	WriteNNInteractionType(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteNNInteractionType(){}

    virtual void writeStream(std::ostream& strm);
};


/**
 * @brief Executes the routine to write \b !nn_interaction_shell_type.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteNNInteractionType<IngredientsType>::writeStream(std::ostream& stream)
{

	stream<<"## nearest neighbor interactions shell type: NNShell (Hoffmann) or EShell (Rampf)\n";
	stream<<"!nn_interaction_shell_type=NNShell\n";

	stream<<"\n\n";

}

/**
 * @class WriteNNInteraction
 * @brief Handles BFM-file write command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteNNInteractionEnergy:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteNNInteractionEnergy(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}

    virtual ~WriteNNInteractionEnergy(){}

    virtual void writeStream(std::ostream& strm);
};


/**
 * @brief Executes the routine to write \b !nn_interaction.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteNNInteractionEnergy<IngredientsType>::writeStream(std::ostream& stream)
{

  stream<< std::endl << "## internal energy due to the next neighbor shell" << std::endl;

  stream<<"#!energy=" << this->getSource().getInternalEnergyCurrentConfiguration(this->getSource())<< std::endl;

  stream << std::endl;

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
void FeatureWangLandauNextNeighbor<LatticeClassType>::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  typedef FeatureWangLandauNextNeighbor<LatticeClassType> my_type;
  fileWriter.registerWrite("!nn_interaction",new WriteNNInteraction<my_type>(*this));
  fileWriter.registerWrite("!nn_interaction_shell_type",new WriteNNInteractionType<my_type>(*this));
  fileWriter.registerWrite("#!energy",new WriteNNInteractionEnergy<IngredientsType>(fileWriter.getIngredients_()));
}




#endif

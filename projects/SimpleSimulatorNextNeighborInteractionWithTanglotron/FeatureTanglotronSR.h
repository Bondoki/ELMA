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

#ifndef LEMONADE_FEATURE_TANGLOTRON_SR_H
#define LEMONADE_FEATURE_TANGLOTRON_SR_H
/**
* @file
*
* @class FeatureTanglotronSR
*
* @brief apply potential for one directional rotation of tanglotron motor
* 
* @details If one of four rotor monomers of the tanglotron is moved by simulation
* Metropolis-criterion is calculated by torque and angle of the move.
*
* @tparam IngredientsType
*
**/

#include <iostream>
#include <cmath>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include "TanglotronMotor.h"
#include "WriteTorque.h"
#include "WriteTanglotron.h"
#include "ReadTorque.h"
#include "ReadTanglotron.h"


class FeatureTanglotronSR: public Feature
{
public:
    
    typedef LOKI_TYPELIST_1(FeatureAttributes) required_features_front;
    typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;
    
    //! monomer_extensions, gives information if monomer is part of tanglotron motors and special which part
    typedef LOKI_TYPELIST_1(TanglotronAttributeTag) monomer_extensions;
    
    /**
    * constructor
    *
    * @param variableCalculatePotential to decide if tanglotron potential should
    * be calculated or not (default true)
    */
    FeatureTanglotronSR():variableCalculatePotential(true) {}
    
    virtual ~FeatureTanglotronSR(){}
    
    /**
    * synchronize
    *
    * @details no duty yet
    *
    * @param ingredients A reference to the IngredientsType - mainly the system.
    */
    template<class IngredientsType> void synchronize(IngredientsType& ingredients)
    {
        //std::cout << "Number of Tanglotrons in Vector: " << getTanglotronMotors().size() << std::endl;
    }
    
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move);
    
    template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveBase& move){ return true; }
	
    /**
    * uses the private variable variableCalculatePotential
    */
	bool getCalculatePotential() const {return variableCalculatePotential;}
	
    /**
    * @brief set bool which decides if tanglotron potential is calculated or not
    *
    * use it in updatern to build ingredients without additional potential
    * 
    * @param yesorno bool
    */
    void setCalculatePotential(bool yesorno)
    {
        variableCalculatePotential=yesorno;
    }
    
    /**
    * uses the private variable Torque
    */
    double getTorque() const
    {
        return Torque;
    }
    
    /**
    * sets the private variable Torque
    * 
    * @param torque_ value of torque
    */
    void setTorque(double torque_)
    {
        Torque = torque_;
    }
    
    /**
    * gets indices of all tanglotron motors
    * uses private vector TanglotronMotorUnits
    */
    const std::vector<TanglotronMotor>& getTanglotronMotors() const
    {
        return TanglotronMotorUnits;
    }
    
    /**
    * add one tanglotron motor to TanglotronMotorUnits
    * 
    * @param oneTanglotronMotor_ one TanglotronMotor unit
    * contains indices of one stator and four rotors
    */
    void addTanglotronMotor(TanglotronMotor oneTanglotronMotor_)
    {
        TanglotronMotorUnits.push_back(oneTanglotronMotor_);
    }
    
    /**
    * like getTanglotronMotors() but without const
    */
    std::vector<TanglotronMotor>& modifyTanglotronMotors()
    {
        return TanglotronMotorUnits;
    }
    
    /**
    * returns size() of TanglotronMotorUnits
    */
    const size_t getTanglotronMotorsSize() const
    {
        return TanglotronMotorUnits.size();
    }
    
    template<class IngredientsType> 
    void exportRead(FileImport<IngredientsType>& fileReader);
    
    template<class IngredientsType> 
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
    
    
private:
    
    //! to decide if potential should be calculated or not
    bool variableCalculatePotential;
    
    //! torque
    double Torque;
    //! indices of stators and rotors of all tanglotron motors in ingredients
    std::vector<TanglotronMotor> TanglotronMotorUnits;
    
    
    template<class IngredientsType> 
    double calculatePotential(const IngredientsType& ingredients, MoveLocalSc& move);
    
    int sgn(double value);
        
};


/**
* @brief checkMove by probability from calculatePotential
* 
* checks if move is allowed by the Metropolis-criterion.
* 
* @param move MoveLocalSc move
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template<class IngredientsType>
bool FeatureTanglotronSR::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) 
{	
    if((ingredients.getMolecules()[move.getIndex()].getIsTanglotron()) && (ingredients.getMolecules()[move.getIndex()].getTanglotronType()!=TanglotronAttributeTag::stator)){
        if(ingredients.getCalculatePotential())
            move.multiplyProbability(calculatePotential(ingredients, move));
    }
    
	return true;
}


/**
* @brief calculate potential for Metropolis-criterion
* 
* Look if moved monomer is a rotor monomer and save position of it as well as
* position of stator and the corresponding second rotor monomer.
* Then calculate vectors between stator and each rotor, axis, radius before and
* after move in plane perpendicular to axis and angle between this two radii.
* Then calculate Metropolis-criterion with angle and torque and return it.
* 
* @param move MoveLocalSc move
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template<class IngredientsType> 
double FeatureTanglotronSR::calculatePotential(const IngredientsType& ingredients, MoveLocalSc& move)
{
    //std::cout << "begin to calcuate potential" << std::endl;
    
    //index of moved monomer
    uint32_t movedMonomerIndex(move.getIndex());
    //ID of moved tanglotron motor
    uint32_t motorID(ingredients.getMolecules()[movedMonomerIndex].getTanglotronID());
    //we need two rotors -> use a vector
    std::vector <VectorInt3> rotorsPosition;
    //we need one stator
    VectorInt3 statorPosition(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexStator()]);
    
    
    //add first moved rotor monomer to rotorsPosition
    rotorsPosition.push_back(ingredients.getMolecules()[movedMonomerIndex]);
    
    //add the second rotor monomer belongig to the moved rotor monomer
    switch(ingredients.getMolecules()[movedMonomerIndex].getTanglotronType()) // either 1=rotorA, 2=rotorB, 3=rotorC, 4=rotorD
    {
        case TanglotronAttributeTag::rotorA:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorB()]);
            break;
        
        case TanglotronAttributeTag::rotorB:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorA()]);
            break;
        
        case TanglotronAttributeTag::rotorC:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorD()]);
            break;
        
        case TanglotronAttributeTag::rotorD:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorC()]);
            break;
        
        default:
            throw std::runtime_error("FeatureTanglotronSR: moved monomer is no rotor monomer!");
    }
    

    //if still here we have two rotor positions in rotorsPosition and one statorPosition and can calculate the potential
    
    //vectors stator->rotor
    //moved vector (according to vector b)
    VectorDouble3 statorRotor1(rotorsPosition.at(0)-statorPosition);
    //non-moved vector (according to vector a)
    VectorDouble3 statorRotor2(rotorsPosition.at(1)-statorPosition);
    
    //axis along statorRotor2 according to normalized vector a
    VectorDouble3 axis(statorRotor2.normalize());
    
    //radius of rotation of vector statorRotor1 according to vector r
    VectorDouble3 radius(statorRotor1-(statorRotor1*axis)*axis);
    //new radius according to vector r' (always constant rotational axis along a)
    VectorDouble3 radiusNew((statorRotor1+move.getDir())-((statorRotor1+move.getDir())*axis)*axis);
    
    //angle between old and new radius
    double angleCos((radius*radiusNew)/(radius.getLength()*radiusNew.getLength()));
    double angle;
    if(radiusNew.getLength()==0) //case for b'||omega => b' comes up to the constant rotationaly axis, but didn't circle in the plane perpendicular to the axis
        return 1.0; //angle = 0;
    else if(angleCos>1.0) //makes sure, that arccos isn't higher than 1 for rounding errors in case of r=r'
            return 1.0; //angle = 0; //acos(1.0);
    else if(angleCos<-1.0) //makes sure, that arccos isn't higher than 1 for rounding errors in case of r=-r'
            angle=M_PI; //acos(-1.0);
    else
        angle=acos(angleCos); //in every non-special case calculate angle according to Delta phi
    
    //normal on plane made of vectors statorRotor2 and statorRotor1 according to vector n=(axb)
    VectorDouble3 normal(crossProduct(statorRotor2,statorRotor1).normalize());
    
    //calculate potential according e^-Delta U
    return (exp(-(sgn(move.getDir()*normal) * Torque * angle)));
    
}


/**
* calculate signum function
* 
* @param value value to decide if it is above or below 0
*/
int FeatureTanglotronSR::sgn(double value)
{
    if(value>0)
        return 1;
    else if(value==0)
        return 0;
    else
        return -1;
}


/**
* read value of torque and indices of tanglotron motors from bfm file
* by function ReadTorque and ReadTanglotron
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template < class IngredientsType >
void FeatureTanglotronSR::exportRead(FileImport<IngredientsType>& fileReader)
{
    IngredientsType& destination=fileReader.getDestination();
    fileReader.registerRead("#!torque", new  ReadTorque< IngredientsType > (destination));
    fileReader.registerRead("#!tanglotrons", new  ReadTanglotron< IngredientsType > (destination));

}


/**
* write value of torque and indices of tanglotron motors in bfm file
* by function WriteTanglotronTorque
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template < class IngredientsType >
void FeatureTanglotronSR::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
    const IngredientsType& source=fileWriter.getIngredients_();
    fileWriter.registerWrite("#!torque", new WriteTorque <IngredientsType> (source));
    fileWriter.registerWrite("#!tanglotrons", new WriteTanglotron <IngredientsType> (source));
}


#endif //LEMONADE_FEATURE_TANGLOTRON_SR_H

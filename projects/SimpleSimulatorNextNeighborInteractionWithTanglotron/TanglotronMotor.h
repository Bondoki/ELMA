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

#ifndef LEMONADE_TANGLOTRON_MOTOR_H
#define LEMONADE_TANGLOTRON_MOTOR_H

#include <iostream>
#include <cmath>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureAttributes.h>


/**
* @class TanglotronAttributeTag
* 
* @brief Extension for tanglotrons - isTanglotron and tanglotronType
* 
* @details Extends tanglotrons by a bool as tag isTanglotron along with getter and
* setter. Initially this tag is set to false.
* Besides another enum tag tanglotronType is set, which holds information about
* stator or rotorA,B,C,D (int 0-4). Also this goes along with getter and setter
* functions. This tag is initially set to 5.
**/
class TanglotronAttributeTag
{
public:

    //! enum variable holding the tanglotronType
    enum tanglotronType{
        stator=0, 
        rotorA=1, 
        rotorB=2,
        rotorC=3,
        rotorD=4
    };
   
    //! Standard constructor- initially the tag is set to 5 and isTanglotron is set to false
    TanglotronAttributeTag():isTanglotron(false),tag(5),tanglotronID(){}

    /**
    * @brief setting boolean tag isTanglotron
    *
    * @param isTanglotron_ bool which tells if monomer is part of a tanglotron
    * motor or not
    **/
    void setIsTanglotron(bool isTanglotron_){ isTanglotron=isTanglotron_; }
	
    //! getting boolean tag isTanglotron
    bool getIsTanglotron() const {return isTanglotron;}
    
    /**
    * @brief setting enum tag tanglotronType
    *
    * @param attribute role of tanglotron part as string (stator, rotorA,B,C,D) or
    * integer (0,1,2,3,4)
    **/
    void setTanglotronType(int attribute){ tag=attribute; }
    
    //! getting enum tag tanglotronType
    int32_t getTanglotronType() const {if(isTanglotron) return tag;}
	
    /**
    * @brief setting size_t tag tanglotronID
    *
    * @param id number of motor in vector TanglotronMotorUnits
    **/
    void setTanglotronID(size_t id){ tanglotronID=id; }
    
    //! getting size_t tag tanglotronID
    size_t getTanglotronID() const {if(isTanglotron) return tanglotronID;}
	
    
private:
    //! Private Variable holding the bool if monomer is part of tanglotron
    bool isTanglotron;
    //! Private variable holding the tanglotronType. Default is 5.
    int32_t tag;
    //! private variable holding the number of one tanglotron motor in the vector holding all motors
    size_t tanglotronID;
};



/**
* @file
*
* @class TanglotronMotor
*
* @brief set and get indices of stator and rotors of a tanglotron motor
* 
* @details One tanglotron motor is built out of one central stator monomer and
* four surrounding rotor monomers which are connected each with the stator
* and in pairs with each other (A and B as well as C and D).
*
**/
class TanglotronMotor: public Feature
{
public:
    
    //! Standard constructor- initially all entries are set to NULL
    TanglotronMotor():StatorIndex(0), RotorAIndex(0), RotorBIndex(0), RotorCIndex(0), RotorDIndex(0), LK(0), LK_progress(), angleSum(0.0), angle_progress(), angleSumAB(0.0), angleSumCD(0.0) {}


    void setTanglotronMotor(uint32_t idxStator, uint32_t idxA, uint32_t idxB, uint32_t idxC, uint32_t idxD);

    uint32_t getIndexStator() const;
    uint32_t getIndexRotorA() const;
    uint32_t getIndexRotorB() const;
	uint32_t getIndexRotorC() const;
    uint32_t getIndexRotorD() const;
    
    void addLK(int8_t halfLKtoAddOrSubstract);
    
    int64_t getLK() const;
    std::vector<int8_t> getLKprogress() const;
    
    void addAngle(double angle);
    void addAngleAB(double angle);
    void addAngleCD(double angle);
    
    double getAngleSum() const;
    std::vector<double> getAngleSumProgress() const;
    double getAngleSumAB() const;
    double getAngleSumCD() const;
    
    
private:
    
    //! Index of the stator monomer of the tanglotron motor
    uint32_t StatorIndex;
    //! Index of the first (A) rotor monomer of the tanglotron motor
    uint32_t RotorAIndex;
    //! Index of the second (B) rotor monomer of the tanglotron motor
    uint32_t RotorBIndex;
    //! Index of the third (C) rotor monomer of the tanglotron motor
    uint32_t RotorCIndex;
    //! Index of the fourth (D) rotor monomer of the tanglotron motor
    uint32_t RotorDIndex;
    
    //! Linking Number, which will be counted on the run
    int64_t LK;
    //! Vector, which holds all added or substracted half Linking Numbers stepwise
    std::vector<int8_t> LK_progress;
    
    //! Linking Number, which will be counted on the run
    double angleSum;
    //! Vector, which holds current angleSum stepwise
    std::vector<double> angle_progress;
    //! motor angle, which will be counted on the run
    double angleSumAB;
    //! motor angle, which will be counted on the run
    double angleSumCD;
       
};



/**
* function to set the indices of the tanglotron motor unit
*/
void TanglotronMotor::setTanglotronMotor(uint32_t idxStator, uint32_t idxA, uint32_t idxB, uint32_t idxC, uint32_t idxD)
{
    StatorIndex = idxStator;
    RotorAIndex = idxA;
    RotorBIndex = idxB;
    RotorCIndex = idxC;
    RotorDIndex = idxD;
}


/**
* functions to get the indices of the parts of the tanglotron motor
*/
uint32_t TanglotronMotor::getIndexStator() const {
	return StatorIndex;
}

uint32_t TanglotronMotor::getIndexRotorA() const {
	return RotorAIndex;
}

uint32_t TanglotronMotor::getIndexRotorB() const {
	return RotorBIndex;
}

uint32_t TanglotronMotor::getIndexRotorC() const {
	return RotorCIndex;
}

uint32_t TanglotronMotor::getIndexRotorD() const {
	return RotorDIndex;
}


/**
* function to add or substract half Linking Numbers
*/
void TanglotronMotor::addLK(int8_t halfLKtoAddOrSubstract)
{   std::cout << "in addLK with " << halfLKtoAddOrSubstract << std::endl;
	LK += halfLKtoAddOrSubstract;
    LK_progress.push_back(halfLKtoAddOrSubstract);
}


/**
* functions to get the Linking Number of the tanglotron motor and its progress
*/
int64_t TanglotronMotor::getLK() const {
	return LK;
}

std::vector<int8_t> TanglotronMotor::getLKprogress() const {
	return LK_progress;
}


/**
* function to add or substract angle from angleSum
*/
void TanglotronMotor::addAngle(double angle)
{
	angleSum += angle;
    angle_progress.push_back(angle);
}


/**
* function to add or substract angleAB from angleSumAB
*/
void TanglotronMotor::addAngleAB(double angle)
{
	angleSumAB += angle;
}


/**
* function to add or substract angleCD from angleSumCD
*/
void TanglotronMotor::addAngleCD(double angle)
{
	angleSumCD += angle;
}


/**
* functions to get the angleSum of the tanglotron motor and its progress
*/
double TanglotronMotor::getAngleSum() const {
	return angleSum;
}

std::vector<double> TanglotronMotor::getAngleSumProgress() const {
	return angle_progress;
}


/**
* function to get the angleSumAB of the tanglotron motor
*/
double TanglotronMotor::getAngleSumAB() const {
	return angleSumAB;
}


/**
* function to get the angleSumCD of the tanglotron motor
*/
double TanglotronMotor::getAngleSumCD() const {
	return angleSumCD;
}


#endif //LEMONADE_TANGLOTRON_MOTOR_H

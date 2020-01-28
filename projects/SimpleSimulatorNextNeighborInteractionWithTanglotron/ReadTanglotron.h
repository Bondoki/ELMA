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

#ifndef LEMONADE_READ_TANGLOTRON_H
#define LEMONADE_READ_TANGLOTRON_H
/**
* @file
*
* @class ReadTanglotron
*
* @brief reading #!tanglotrons from bfm file
* 
* @details Form:
* #!tanglotrons
* indexStator indexRotorA indexRotorB indexRotorC indexRotorD
*
* @tparam IngredientsType
**/

#include <iostream>
#include <sstream>

#include <LeMonADE/io/AbstractRead.h>


template<class IngredientsType>
class ReadTanglotron: public ReadToDestination<IngredientsType>
{
    public:
        //! constructor
        ReadTanglotron(IngredientsType& destination): ReadToDestination<IngredientsType>(destination){}

        void execute();
};


/**
* ReadTanglotron::execute()
* 
* @brief reading tanglotrons and add them to private vector TanglotronMotorUnits
*/
template<class IngredientsType>
void ReadTanglotron<IngredientsType>::execute()
{
    std::cout << "reading #!tanglotrons ...\n";
    
    uint32_t idxStat, idxA, idxB, idxC, idxD;
    std::string line;
    std::streampos previous;
    
    std::getline(this->getInputStream(),line);
    previous=(this->getInputStream()).tellg();
    getline(this->getInputStream(),line);
    
    while(!line.empty() && !((this->getInputStream()).fail())){

        //stop at next Read and set the get-pointer to the position before the Read
        if(this->detectRead(line)){
            (this->getInputStream()).seekg(previous);
            break;
        }
  
        std::stringstream stream(line);
        
        //read index of stator, throw exception if extraction fails
        stream >> idxStat >> idxA >> idxB >> idxC >> idxD;
        if(stream.fail()){
            std::stringstream messagestream;
            messagestream << "ReadTanglotron<IngredientsType>::execute() -> Could not read indices";
            throw std::runtime_error(messagestream.str());
        }
        

        //if streaming worked set the TanglotronMotor with indices -1 because of
        //different index definitions of bfm and LeMonADE and add ist to private
        //vector TanglotronMotorUnits
        if(!stream.fail()){
            TanglotronMotor newTanglotron;
            newTanglotron.setTanglotronMotor(idxStat-1, idxA-1, idxB-1, idxC-1, idxD-1);
            this->getDestination().addTanglotronMotor(newTanglotron);
            
            this->getDestination().modifyMolecules()[idxStat-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronType(TanglotronAttributeTag::stator);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxA-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronType(TanglotronAttributeTag::rotorA);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxB-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronType(TanglotronAttributeTag::rotorB);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxC-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronType(TanglotronAttributeTag::rotorC);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxD-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronType(TanglotronAttributeTag::rotorD);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            
            std::getline(this->getInputStream(),line);
        }
        
        //otherwise throw an exception
        else{
            std::stringstream messagestream;
            messagestream << "ReadTanglotron<IngredientsType>::execute()\n" << "Could not read indices in #!tanglotrons";
            throw std::runtime_error(messagestream.str());
    
        }
    }    
}



#endif //LEMONADE_READ_TANGLOTRON_H

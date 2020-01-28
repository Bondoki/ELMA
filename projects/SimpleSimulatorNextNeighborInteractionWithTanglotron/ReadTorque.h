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

#ifndef LEMONADE_READ_TORQUE_H
#define LEMONADE_READ_TORQUE_H
/**
* @file
*
* @class ReadTorque
*
* @brief reading #!torque from bfm file
* 
* @details Form: #!torque=10
* 
* @tparam IngredientsType
**/

#include <iostream>
#include <sstream>

#include <LeMonADE/io/AbstractRead.h>


template<class IngredientsType>
class ReadTorque: public ReadToDestination<IngredientsType>
{
    public:
        //! constructor
        ReadTorque(IngredientsType& destination):ReadToDestination<IngredientsType>(destination){}
        
        void execute();
};


/**
* ReadTorque::execute()
* 
* @brief reading torque and set private variable Torque
*/
template<class IngredientsType>
void ReadTorque<IngredientsType>::execute()
{
    std::cout << "reading #!torque...";
    
    double torque;
    
    if(this->getInputStream() >> torque)
    { 
        this->getDestination().setTorque(torque);
        std::cout << torque << std::endl;
    }
    else
        throw std::runtime_error("ReadTorque<IngredientsType>::execute()\n Could not read torque");    
}


#endif //LEMONADE_READ_TORQUE_H

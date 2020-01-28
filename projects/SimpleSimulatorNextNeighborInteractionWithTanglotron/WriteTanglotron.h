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

#ifndef LEMONADE_WRITE_TANGLOTRON_H
#define LEMONADE_WRITE_TANGLOTRON_H
/**
* @file
*
* @class WriteTanglotron
*
* @brief writing #!tanglotrons to bfm file
* 
* @details Form:
* #!tanglotrons
* indexStator indexRotorA indexRotorB indexRotorC indexRotorD
*
* @tparam IngredientsType
**/

#include <iostream>

#include <LeMonADE/io/AbstractWrite.h>


template<class IngredientsType>
class WriteTanglotron: public AbstractWrite<IngredientsType>
{
    public:
        //! constructor
        WriteTanglotron(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}
  
        /**
        * writeStream
        * 
        * @brief getTanglotronMotors and write it to bfm file
        */ 
        void writeStream(std::ostream& strm)
        {
            const IngredientsType& ingredients=(this->getSource());
            
            //add 1 to indices of tanglotron monomers because of different index
            //definitions of bfm and LeMonADE
            strm << "#!tanglotrons\n";
            for(uint n=0; n<(this->getSource()).getTanglotronMotors().size(); n++){
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexStator()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorA()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorB()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorC()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorD()+1 << "\n";
            }
            strm << "\n";
        }
};


#endif //LEMONADE_WRITE_TANGLOTRON_H

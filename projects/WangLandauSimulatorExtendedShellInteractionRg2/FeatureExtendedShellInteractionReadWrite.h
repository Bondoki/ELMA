/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2020 by
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

#ifndef FEATURE_EXTENDED_SHELL_INTERACTION_READ_WRITE_H
#define FEATURE_EXTENDED_SHELL_INTERACTION_READ_WRITE_H

/**
 * @file
 * @date 2016/06/18, 2020/04/01
 * @author Hauke Rabbel, Ron Dockhorn
 * @brief Def. and impl. of class templates ReadExtendedShellInteraction and WriteExtendedShellInteraction
**/

#include <iostream>
#include <string>

#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>

/**
 * @class ReadExtendedShellInteraction
 * @brief Handles BFM-file read command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadExtendedShellInteraction: public ReadToDestination<IngredientsType>
{
public:
    ReadExtendedShellInteraction(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadExtendedShellInteraction(){}
    virtual void execute();
};

/**
 * @class WriteExtendedShellInteraction
 * @brief Handles BFM-file write command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteExtendedShellInteraction:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteExtendedShellInteraction(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteExtendedShellInteraction(){}

    virtual void writeStream(std::ostream& strm);
};

/**
 * @class ReadExtendedShellInteractionRadius
 * @brief Handles BFM-file read command !nn_interaction_shell_radius
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadExtendedShellInteractionRadius: public ReadToDestination<IngredientsType>
{
public:
    ReadExtendedShellInteractionRadius(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadExtendedShellInteractionRadius(){}
    virtual void execute();
};

/**
 * @class WriteExtendedShellInteractionRadius
 * @brief Handles BFM-file write command !nn_interaction_shell_radius
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteExtendedShellInteractionRadius:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteExtendedShellInteractionRadius(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteExtendedShellInteractionRadius(){}

    virtual void writeStream(std::ostream& strm);
};


/**
 * @class ReadExtendedShellInteractionType
 * @brief Handles BFM-file read command !nn_interaction_shell_type
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadExtendedShellInteractionType: public ReadToDestination<IngredientsType>
{
public:
    ReadExtendedShellInteractionType(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadExtendedShellInteractionType(){}
    virtual void execute();
};

/**
 * @class WriteExtendedShellInteractionType
 * @brief Handles BFM-file write command !nn_interaction_shell_type
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteExtendedShellInteractionType:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteExtendedShellInteractionType(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteExtendedShellInteractionType(){}

    virtual void writeStream(std::ostream& strm);
};

/////////////MEMBER IMPLEMENTATIONS ////////////////////////////////////////////

/**
 * @brief Executes the reading routine to extract \b !nn_interaction.
 *
 * @throw <std::runtime_error> fail to read monomer types or interaction energy.
 **/
template<class IngredientsType>
void ReadExtendedShellInteraction<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();

    int32_t typeA,typeB;
    double interactionConstant;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);
    
    try
      {
	file>>typeA;
	file>>typeB;
	file>>interactionConstant;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadExtendedShellInteraction::execute().\n";
	errormessage<<"Could not read interaction from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }
    
    //now save the interaction tuple just read from the file
    ingredients.setExtendedShellInteraction(typeA,typeB,interactionConstant);

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));
}

/**
 * @brief Executes the reading routine to extract \b !nn_interaction_shell_radius.
 *
 * @throw <std::runtime_error> fail to read monomer types or interaction energy.
 **/
template<class IngredientsType>
void ReadExtendedShellInteractionRadius<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();

    double nn_interactionRadius;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);

    try
      {
	file>>nn_interactionRadius;
	std::cout << "nn_interaction_radius: " << nn_interactionRadius <<std::endl;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadExtendedShellInteractionRadius::execute().\n";
	errormessage<<"Could not read interaction radius from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }

    //now save the interaction tuple just read from the file
    ingredients.setExtendedShellInteractionRadius(nn_interactionRadius);

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));
}

/**
 * @brief Executes the routine to write \b !nn_interaction.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteExtendedShellInteraction<IngredientsType>::writeStream(std::ostream& stream)
{
  int32_t nSpecies=255; //number must fit into 8 bit (uint8_t based lattice)
  stream<<"## nearest neighbor interactions between types in kT (default 0.0kT)\n";

  for(int32_t typeA=1;typeA<=nSpecies;typeA++)
    {
      for(int32_t typeB=1;typeB<=typeA;typeB++)
        {
	  if(this->getSource().getExtendedShellInteraction(typeA,typeB)!=0.0)
            {
	      stream<<"!nn_interaction "<<typeB<<" "<<typeA<<" "<<this->getSource().getExtendedShellInteraction(typeB,typeA)<<"\n";
            }

        }

    }
  stream<<"\n\n";

}

/**
 * @brief Executes the routine to write \b !nn_interaction_shell_radius.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteExtendedShellInteractionRadius<IngredientsType>::writeStream(std::ostream& stream)
{
  stream<<"## nearest neighbor interactions shell  2 <= r <= rmax\n";
  stream<<"!nn_interaction_shell_radius="<<this->getSource().getExtendedShellInteractionRadius()<<"\n";

  stream<<"\n\n";

}

/**
 * @class WriteNNInteraction
 * @brief Handles BFM-file write command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteExtendedShellInteractionEnergy:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteExtendedShellInteractionEnergy(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}

    virtual ~WriteExtendedShellInteractionEnergy(){}

    virtual void writeStream(std::ostream& strm);
};

/**
 * @brief Executes the routine to write \b !energy.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteExtendedShellInteractionEnergy<IngredientsType>::writeStream(std::ostream& stream)
{

  stream<< std::endl << "## internal energy due to the extended interaction shell" << std::endl;

  stream<<"#!energy=" << this->getSource().getInternalEnergyCurrentConfiguration(this->getSource())<< std::endl;

  stream << std::endl;

}


/**
 * @brief Executes the reading routine to extract \b !nn_interaction_shell_type.
 *
 * @throw <std::runtime_error> fail to read interaction types.
 **/
template<class IngredientsType>
void ReadExtendedShellInteractionType<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();

    std::string nn_interaction_type;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);

    try
      {
	file>>nn_interaction_type;
	std::cout << "nn_interaction_type: " << nn_interaction_type <<std::endl;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadExtendedShellInteractionType::execute().\n";
	errormessage<<"Could not read interaction type from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));

    //now save the interaction tuple just read from the file
    ingredients.setShellInteractionType(nn_interaction_type);


}

/**
 * @brief Executes the routine to write \b !nn_interaction_shell_type.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteExtendedShellInteractionType<IngredientsType>::writeStream(std::ostream& stream)
{
  stream<<"## nearest neighbor interactions shell type: NNShell (Hoffmann) or EShell (Rampf)\n";
  stream<<"!nn_interaction_shell_type="<<this->getSource().getShellInteractionType()<<"\n";

  stream<<"\n\n";

}

#endif // FEATURE_NN_INTERACTION_READ_WRITE_H

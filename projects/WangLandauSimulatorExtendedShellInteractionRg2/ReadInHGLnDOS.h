#ifndef READ_IN_HGLNDOS__H
#define READ_IN_HGLNDOS__H

#include <iostream>
#include <fstream>

#include <stdint.h>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <cstring>

#include "HistogramGeneralStatistik1D.h"

class ReadInHGLnDOS
{
public:

	ReadInHGLnDOS(double min, double max, size_t nbins, std::string _file_in)
	{

		HG_LnDOS.reset(min, max, nbins);

		file_in =_file_in;

	}

	HistogramGeneralStatistik1D& modifyHGLnDOS() {return HG_LnDOS;}
	const HistogramGeneralStatistik1D& getHGLnDOS() const {return HG_LnDOS;}

	virtual ~ReadInHGLnDOS()
	{
	}

	bool readin()
	{
		std::cout <<  " Reading file HGLnDOS " << std::endl;
		std::ifstream stream(file_in);

		std::string line, Read;
		std::streampos linestart;

		if (stream.is_open())
			while(!stream.eof() && !stream.fail()){
				//cout<<"findRead"<<endl;
				linestart=stream.tellg();
				getline(stream,line);
				if(!line.empty() && line.size()>1){
					bool ReadFound=(line.at(0)=='#' || (line.at(0)=='!'));

					// only comments or commands
					if(ReadFound){
						std::cout <<  line << std::endl;
						/*size_t length;
			//look for = sign
			length=line.find("=");
			if (length==std::string::npos){
				stream.seekg(linestart);
				stream>>Read;
				return Read;
			}
			else{
				stream.seekg(linestart);
				getline(stream,Read,'=');
				return Read;
			}
						 */

					}
					else
					{
						std::cout <<  line << std::endl;

						std::vector<std::string> entry = tokenize2Parameter(line, ' ', '\t');

						//entry has  2 tokens
						if((entry.size() != 2) )
							throw std::runtime_error("File HGLnDOS corrupted. EXITing\n");

						double energy = 0.0;
						double lnDOS = 0.0;

						std::istringstream ccenergy(entry[0]);
						ccenergy >> energy;

						//first entry is wrong
						if(ccenergy.fail()){
							throw std::runtime_error("Error in Energy in HGLnDOS file");
						}

						std::istringstream cclnDOS(entry[1]);
						cclnDOS >> lnDOS;

						//second entry is wrong
						if(cclnDOS.fail()){
							throw std::runtime_error("Error in LnDOS in HGLnDOS file");
						}

						//std::cout << energy << " -> " << lnDOS << std::endl;
						HG_LnDOS.addValue(energy, lnDOS);
					}
				}
			}

		stream.close();

		//if still here, the end of the file has been reached
		std::cout << "endoffile HGLnDOS" << std::endl ;

	}

	std::vector<std::string> tokenize1Parameter(const std::string& str,
				char delim) {
			std::vector<std::string> tokens;
			std::stringstream mySstream(str);
			std::string temp;

			while (getline(mySstream, temp, delim)) {
				tokens.push_back(temp);
			}

			return tokens;

		}

	std::vector<std::string> tokenize2Parameter(const std::string& str,
	  		char delim1, char delim2) {
	  	std::vector<std::string> tokens;
	  	std::stringstream mySstream(str);
	  	std::string temp;

	  	while (getline(mySstream, temp, delim1)) {

	  		if(temp.size() != 0)
	  		{
	  			std::stringstream mySstream2(temp);
	  			std::string temp2;
	  			while (getline(mySstream2, temp2, delim2)) {
	  				//printf("%s \n", temp2.c_str());
	  				if(temp2.size() != 0)
	  					tokens.push_back(temp2);
	  			}
	  		}

	  	}

	  	return tokens;

	  }

private:

	HistogramGeneralStatistik1D HG_LnDOS;

	std::string file_in;

};

#endif //READ_IN_HGLNDOS__H

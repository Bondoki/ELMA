#ifndef Histogram1D_H
#define Histogram1D_H

#include <iostream>

#include<stdint.h>
#include<stdexcept>
#include<vector>
#include<cmath>

class Histogram1D
{
public:

	Histogram1D()
	{
		nBins = 0;
		nValues= 0.0;
		valueSum=0.0;
		minVal=0.0;
		maxVal = 0.0;
		binsize = 0.0;
	}
	Histogram1D(double min, double max, int nbins)
	{
		
		nBins = nbins;
		nValues= 0.0;
		valueSum=0.0;
		minVal=min;
		maxVal = max;
		
		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins,0.0);
		
	}

	virtual ~Histogram1D()
	{
	}
	
	void reset(double min=0.0,double max=0.0,size_t nbins=0)
	{
		histogram.clear();
		nBins = nbins;
		nValues= 0.0;
		valueSum=0.0;
		minVal=min;
		maxVal = max;
		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins,0.0);
		
	}
	
	size_t getBinNo(double x) const
	{
		if(x<minVal)
			throw std::runtime_error("Histogram1D::getBinNo()...x value too low\n");
		else if(x>maxVal)
			throw std::runtime_error("Histogram1D::getBinNo()...x value too high\n");
		
		return uint32_t( std::floor((x-minVal)/binsize));
	}
	
	double getClosestBinCenter(double x) const
	{
		return (minVal+getBinNo(x)*binsize + binsize/2.0);
	}
	
	
	double getCountInBin(size_t bin) const
	{
		return histogram[bin];
	}
	
	double getCountAt(double x) const
	{
		return histogram[getBinNo(x)];
	}
	
	double getCountInBinNormalized(size_t bin) const
	{
		return histogram[bin]/valueSum;
	}
	
	double getCenterOfBin(int bin) const
	{
		return (minVal+(bin)*binsize + binsize/2.0);
	}
	
	int getNBins() const
	{
		return nBins;
	}
	
	double getNValues() const
	{
		return nValues;
	}
	
	double getBinwidth() const
	{
		return binsize;
	}
	
	double getMinCoordinate() const{return minVal;}
	double getMaxCoordinate() const{return maxVal;}

	void addValue(double x, double statisticalWeight=1.0) 
	{
		
		if( ( x > maxVal) || (x<minVal))
		{
			std::stringstream errormessage;
			errormessage<<"Histogram1D::addValue(). Value "<<x<<" outside boundaries "<<minVal<<" "<<maxVal<<std::endl;
			throw std::runtime_error(errormessage.str());	
		}
		
		histogram[getBinNo(x)]+=statisticalWeight;
		
		valueSum+=statisticalWeight;
		nValues+=1.0;
	}
	
	const std::vector<double>& getVectorValues() const
	{
		return histogram;
	}
	
	std::vector<double> getVectorValuesNormalized() const
	{
		
		std::vector<double> result=histogram;
		if(nValues>0.0)
		{
			for(size_t n=0;n<result.size();n++)
				result[n]/=valueSum;
		}
		return result;
	}
	
	
	std::vector<double> getVectorBins() const
	{
		std::vector<double> positions;
		positions.resize(nBins);
		
		for(size_t n=0;n<nBins;n++)
		{
			positions[n]=getCenterOfBin(n);
		}
		
		return positions;
	}
	
	
private:


	size_t nBins; 
	double nValues;
	double valueSum;
	double minVal; // untere Grenze
	double maxVal; // obere Grenze
	double binsize; //Intervalleinteilung
	std::vector<double> histogram;
};

#endif //Histogram1D.h

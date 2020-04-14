#ifndef HistogramGeneralStatistik1D_H
#define HistogramGeneralStatistik1D_H

#include <iostream>

#include <stdint.h>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "StatisticMoment.h"

class HistogramGeneralStatistik1D
{
public:

	HistogramGeneralStatistik1D()
	{
		nBins = 0;
		nValues= 0.0;

		minVal=0.0;
		maxVal = 0.0;
		binsize = 0.0;
	}
	HistogramGeneralStatistik1D(double min, double max, int nbins)
	{
		
		nBins = nbins;
		nValues= 0.0;

		minVal=min;
		maxVal = max;
		
		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins);
		
	}

	virtual ~HistogramGeneralStatistik1D()
	{
	}
	
	void reset(double min=0.0,double max=0.0,size_t nbins=0)
	{
		histogram.clear();

		nBins = nbins;
		nValues= 0.0;

		minVal=min;
		maxVal = max;

		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins);

		
	}
	
	size_t getBinNo(double x) const
	{
		if(x<minVal)
			throw std::runtime_error("HistogramGeneralStatistik1D::getBinNo()...x value too low\n");
		else if(x>maxVal)
			throw std::runtime_error("HistogramGeneralStatistik1D::getBinNo()...x value too high\n");
		
		return uint32_t( std::floor((x-minVal)/binsize));
	}
	
	double getClosestBinCenter(double x) const
	{
		return (minVal+getBinNo(x)*binsize + binsize/2.0);
	}
	
	
	double getFirstMomentInBin(size_t bin) const
	{
		return histogram[bin].ReturnM1();
	}
	
	double getSecondMomentInBin(size_t bin) const
		{
			return histogram[bin].ReturnM2();
		}

	double getNumCountInBin(size_t bin) const
			{
				return histogram[bin].ReturnN();
	}

	double getNumCountAt(double x) const
	{
		return histogram[getBinNo(x)].ReturnN();
	}
	
	double getFirstMomentAtKey(double x) const
		{
			return histogram[getBinNo(x)].ReturnM1();
		}

	double getCountAt(double x) const
		{
			return histogram[getBinNo(x)].ReturnM1();
		}
	
	double getCenterOfBin(int bin) const
	{
		return (minVal+(bin)*binsize + binsize/2.0);
	}
	
	double getLowerBoundaryOfBin(int bin) const
		{
			return (minVal+(bin)*binsize);
		}

	double getHigherBoundaryOfBin(int bin) const
			{
				return (minVal+(bin+1)*binsize);
			}

	int getNBins() const
	{
		return nBins;
	}
	
	
	double getBinwidth() const
	{
		return binsize;
	}
	
	double getMinCoordinate() const{return minVal;}
	double getMaxCoordinate() const{return maxVal;}

	//double getValueSum() const{return valueSum;}

	void addValue(double x, double statisticalWeight=1.0) 
	{
		
		if( ( x > maxVal) || (x<minVal))
		{
			std::stringstream errormessage;
			errormessage<<"HistogramGeneralStatistik1D::addValue(). Value "<<x<<" outside boundaries "<<minVal<<" "<<maxVal<<std::endl;
			throw std::runtime_error(errormessage.str());	
		}
		
		histogram[getBinNo(x)].AddValue(statisticalWeight);
		
		nValues+=1.0;
	}
	
	void resetValue(double x, double statisticalWeight=1.0)
	{

		if( ( x > maxVal) || (x<minVal))
		{
			std::stringstream errormessage;
			errormessage<<"HistogramGeneralStatistik1D::addValue(). Value "<<x<<" outside boundaries "<<minVal<<" "<<maxVal<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
		histogram[getBinNo(x)].clear();

		histogram[getBinNo(x)].AddValue(statisticalWeight);


	}

	
	const std::vector<StatisticMoment>& getVectorValues() const
		{
			return histogram;
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
	double minVal; // untere Grenze
	double maxVal; // obere Grenze
	double binsize; //Intervalleinteilung

	std::vector<StatisticMoment> histogram;
};

#endif //HistogramGeneralStatistik1D.h

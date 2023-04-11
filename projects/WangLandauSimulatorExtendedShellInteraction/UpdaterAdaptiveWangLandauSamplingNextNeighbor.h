#ifndef ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR__H
#define ADAPTIVE_UPDATER_WANG_LANDAU_SAMPLING_NEXT_NEIGHBOR__H

#include <limits>
#include <iomanip>

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/UpdaterAbstractCreate.h>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include "Histogram1D.h"

/* ****************************************************************************
 * Simulator for adaptive umbrella sampling along a (1D) reaction coordinate.
 * The simulator runs standart BFM, while during  the simulation recording a 
 * histogram of the reaction coordinate. 
 * 
 * When the histogram has converged by three given criteria (described below),
 * a bias potential is constructed from the histogram, the histogram is 
 * reset to zero, and the procedure is repeated with the bias potential. 
 * This procedure is carried out with updates on the bias potential based on 
 * collected histograms, until the histogram becomes flat. Once this is the case
 * the bias potential is the inverted free energy profile along the coordinate.
 * 
 * The coordinate itsself is not specified in this updater. It requires 
 * Ingredients to have an interface with the signature 
 * double getCurrentReactionCoordinate()
 * Such an interface is provided by the FeatureBilayerDistancePotential. In that
 * case the coordinate is the distance to the bilayer. Other features could be 
 * used, provided they have the above described interface.
 * 
 * For the convergence of the histogram, three criteria are used to judge whether
 * the histogram has (locally) converged, i.e. does not change much any more.
 * These criteria are:
 * 1) relative change of the mean value of the expectation value of the coordinate.
 * <x>=sum(histogram(x) * x).
 * 2) relative change of the variance of the histogram coordinate
 * <dx**2>=<x**2>-<x>**2
 * 3) relative change of the square deviation from the mean of the histogram
 * itsself   <dH**2>=< (h(x) - <h>)**2 >
 * For these checks, the normalized histogram is used, and it is checked whether
 * these values changed in the last two update trials by more than some threshold
 * percentage.
 * 
 * For the absolute convergence, the square deviation of the histogram from the
 * mean is required to be less than some absolute value ( e.g. 0.1)
 * 
 * The threshold values can be set by the appropriate interfaces or the 
 * constructor.
 * ****************************************************************************/

/**
 * @file
 *
 * @class UpdaterAdaptiveWangLandauSamplingNextNeighbor
 *
 * @brief Adaptive Umbrella Sampling updater along a reaction coordinate
 *
 * @details requires FeatureBilayerDistancePotential or other feature handling a 
 * possibly discrete potential along a reaction coordinate. this feature needs to
 * provide the functions: get/modifyPotential()  returning the potential object, 
 * getCurrentReactionCoordinate() returning the current value of the reaction
 * coordinate. also, the potential object needs to provide the functions 
 * addValue(coordinate, value) and getPMFMaxCoordinate(),getPMFMinCoordinate(),
 * getPMFValues() returning a std::vector<double> of the pmf values
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam <MoveType> name of the specialized move.
 */
template<class IngredientsType,class MoveType>
class UpdaterAdaptiveWangLandauSamplingNextNeighbor: public UpdaterAbstractCreate<IngredientsType>
{
	
public:

	 typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
	/**
	 * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
	 *
	 * @param ing a reference to the IngredientsType - mainly the system
	 * @param steps MCS per cycle to performed by execute()
	 * @param stepsBeforeBiasCheck interval between checks if histogram has converged im MCS
	 * @param stepsBeforeHistogramUpdate interval between updates on coordinate histogram in MCS
	 * @param convergenceThreshold threshold for relative convergence of histogram
	 * @param fullFlatnessThreshold threshold for absolute convergence of histogram
	 * @param maxAge maximum MCS limit. If the histogram converges earlier, the simulation also stops
	 */
	UpdaterAdaptiveWangLandauSamplingNextNeighbor(IngredientsType& ing
	,uint32_t steps
	,uint32_t stepsBeforeBiasCheck
	,double _flatness=0.85
	,double _initialModificationFactor=1.01
	,uint64_t maxAge=1000000000
	,double minWindow = -100.0
	,double maxWindow = +100.0
	,int numWindow=0
	,double _modificationFactorThreshold=std::exp(std::pow(10,-8))
	,double _modificationFactorThesholdUsing1t=0.0
	);
	
	/**
	 * @brief This checks all used Feature and applies all Feature if all conditions are met.
	 *
	 * @details This function runs over \a steps MCS and performs the moves.
	 * It setting the age of the system and prints a simple simple simulation speed
	 * in the number of attempted monomer moves per s (tried and performed monomer moves/s).
	 *
	 * @return True if function is done.
	 */
	bool execute();

	bool executeFindWindow();
	
	
	/**
	 * @brief This function is called \a once in the beginning of the TaskManager.
	 *
	 * @details It´s a virtual function for inheritance.
	 * Use this function for initializing tasks (e.g. init SDL)
	 *
	 **/
	virtual void initialize();
	
	/**
	 * @brief This function is called \a once in the end of the TaskManager.
	 *
	 * @details writes the final state of histogram, bias potential and the
	 * history of the convergence criteria to output files.
	 *
	 **/
	virtual void cleanup(){

		std::stringstream itr;
		itr << std::setw(2) << std::setfill('0') << iteration;

		dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_histogram.dat"));
		dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_totalhistogram.dat"));
		dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_iteration" + itr.str() +  "_final_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
		dumpConvergenceProgress();
		//updateModificationFactor();

		std::stringstream ssprefixWindowBFM;
		ssprefixWindowBFM << "_idxWin" << std::setw(4) << std::setfill('0') << numberIdxWindow;
		writeBFM_File(std::string(ingredients.getName() + ssprefixWindowBFM.str() + "_final.bfm"));
	};
	
	
	//! set the threshold for relative convergence of histogram (for update of bias)
	void setRelativeConvergenceThreshold(double threshold){convergenceThreshold=threshold;}
	
	//! return the threshold for relative convergence of histogram (for update of bias)
	double getRelativeConvergenceThreshold() const{return convergenceThreshold;}
	
	//! set the threshold for complete convergence of histogram (for final simulation convergence)
	void setFullConvergenceThreshold(double threshold){fullFlatnessThreshold=threshold;}
	
	//! return the threshold for complete convergence of histogram (for final simulation convergence)
	double getFullConvergenceThreshold() const {return fullFlatnessThreshold;}
	
	//! get the number of pmf updates performed up to now on the bias potential
	uint32_t getNPerformedUpdates()const{return nUpdatesPerformed;}
	
	//! check if the histogram has converged according to chosen criteria
	bool histogramConverged();

	void doResetForNextIteration()
	{
		updateModificationFactor();
		iteration++;
		std::cout << "Start Iteration" << iteration << std::endl;

		iterationconverged = false;
		iterationfirstconverged = false;
		//if the simulation has converged, end the simulation
		//if(simulationConverged==true) return false;

		//set histogram to 0
		resetHistogram();
	}

	void outputConvergedIteration()
	{

		//write some file output
		dumpConvergenceProgress();
		if(writeHistogramProgress==true){
			std::stringstream filename;
			filename<<ingredients.getName() << prefixWindow  << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_histogram_mcs"
					<<ingredients.getMolecules().getAge()
					<<".dat";
			dumpHistogram(filename.str());

			std::stringstream filenameTotal;
			filenameTotal<<ingredients.getName() << prefixWindow << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_totalhistogram_mcs"<<ingredients.getMolecules().getAge()<<".dat";
			dumpTotalHistogram(filenameTotal.str());
		}


		//write some more output
		if(writePMFProgress==true)
		{
			std::stringstream filenametmp;
			filenametmp<<ingredients.getName() << prefixWindow << "_iteration" << std::setw(2) << std::setfill('0') << iteration << "_HGLnDOS" << "_mcs"<<ingredients.getMolecules().getAge() << ".dat";

			dumpHGLnDOS(filenametmp.str(), ingredients.getMinWin(), ingredients.getMaxWin());

		}

		std::stringstream ssprefixWindowBFM;
		ssprefixWindowBFM << "_idxWin" << std::setw(4) << std::setfill('0') << numberIdxWindow;
		writeBFM_File(std::string(ingredients.getName() + ssprefixWindowBFM.str() + ".bfm"));


	}

	bool isIterationConverged() {return iterationconverged;};

	bool isFirstConverged() {return iterationfirstconverged;};

	void unsetFirstConverged() {iterationfirstconverged=false;};

	void checkForFirstWindowAppearance();
	
	int getIteration() {return iteration;};
	
	void setReachedFinalState() {reachedFinalState = true;};
	
	bool getReachedFinalState() {return reachedFinalState;};

private:
	using BaseClass::ingredients;

	using BaseClass::randomBondvector;

	//! write the currently collected histogram to a file named histogram_mcs.dat
	void dumpHistogram(std::string filePrefix);
	
	//! write the currently collected histogram to a file named histogram_mcs.dat
	void dumpTotalHistogram(std::string filePrefix);

	//! writes the time development of the mean,variance, and square deviation from average to file names convergence.dat
	void dumpConvergenceProgress();
	
	//! write the currently state to bfm file
	void writeBFM_File(std::string filePrefix);

	//! update the used bias potential according to collected histogram
	
	void updateModificationFactor();

	//! write the currently logarithm of Density of States DOS i.e. file named HGLnDOS.dat
	void dumpHGLnDOS(std::string prefix="", double min= std::numeric_limits<double>::min(), double max= std::numeric_limits<double>::max());

	
	//! reset the histogram to empty state
	void resetHistogram();
	
	//! check the mean square deviation from a flat histogram to check for full conversion of the simulation
	double histogramFlattness() const;
	
	//! check if mean and variance of histogram have converged. called as part of histogramConverged()
	bool distributionCheck();
	
	//! check if histogram flatness has converged, called as part of histogramConverged()
	bool flatnessCheck();
	
	//! set bin width of coordinate histogram and external potential
	bool setBinWidth();
	
	//! shift the histogram be its maximum entry
	void shiftHGLnDOS();
	
	//! System information
	//IngredientsType& ingredients;
	
	//! Specialized move to be used
	MoveType move;
	
	//! Max number of mcs to be executed
	uint32_t nsteps;
	
	//! Number of steps between checks for update of bias
	const uint32_t nStepsBeforeBiasCheck;
	
	//! Number of checks before update of histogram
	uint32_t nStepsBeforeHistogramUpdate;
	
	//! current counter to decide if update check is needed
	uint32_t counter_nStepsBeforeBiasCheck;
	
	//! current counter to decide if update check is needed
	uint32_t counter_nStepsBeforeHistogramUpdate;
	
	//! counts the number of bias updates performed so far
	uint32_t nUpdatesPerformed;
	
	//! histogram of reaction coordinate
	//Histogram1D histogram;
	
	//values for checking convergence
	double oldMean,oldOldMean;
	double oldVariance,oldOldVariance;
	double oldFlatnesValue,oldOldFlatnessValue;
	
	//thresholds for convergence
	double convergenceThreshold;
	double fullFlatnessThreshold; //for flatness
	double maxBinsThreshold;
	
	//time history of convergence criteria
	std::vector<double> varianceSeries;
	std::vector<double> meanSeries;
	std::vector<double> flatnessSeries;
	
	//! set to true if full convergence threshold has been reached
	bool simulationConverged;
	
	//! flags for setting or unsetting extra output of histogram
	bool writeHistogramProgress;
	
	//! flags for setting or unsetting extra output of bias potential
	bool writePMFProgress;
	
	//! maximum age of the system
	uint64_t maxSystemAge;
	
	double minWindow;
	double maxWindow;
	std::string prefixWindow;
	int iteration;

	double initialModificationFactor;

	bool iterationconverged;

	bool iterationfirstconverged;

	int numberIdxWindow;

	double flatness;
	
	double modificationFactorThreshold;
	
	bool reachedFinalState;
};


////////////////////////////////////////////////////////////////////////////////
// member implementations
///////////////////////////////////////////////////////////////////////////////


// constructor
template<class IngredientsType,class MoveType>
UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::
UpdaterAdaptiveWangLandauSamplingNextNeighbor(IngredientsType& ing,
				uint32_t steps, 
				uint32_t stepsBeforeBiasCheck,
				double _flatness,
				double _initialModificationFactor,
				uint64_t maxAge,
				double _minWindow,
				double _maxWindow,
				int numberWindow,
				double _modificationFactorThreshold,
				double _modificationFactorThesholdUsing1t
				)
://ingredients(ing),
nsteps(steps)
,counter_nStepsBeforeBiasCheck(0)
,counter_nStepsBeforeHistogramUpdate(0)
,nStepsBeforeBiasCheck(stepsBeforeBiasCheck)
,writeHistogramProgress(true)
,writePMFProgress(true)
,nUpdatesPerformed(0)
,oldMean(0.0)
,oldOldMean(0.0)
,oldVariance(0.0)
,oldOldVariance(0.0)
,oldFlatnesValue(10.0)
,oldOldFlatnessValue(10.0)
,flatness(_flatness)
,initialModificationFactor(_initialModificationFactor)
,maxSystemAge(maxAge)
,minWindow(_minWindow)
,maxWindow(_maxWindow)
,iteration(0)
,simulationConverged(false)
,numberIdxWindow(numberWindow)
,modificationFactorThreshold(_modificationFactorThreshold)
,BaseClass(ing)
{
	nStepsBeforeHistogramUpdate=1000;
	convergenceThreshold=0.02;
	fullFlatnessThreshold=0.1;
	maxBinsThreshold=10000.0;

	varianceSeries.resize(0);
	meanSeries.resize(0);
	flatnessSeries.resize(0);
	ingredients.setModificationFactor(ingredients, initialModificationFactor);
	ingredients.setWindowState(false, ingredients.getVisitsEnergyStates().getMinCoordinate(), ingredients.getVisitsEnergyStates().getMaxCoordinate());

	iterationconverged=false;
	iterationfirstconverged = false;

	std::stringstream ssprefixWindow;
	ssprefixWindow << "_idxWin" << std::setw(4) << std::setfill('0') << numberWindow;
	prefixWindow = ssprefixWindow.str();

	reachedFinalState = false;
	
	ingredients.setModificationFactorUsing1t(_modificationFactorThesholdUsing1t);
}

//execute simulation
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::executeFindWindow()
{
	//some output for simulation performance and state
	time_t startTimer = std::time(NULL); //in seconds
	uint64_t systemAge=ingredients.getMolecules().getAge();
	std::cout<<"mcs "<<systemAge 
		 << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) <<std::endl;
	
	//simulation loop
	if(writeHistogramProgress==true)
				{
					//shiftHGLnDOS();

					 dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_histogram.dat"));
					 //dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
					 dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_totalhistogram.dat"));
					 //std::cout << "dump HGLnDOS " << ingredients.getMinWin() << "   " << ingredients.getMaxWin() << std::endl;
				}
				
	//shiftHGLnDOS();
	
	resetHistogram();
	
	
	for(int n=0;n<nsteps;n++){
		
		//update bias potential if the histogram has converged well enough
		if(counter_nStepsBeforeBiasCheck==nStepsBeforeBiasCheck){
			counter_nStepsBeforeBiasCheck=0;
			if(writeHistogramProgress==true)
				{
					//shiftHGLnDOS();

					 //dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_histogram.dat"));
					 dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
					 std::cout << "dump HGLnDOS " << ingredients.getMinWin() << "   " << ingredients.getMaxWin() << std::endl;
				}
			

		}
		
		
		//std::cout << "1MCS";
		//do one monte carlo sweep
		for(int m=0;m<ingredients.getMolecules().size();m++)
		{
			move.init(ingredients);
			
			if(move.check(ingredients)==true)
			{
				move.apply(ingredients);
			}
			else 
			{ 
				ingredients.rejectMove(ingredients);
			}
		}

		checkForFirstWindowAppearance();
		


		//get update on histogram after every sweep. more often would 
		//probably not make much sense, because the chain diffuses slowly
		try{
			if(counter_nStepsBeforeHistogramUpdate==nStepsBeforeHistogramUpdate){
				counter_nStepsBeforeHistogramUpdate=0;
				//histogram.addValue(ingredients.getInternalEnergy(ingredients));
			}
			
		}
		catch(std::runtime_error& e){
			std::stringstream errormessage;
			errormessage<<"UpdaterAdaptiveWangLandauSamplingNextNeighbor: error while altering histogram.\n"
			<<"original message was:\n"<<e.what()<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
		
		//increment counter
		counter_nStepsBeforeBiasCheck++;
		counter_nStepsBeforeHistogramUpdate++;
		
	}
	ingredients.updateHGLnDOSBLENDER();
	
	
	//update age of system and give some more output on simulation progress
	ingredients.modifyMolecules().setAge(ingredients.modifyMolecules().getAge()+nsteps);
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with " << (((1.0*nsteps)*ingredients.getMolecules().size())/(std::difftime(std::time(NULL), startTimer)) ) << " [attempted moves/s]" <<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) << " with " << nsteps << " MCS "<<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with energy " << ingredients.getInternalEnergyCurrentConfiguration(ingredients) <<std::endl;

	return true;
}

 
//execute simulation
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::execute()
{
	//some output for simulation performance and state
	time_t startTimer = std::time(NULL); //in seconds
	uint64_t systemAge=ingredients.getMolecules().getAge();
	std::cout<<"mcs "<<systemAge 
		 << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) <<std::endl;
	
	//if max age has reached, stop
	if(systemAge>=maxSystemAge) {
		std::cout<<"max age has been reached. maxAge "<<maxSystemAge<<" system age "<<systemAge<<std::endl;
		return false;
	}
	
	if(simulationConverged)
		{
			std::cout<<"Simulation converged! " << std::endl;
			return false;
		}
	
	//simulation loop
	if(writeHistogramProgress==true)
				{
					//shiftHGLnDOS();

					 dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_histogram.dat"));
					 //dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
					 dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_totalhistogram.dat"));
					 //std::cout << "dump HGLnDOS " << ingredients.getMinWin() << "   " << ingredients.getMaxWin() << std::endl;
				}
				
	//shiftHGLnDOS();
	
	resetHistogram();
	//ingredients.initBLENDER();
	
	for(int n=0;n<nsteps;n++){
		
		//update bias potential if the histogram has converged well enough
		if(counter_nStepsBeforeBiasCheck==nStepsBeforeBiasCheck){
			counter_nStepsBeforeBiasCheck=0;
			if(writeHistogramProgress==true)
				{
					//shiftHGLnDOS();

					 //dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_tmp_histogram.dat"));
					 dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_tmp_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());
					 std::cout << "dump HGLnDOS " << ingredients.getMinWin() << "   " << ingredients.getMaxWin() << std::endl;
				}
			

		}
		
		
		//std::cout << "1MCS";
		//do one monte carlo sweep
		for(int m=0;m<ingredients.getMolecules().size();m++)
		{
			move.init(ingredients);
			
			if(move.check(ingredients)==true)
			{
				move.apply(ingredients);
			}
			else 
			{ 
				ingredients.rejectMove(ingredients);
			}
		}
		//std::cout << " done" << std::endl;
		
		//if( (ingredients.isEnergyInWindow() == false ) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) < maxWindow+2*ingredients.getHGLnDOS().getBinwidth()) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) > minWindow-2*ingredients.getHGLnDOS().getBinwidth()) )
		//checkForFirstWindowAppearance();





		//get update on histogram after every sweep. more often would 
		//probably not make much sense, because the chain diffuses slowly
		try{
			if(counter_nStepsBeforeHistogramUpdate==nStepsBeforeHistogramUpdate){
				counter_nStepsBeforeHistogramUpdate=0;
				//histogram.addValue(ingredients.getInternalEnergy(ingredients));
			}
			
		}
		catch(std::runtime_error& e){
			std::stringstream errormessage;
			errormessage<<"UpdaterAdaptiveWangLandauSamplingNextNeighbor: error while altering histogram.\n"
			<<"original message was:\n"<<e.what()<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
		
		//increment counter
		counter_nStepsBeforeBiasCheck++;
		counter_nStepsBeforeHistogramUpdate++;
		
	}
	//ingredients.updateHGLnDOSBLENDER();
	
	checkForFirstWindowAppearance();
	//update age of system and give some more output on simulation progress
	ingredients.modifyMolecules().setAge(ingredients.modifyMolecules().getAge()+nsteps);
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with " << (((1.0*nsteps)*ingredients.getMolecules().size())/(std::difftime(std::time(NULL), startTimer)) ) << " [attempted moves/s]" <<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((std::difftime(std::time(NULL), startTimer)) ) << " with " << nsteps << " MCS "<<std::endl;
	std::cout<<"mcs "<<ingredients.getMolecules().getAge() << " with energy " << ingredients.getInternalEnergyCurrentConfiguration(ingredients) <<std::endl;


	if(ingredients.getModificationFactor(ingredients) < modificationFactorThreshold )
	{
		std::cout<<"Simulation converged! " << std::endl;
		return false;
	}
	else
	{
		std::cout<<"Simulation approach convergence.. " << std::endl;
		return true;
	}
}



// initialize only sets the histogram and bias potential to 0
template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::initialize()
{
	resetHistogram();
	dumpHistogram(std::string(ingredients.getName() + prefixWindow +  "_initial_histogram.dat"));
	dumpTotalHistogram(std::string(ingredients.getName() + prefixWindow +  "_initial_totalhistogram.dat"));
	//dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_initial"), ingredients.getHGLnDOS().getMinCoordinate(), ingredients.getHGLnDOS().getMaxCoordinate());
	dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_initial"), ingredients.getHGLnDOS().getMinCoordinate(), ingredients.getHGLnDOS().getMaxCoordinate());

	writeBFM_File(std::string(ingredients.getName() + prefixWindow + "_initial.bfm"));
}

// check if system is in window and reset all histograms
template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::checkForFirstWindowAppearance()
{
	if( (ingredients.isEnergyInWindow() == false ) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) <= maxWindow) && (ingredients.getInternalEnergyCurrentConfiguration(ingredients) >= minWindow) )
			{
				std::cout << "RW is in energy window: [" << minWindow << " ; " << maxWindow << "] with " <<  ingredients.getInternalEnergyCurrentConfiguration(ingredients) << std::endl;
				// RW is in energy window
				ingredients.setWindowState(true, minWindow, maxWindow);

				std::stringstream ssprefixWindow;
				ssprefixWindow << "_idxWin" << std::setw(4) << std::setfill('0') << numberIdxWindow << "_minWin" << minWindow << "_maxWin" <<  maxWindow;

				prefixWindow = ssprefixWindow.str();

				std::stringstream ssprefixWindowBFM;
				ssprefixWindowBFM << "_idxWin" << std::setw(4) << std::setfill('0') << numberIdxWindow;
				writeBFM_File(std::string(ingredients.getName() + ssprefixWindowBFM.str() + ".bfm"));

				//dump histograms in the start of iteration procedure
				dumpHistogram(std::string(ingredients.getName() + prefixWindow + "_startIteration_histogram.dat"));
				dumpHGLnDOS(std::string(ingredients.getName() + prefixWindow + "_startIteration_HGLnDOS.dat"), ingredients.getMinWin(), ingredients.getMaxWin());

				// delete all histograms

				ingredients.modifyVisitsEnergyStates().reset(ingredients.getVisitsEnergyStates().getMinCoordinate(),ingredients.getVisitsEnergyStates().getMaxCoordinate(),ingredients.getVisitsEnergyStates().getNBins());
				ingredients.modifyTotalVisitsEnergyStates().reset(ingredients.getTotalVisitsEnergyStates().getMinCoordinate(),ingredients.getTotalVisitsEnergyStates().getMaxCoordinate(),ingredients.getTotalVisitsEnergyStates().getNBins());

				// rest HGLnDOS to avoid overshoot at boundaries
				// as the WL-BLENDER algorithm converges to self-consitent solution you need to clear the initial guess
				ingredients.modifyHGLnDOS().reset(ingredients.getHGLnDOS().getMinCoordinate(),ingredients.getHGLnDOS().getMaxCoordinate(),ingredients.getHGLnDOS().getNBins());

			}
}

template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::resetHistogram()
{
	size_t nbins=ingredients.getVisitsEnergyStates().getNBins();
	double min = ingredients.getVisitsEnergyStates().getMinCoordinate();
	double max = ingredients.getVisitsEnergyStates().getMaxCoordinate();
	//histogram.reset(-10000,0.4,nbins);
	ingredients.modifyVisitsEnergyStates().reset(min,max,nbins);
}

template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::shiftHGLnDOS()
{
	std::cout << "Shift HgLnDOS ";

	double eln = ingredients.getHGLnDOS().getFirstMomentInBin(0);
	//size_t bin = 0;
			//for (size_t n=1; n < ingredients.getHGLnDOS().getNBins(); n++) {
				//if ((ingredients.getHGLnDOS().getFirstMomentInBin(n) < eln && ingredients.getHGLnDOS().getFirstMomentInBin(n) != 0 ) || eln == 0)
				//{
					//eln = ingredients.getHGLnDOS().getFirstMomentInBin(n);
					//bin = n;
				//}
			//}

	//std::cout << " by " << eln << " (" << ingredients.getHGLnDOS().getFirstMomentInBin(bin) << ") at bin " << bin << " with energy "  << ingredients.getHGLnDOS().getCenterOfBin(bin) << std::endl;

	// this is an arbitrary shift to avoid extrem high numbers
	// the shift value can be set e.g. minimal energy entry
	eln = ingredients.getHGLnDOS().getCountAt(ingredients.getMinWin()+ingredients.getHGLnDOS().getBinwidth());
	std::cout << " by " << eln << " with energy "  << (ingredients.getMinWin()+ingredients.getHGLnDOS().getBinwidth()) << std::endl;
	
	
	for (size_t n=0; n < ingredients.getHGLnDOS().getNBins(); n++)
	{
		if(ingredients.getHGLnDOS().getVectorValues()[n].ReturnN() != 0)
		{
			double centerEnergy = ingredients.getHGLnDOS().getCenterOfBin(n);
			ingredients.modifyHGLnDOS().resetValue(centerEnergy, ingredients.getHGLnDOS().getCountAt(centerEnergy)-eln);
		}

	}

}



//checks if mean, variance and deviation from average of the histogram have
//changed more than the percentage threshold. If all has converged, the function
//returns true
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::histogramConverged()
{
	// only check for flat histogram if RW is in desired energy space
	if(ingredients.isEnergyInWindow() == false)
		return false;

	if(iterationconverged == true)
		return true;
	
	return false;

	std::cout<<"Check histogram flatness" << std::endl;

	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValues();//histogram.getVectorValues();

	// check if we have an entry for every element within the window:
	// avoiding 'artifical' jumps in LnDOS do to inpropriate sampling
	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		if( (ingredients.getVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
		{
			if(currentHistogramState.at(n) == 0.0)
			{
				return false;
			}
		}
	}
	
	// otherwise calculate average flatness
	double mean= 0.0;
	int entries = 0;
	//double min =currentHistogramState.at(n)

	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		if( (ingredients.getVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
		{
			if(currentHistogramState.at(n) != 0.0)
			{
				mean += currentHistogramState.at(n);
				entries++;
			}
		}
	}

	mean = mean/(1.0*entries);
	
	// this can only happen if the energy window is to restrictive - auto-adjust and broaden the window
	/*std::vector<double> currentTotalHistogramState=ingredients.getTotalVisitsEnergyStates().getVectorValues();//histogram.getVectorValues();
	int entriesTotal = 0;
	for(size_t n=0;n<currentTotalHistogramState.size();n++)
		{
			if( (ingredients.getTotalVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getTotalVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
			{
				if(currentTotalHistogramState.at(n) != 0.0)
				{
					entriesTotal++;
				}
			}
		}
	*/

	/*if(entriesTotal <= 10)
	{
		// broaden the window in 5% in every direction
		double bordershift = (ingredients.getMaxWin()-ingredients.getMinWin())/20.0;

		if(bordershift < ingredients.getVisitsEnergyStates().getBinwidth())
			bordershift = ingredients.getVisitsEnergyStates().getBinwidth();

		ingredients.setMinMaxWin(ingredients.getMinWin()-bordershift, ingredients.getMaxWin()+bordershift);

		//check new boundary
		if(ingredients.getMaxWin() >= ingredients.getVisitsEnergyStates().getMaxCoordinate())
			ingredients.setMinMaxWin(ingredients.getMinWin(), ingredients.getVisitsEnergyStates().getMaxCoordinate());

		if(ingredients.getMinWin() <= ingredients.getVisitsEnergyStates().getMinCoordinate())
			ingredients.setMinMaxWin(ingredients.getVisitsEnergyStates().getMinCoordinate(), ingredients.getMaxWin());

		std::cout << "auto-adjusting window to interval [" << ingredients.getMinWin() << " ; " << ingredients.getMaxWin() << "]" << std::endl;

		return false;
	}
	*/

	bool isFlat = true;

	int entriesTotal = 0;

	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		if( (ingredients.getVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
		{
			if( currentHistogramState.at(n) != 0.0)
			{
				entriesTotal++;
				//	if(std::abs((currentHistogramState.at(n)-mean)/mean) > 0.33)
				if(currentHistogramState.at(n)/mean < flatness)
				{
					isFlat = false;

					return false;
				}
			}
		}
	}

	// to care that we have at least 10 entries avoid stucked conformations
	if(entriesTotal <= 5)
		return false;

	// still here: histogramm converged
	iterationconverged = true;

	iterationfirstconverged = true;

	return isFlat;

	/*
	bool distributionConverged=distributionCheck();
	bool flatnessConverged=flatnessCheck();
	
	bool converged= distributionConverged && flatnessConverged;
	
	std::cout<<"convergence checks: "<<converged<<std::endl;
		
	return converged;
	*/
	
}

//checks if the mean and variance of the reaction coordinate have changed more
//than percentage threshold. also records these two coordinates in the
//histogram convergence vectors for output in cleanup
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::distributionCheck()
{
	//calculate mean and variance
	double mean=0.0;
	double variance=0.0;
	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValuesNormalized();//histogram.getVectorValuesNormalized();
	
	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		mean+=currentHistogramState.at(n)*double(n);
		variance+=currentHistogramState.at(n)*double(n)*double(n);
	}
	variance-=mean*mean;
	
	//calculate the relative change with respect to the last check
	double meanValueChange=std::fabs((mean-oldMean)/mean);
	double varianceChange=std::fabs((variance-oldVariance)/variance);
	
	//calculate the relative change with respect to the second to last check
	double meanValueChangeSecond=std::fabs((mean-oldOldMean)/mean);
	double varianceChangeSecond=std::fabs((variance-oldOldVariance)/variance);
	
	if(meanValueChange == std::numeric_limits<double>::infinity()
		||meanValueChange == -std::numeric_limits<double>::infinity()
		||meanValueChange == std::numeric_limits<double>::signaling_NaN()
		||varianceChange == std::numeric_limits<double>::infinity()
		||varianceChange == -std::numeric_limits<double>::infinity()
		||varianceChange == std::numeric_limits<double>::signaling_NaN()
		||meanValueChangeSecond == std::numeric_limits<double>::infinity()
		||meanValueChangeSecond == -std::numeric_limits<double>::infinity()
		||meanValueChangeSecond == std::numeric_limits<double>::signaling_NaN()
		||varianceChangeSecond == std::numeric_limits<double>::infinity()
		||varianceChangeSecond == -std::numeric_limits<double>::infinity()
		||varianceChangeSecond == std::numeric_limits<double>::signaling_NaN())
	{
		std::cerr<<"WARNING: UpdaterAdaptiveWangLandauSamplingNextNeighbor finds NaN or inf in convergence check in mcs "
		<<ingredients.getMolecules().getAge()<<"\n";
		return false;
	}
	else if(meanValueChange>convergenceThreshold 
		||varianceChange>convergenceThreshold
		||meanValueChangeSecond>convergenceThreshold
		||varianceChangeSecond>convergenceThreshold){
		std::cout<<"distribution not converged!...\n"
		<<"mean change "<<meanValueChange<<" previous mean change "<<meanValueChangeSecond
		<<"\nvariance change "<<varianceChange<<" previous variance change "<<varianceChangeSecond
		<<"\nthreshold "<<convergenceThreshold<<std::endl;
	        oldOldMean=oldMean;
	        oldOldVariance=oldVariance;
		oldMean=mean;
		oldVariance=variance;
	        varianceSeries.push_back(variance);
	        meanSeries.push_back(mean);
		return false;
	}
	else{
		std::cout<<"distribution converged!...\n"
		<<"mean change "<<meanValueChange<<" previous mean change "<<meanValueChangeSecond
		<<"\nvariance change "<<varianceChange<<" previous variance change "<<varianceChangeSecond
		<<"\nthreshold "<<convergenceThreshold<<std::endl;
		
		oldOldMean=oldMean;
	        oldOldVariance=oldVariance;
		oldMean=mean;
		oldVariance=variance;
	        varianceSeries.push_back(variance);
	        meanSeries.push_back(mean);
		
		return true;
	}
	
}

//checks if the flatness of the histogram, measured by square deviation from 
//mean reached some limiting value (fullConvergenceThreshold), and if it has
//changed more than the percentage threshold (convergenceThreshold) in the
//last two checks. 
template<class IngredientsType, class MoveType>
bool UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::flatnessCheck()
{
	double currentFlatness=histogramFlattness();
	flatnessSeries.push_back(currentFlatness);
	
	
	double flatnessChange=std::sqrt(std::pow((currentFlatness-oldFlatnesValue)/currentFlatness,2));
	double flatnessChangeSecond=std::sqrt(std::pow((currentFlatness-oldOldFlatnessValue)/currentFlatness,2));
	
	std::cout<<"histogram flatness: "<<currentFlatness<<std::endl;
	std::cout<<"histogram flatness change: "<<flatnessChange<<" old flatness change: "<<flatnessChangeSecond<<std::endl;
	
	oldOldFlatnessValue=oldFlatnesValue;
	oldFlatnesValue=currentFlatness;
	
	if(currentFlatness==0)
	{
		std::cout<<"WARNING...found NaN or Inf in flatness check\n";
	}
		
		if(flatnessChange<=convergenceThreshold && flatnessChangeSecond<=convergenceThreshold && currentFlatness<=fullFlatnessThreshold){
			simulationConverged=true;
		}
// 		else{
// 			std::vector<double> currentHistogram=histogram.getVectorValues();
// 			bool binsHaveReachedMax=true;
// 			for(size_t n=0;n<currentHistogram.size();n++){
// 				binsHaveReachedMax &= (currentHistogram[n]>maxBinsThreshold);
// 			}
// 			
// 			if (binsHaveReachedMax) simulationConverged=true;
// 		}
		
	
	
	
	if(flatnessChange<=convergenceThreshold && flatnessChangeSecond<=convergenceThreshold){
		
		std::cout<<"flatness converged...\n";
		return true;
	}
	else return false;
}

//calculate the measure for the absolute histogram flatness
template<class IngredientsType, class MoveType>
double UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::histogramFlattness()const
{
	
	std::vector<double> currentHistogramState=ingredients.getVisitsEnergyStates().getVectorValuesNormalized();//.histogram.getVectorValuesNormalized();
	
	//calculate mean
	double mean=1.0/double(currentHistogramState.size());
	double squareDeviation=0.0;
	
	for(size_t n=0;n<currentHistogramState.size();n++)
	{
		squareDeviation+=(currentHistogramState.at(n)-mean)*(currentHistogramState.at(n)-mean)*double(currentHistogramState.size());
	}
	
	squareDeviation=std::sqrt(squareDeviation);
	
	
	return squareDeviation;
	
}

template<class IngredientsType, class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::updateModificationFactor()
{
	double newModificationFactor = std::pow(ingredients.getModificationFactor(ingredients), 0.5);
	//double newModificationFactor = std::exp(nsteps*1.0/ingredients.modifyMolecules().getAge());
	ingredients.setModificationFactor(ingredients, newModificationFactor);

	std::cout << "new modification factor :" << std::setprecision(15) << newModificationFactor << std::endl;
}


//write current state to bfm file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::writeBFM_File(std::string filename)
{
	AnalyzerWriteBfmFile<IngredientsType> ABFM(filename, ingredients);
	ABFM.initialize();
	ABFM.execute();
	ABFM.cleanup();
}

//write current histogram to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpHistogram(std::string filename)
{

	std::ofstream file(filename.c_str());	

	std::vector<double> currentHistogram=ingredients.getVisitsEnergyStates().getVectorValues();//.histogram.getVectorValues();
	std::vector<double> bins=ingredients.getVisitsEnergyStates().getVectorBins();//histogram.getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	//file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	//file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor(ingredients) << (ingredients.using1tMethod() ? " with 1/t sampling" : "") << std::endl;
	file << "# C0-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterCZeroBLENDER() << std::endl;
	file << "# 1N-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterOneOverNBLENDER() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getVisitsEnergyStates().getMinCoordinate() << " ; " << ingredients.getVisitsEnergyStates().getMaxCoordinate() << " ; " << ingredients.getVisitsEnergyStates().getNBins() << " ]" << std::endl;
	file << "# used histogram: [" << ingredients.getMinWin() << " ; " << ingredients.getMaxWin() << " ]" << std::endl;
	file << "# " << std::endl;


	for(size_t n=0;n<currentHistogram.size();n++){
		if( (ingredients.getVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
			if(currentHistogram[n] != 0)
				file<< std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) << currentHistogram[n]<<"\n";
	}
	file.close();
	
}

//write current histogram to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpTotalHistogram(std::string filename)
{

	std::ofstream file(filename.c_str());

	std::vector<double> currentHistogram=ingredients.getTotalVisitsEnergyStates().getVectorValues();//.getVectorValuesNormalized();//.histogram.getVectorValues();
	std::vector<double> bins=ingredients.getTotalVisitsEnergyStates().getVectorBins();//histogram.getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	//file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	//file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor(ingredients) << (ingredients.using1tMethod() ? " with 1/t sampling" : "") << std::endl;
	file << "# C0-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterCZeroBLENDER() << std::endl;
	file << "# 1N-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterOneOverNBLENDER() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getTotalVisitsEnergyStates().getMinCoordinate() << " ; " << ingredients.getTotalVisitsEnergyStates().getMaxCoordinate() << " ; " << ingredients.getTotalVisitsEnergyStates().getNBins() << " ]" << std::endl;
	file << "# used histogram: [" << ingredients.getMinWin() << " ; " << ingredients.getMaxWin() << " ]" << std::endl;
	file << "# " << std::endl;

	for(size_t n=0;n<currentHistogram.size();n++){
		if( (ingredients.getTotalVisitsEnergyStates().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getTotalVisitsEnergyStates().getCenterOfBin(n) <= ingredients.getMaxWin()) )
			if(currentHistogram[n] != 0)
				file << std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) << currentHistogram[n]<<"\n";
	}
	file.close();

}


//write current bias potential to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpHGLnDOS(std::string prefix, double _min, double _max)
{
	std::stringstream filename;
	filename<< prefix;// << "_HGLnDOS" << "_mcs"<<ingredients.getMolecules().getAge()<<".dat";
	std::ofstream file(filename.str().c_str());

	std::vector<StatisticMoment> currentHGLnDOS=ingredients.getHGLnDOS().getVectorValues();
	std::vector<double> bins=ingredients.getHGLnDOS().getVectorBins();

	file << "# " << ingredients.getName() << std::endl;
	file << "# MCS: " << std::setprecision(15) << ingredients.getMolecules().getAge() << std::endl;
	//file << "# f0: " << std::setprecision(15) << initialModificationFactor << std::endl;
	//file << "# f: " << std::setprecision(15) << ingredients.getModificationFactor(ingredients) << (ingredients.using1tMethod() ? " with 1/t sampling" : "") << std::endl;
	file << "# C0-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterCZeroBLENDER() << std::endl;
	file << "# 1N-Blender: " << std::setprecision(15) << ingredients.getComputationalParameterOneOverNBLENDER() << std::endl;
	file << "# Iteration" << iteration << std::endl;
	file << "# histogram: [" << ingredients.getHGLnDOS().getMinCoordinate() << " ; " << ingredients.getHGLnDOS().getMaxCoordinate() << " ; " << ingredients.getHGLnDOS().getNBins() << " ]" << std::endl;
	file << "# used histogram: [" << ingredients.getMinWin() << " ; " << ingredients.getMaxWin() << " ]" << std::endl;
	file << "# vs histogram: [" << _min << " ; " << _max << " ]" << std::endl;
	file << "# " << std::endl;

	for(size_t n=0;n<currentHGLnDOS.size();n++){
		//if( (ingredients.getHGLnDOS().getCenterOfBin(n) >= ingredients.getMinWin()) && (ingredients.getHGLnDOS().getCenterOfBin(n) <= ingredients.getMaxWin()) )
		if( (ingredients.getHGLnDOS().getCenterOfBin(n) >= _min) && (ingredients.getHGLnDOS().getCenterOfBin(n) <= _max) )
			if(currentHGLnDOS[n].ReturnN() != 0)
				file << std::setprecision(15) << bins[n] << "\t" << std::setprecision(15) <<currentHGLnDOS[n].ReturnM1()<<"\n";
	}
	file.close();

}

//write convergence history to file
template<class IngredientsType,class MoveType>
void UpdaterAdaptiveWangLandauSamplingNextNeighbor<IngredientsType,MoveType>::dumpConvergenceProgress()
{
	std::stringstream filename;
	filename<<"convergence.dat";
	std::ofstream file;
	file.open(filename.str().c_str(),std::ofstream::out);
       
	
	for(size_t n=0;n<meanSeries.size();n++){
		file << std::setprecision(15) << meanSeries.at(n)<<" "<<varianceSeries.at(n)<<" "<<flatnessSeries.at(n)<<"\n";
	}
	
	file.close();
	
}



#endif


// Author: Jeremy Pike, 16/10/15
// Contains modified script from the Icy Colocalization Simulator plugin (Thibault Lagache)
// See chapter 1 of the PhD thesis of N. Chenouard, Telecom Paris & Institut Pasteur	

/* This script simply uses the Colocalization Simulator 
 plugin to batch produce data with varying levels of colocalization and noise.

*/

importClass(Packages.icy.file.FileUtil)
importClass(Packages.plugins.adufour.blocks.tools.input.File)
importClass(Packages.plugins.tprovoost.scripteditor.uitools.filedialogs.FileDialog)
importClass(Packages.icy.file.Saver)
importClass(Packages.icy.type.DataType)
importClass(Packages.icy.sequence.SequenceUtil)
importClass(Packages.icy.sequence.Sequence)
importClass(Packages.plugins.adufour.ezplug.EzVarDouble)
importClass(Packages.plugins.adufour.ezplug.EzVarBoolean)
importClass(Packages.plugins.adufour.ezplug.EzVarInteger)
importClass(Packages.plugins.lagache.colocSimulator.generator3d)

//////////////////////////////////
////// User set paramters ////////
/////////////////////////////////

// Colocalization levels to simulate
colocLevels = [0, 0.2, 0.4, 0.6, 0.8, 1]

// Number of stacks to generate for each condition
numberToGen = 50
 
// Noise level used for the Poisson noise and Mean Gaussian noise parameters
noiseLevels = [0, 2, 4]

//Standard deviation of Gaussian noise
std = 1

// Minimum and maximum spot intensities
iMin=20
iMax=20

// Number of particles for channel 1 and channel 2
numParticles1 = 100
numParticles2 = 100

// Dimensions of generated stacks
width = 256
height = 256
numSlices = 50

// Mean distance and standard deviation of colocalization
meanDist = 0
stdDist = 0

//////////////////////////
//////// Script //////////
/////////////////////////


// Great Ez Variables for plugin function
seq_width = new EzVarInteger("Sequence width", width, 1, 100000, 1)
seq_height = new EzVarInteger("Sequence height", height, 1, 100000, 1)
frames = new EzVarInteger("Number of frames", numSlices, 1, 100000, 1)
seq_length = new EzVarInteger("Sequence length",1, 0, 1000, 1)
points = new EzVarBoolean("points", false)
stdGaussian = new EzVarDouble("Std Gaussian noise", 1, 0, 100, 1)
meanPercentage = new EzVarDouble("Mean distance of colocalization", meanDist, 0, 10, 0.1)
stdPercentage = new EzVarDouble("Std distance of colocalization", stdDist, 0, 3, 0.01)

// Open dialog to choose an output folder
folder = FileDialog.saveFolder();

// Loop through all noise levels
for (noiseInd = 0; noiseInd <= noiseLevels.length; noiseInd++) { 

	// Generate Ez variables
	meanGaussian = new EzVarDouble("Mean Gaussian noise", noiseLevels[noiseInd], 0, 100, 1)
	poissonNoise = new EzVarDouble("Poisson noise", noiseLevels[noiseInd], 0, 100, 1)

	//Loop though colocalization levels
	for (colocInd = 0; colocInd < colocLevels.length; colocInd++) {
		colocPercentage = colocLevels[colocInd]
		
		// Genereate the specified number of stacks
		for (i = 0; i < numberToGen; i++) {
			
			// Calculate the number of colocalized and randomly distrubted spots for channel 2
			numParticles2_coloc = Math.min(colocPercentage * numParticles2, numParticles1)
			numParticles2_random = Math.max(numParticles2 - numParticles2_coloc, 0)
	
			
			// New sequences for generated data	
			genSeq1 = new Sequence()
			genSeq2 = new Sequence()

			// Generate data using Colocalization Simulator plugin
			generator3d.main(genSeq1, genSeq2, iMin, iMax,numParticles1, numParticles2, numParticles2_coloc, numParticles2_random, seq_width, seq_height, seq_length, frames, meanPercentage, stdPercentage, meanGaussian, stdGaussian, poissonNoise, points)	

			//Convert data to 8bit
			genSeq1 = SequenceUtil.convertToType(genSeq1, DataType.UBYTE, true)
			genSeq2 = SequenceUtil.convertToType(genSeq2, DataType.UBYTE, true)
				

			// Generate sensible file names and save to specified directory
			file1 = FileUtil.createFile(folder.getAbsolutePath() + "//" + "C1_ColocPercentage_" + colocPercentage + "_NoiseLevel_" + noiseLevels[noiseInd] + "_" + i + ".tif")
			file2 = FileUtil.createFile(folder.getAbsolutePath() + "//" + "C2_ColocPercentage_" + colocPercentage + "_NoiseLevel_" + noiseLevels[noiseInd] + "_" + i + ".tif")
			Saver.save(genSeq1, file1, false, true)
			Saver.save(genSeq2, file2, false, true)
			
	
	
		}
	}
}


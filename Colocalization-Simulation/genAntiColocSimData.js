
// Author: Jeremy Pike, 16/10/15
// Contains modified script from the Icy Colocalization Simulator plugin (Thibault Lagache)
// See chapter 1 of the PhD thesis of N. Chenouard, Telecom Paris & Institut Pasteur	

/* This script uses a modified version of the Colocalization Simulator 
 plugin called Anti Colocalization Simulator to batch produce data with varying 
 levels of anti-colocalization and noise.

*/


importClass(Packages.plugins.jpike.anticolocalizationsimulator.generatorAnti3D)
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

// Anti-colocalization levels to simulate
antiColocLevels = [0.2, 0.4, 0.6, 0.8, 1]

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
seqLength = 1

// Minimum distance for anti-colocalized particles (pixels)
distanceBound = 20

//////////////////////////
//////// Script //////////
/////////////////////////

// Open dialog to choose an output folder
folder = FileDialog.saveFolder();


// Loop through all noise levels
for (noiseInd = 0; noiseInd <= noiseLevels.length; noiseInd++) { 
	

	//Loop though colocalization levels
	for (colocInd = 0; colocInd < antiColocLevels.length; colocInd++) {
		antiColocPercentage = antiColocLevels[colocInd]

		// Genereate the specified number of stacks
		for (i = 0; i < numberToGen; i++) {
			
			// Calculate the number of anti-colocalized and randomly distrubted spots for channel 2
			numParticles2_antiColoc = Math.min(antiColocPercentage * numParticles2, numParticles1)
			numParticles2_random = Math.max(numParticles2 - numParticles2_antiColoc, 0)
	
			// New sequences for generated data	
			genSeq1 = new Sequence()
			genSeq2 = new Sequence()

			// Generate data using Anti-Colocalization Simualator plugin
			generatorAnti3D.main(genSeq1, genSeq2, iMin, iMax, numParticles1, numParticles2, numParticles2_antiColoc, numParticles2_random, seq_width, seq_height, seq_length, frames,  distanceBound, noiseLevels[noiseInd], stdGaussian, noiseLevels[noiseInd])

			//Convert data to 8bit
			genSeq1 = SequenceUtil.convertToType(genSeq1, DataType.UBYTE, true)
			genSeq2 = SequenceUtil.convertToType(genSeq2, DataType.UBYTE, true)
				
		
			// Generate sensible file names and save to specified directory
			file1 = FileUtil.createFile(folder.getAbsolutePath() + "//" + "C1_AntiColocPercentage_" + antiColocPercentage + "_NoiseLevel_" + noiseLevels[noiseInd] + "_" + i + ".tif")
			file2 = FileUtil.createFile(folder.getAbsolutePath() + "//" + "C2_AntiColocPercentage_" + antiColocPercentage + "_NoiseLevel_" + noiseLevels[noiseInd] + "_" + i + ".tif")
			Saver.save(genSeq1, file1, false, true)
			Saver.save(genSeq2, file2, false, true)
			
	
	
		}
	}
}


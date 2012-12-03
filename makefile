COMPILER = openmpic++

OBJFUNC = costruct.o ReadInput.o distruct.o computeCoordinate.o initialiseSurface.o LGConfig.o relabel.o leakSearch.o equilibrium.o sign.o initialise.o computeMomenta.o ARolling.o densitiesAtSolidBoundaries.o computeFreeEnergy.o invComputeCoordinate.o collision.o writeDensityFile.o  saveFiles.o propagation.o applyBoundaryConditions.o LBAlgorithm.o makematrix.o commReadData.o exchangeMask.o exchangePhi.o generateGlobalMask.o generatePhiGlobal.o exchangeDensities.o String.o writeZPlanXVelocityFile.o duplicateDrop1.o duplicateDrop2.o generateNGlobal.o generateGlobal.o duplicateArray.o duplicateArray_int.o computeContactArea.o writeInfoFile.o writeJumpFile.o initialiseSurfaceLater.o computeDissipation.o writeDissipationFile.o exchangeVelocities.o exchangeChemPot.o LGConfigRev.o exchangeDensities_ffgg.o Dimensionlessnum.o writehandlefile.o

main: $(OBJFUNC) ./functions/main.cpp 
	$(COMPILER) $(OBJFUNC) -o wetS ./functions/main.cpp
	#rm *.o

costruct.o: ./functions/costruct.cpp
	$(COMPILER) -c ./functions/costruct.cpp

computeCoordinate.o: ./functions/computeCoordinate.cpp
	$(COMPILER) -c ./functions/computeCoordinate.cpp

ReadInput.o: ./functions/ReadInput.cpp
	$(COMPILER) -c ./functions/ReadInput.cpp

initialiseSurface.o: ./functions/initialiseSurface.cpp
	$(COMPILER) -c ./functions/initialiseSurface.cpp

LGConfig.o: ./functions/LGConfig.cpp
	$(COMPILER) -c ./functions/LGConfig.cpp

relabel.o: ./functions/relabel.cpp
	$(COMPILER) -c ./functions/relabel.cpp

leakSearch.o: ./functions/leakSearch.cpp
	$(COMPILER) -c ./functions/leakSearch.cpp

equilibrium.o: ./functions/equilibrium.cpp
	$(COMPILER) -c ./functions/equilibrium.cpp

sign.o: ./functions/sign.cpp
	$(COMPILER) -c ./functions/sign.cpp

initialise.o:  ./functions/initialise.cpp
	$(COMPILER) -c ./functions/initialise.cpp

computeMomenta.o: ./functions/computeMomenta.cpp
	$(COMPILER) -c ./functions/computeMomenta.cpp

ARolling.o: ./functions/ARolling.cpp
	$(COMPILER) -c ./functions/ARolling.cpp

densitiesAtSolidBoundaries.o: ./functions/densitiesAtSolidBoundaries.cpp
	$(COMPILER) -c ./functions/densitiesAtSolidBoundaries.cpp

computeFreeEnergy.o: ./functions/computeFreeEnergy.cpp
	$(COMPILER) -c ./functions/computeFreeEnergy.cpp

invComputeCoordinate.o: ./functions/invComputeCoordinate.cpp
	$(COMPILER) -c ./functions/invComputeCoordinate.cpp

collision.o: ./functions/collision.cpp
	$(COMPILER) -c ./functions/collision.cpp

writeDensityFile.o: ./functions/writeDensityFile.cpp
	$(COMPILER) -c ./functions/writeDensityFile.cpp

saveFiles.o: ./functions/saveFiles.cpp
	$(COMPILER) -c ./functions/saveFiles.cpp

propagation.o: ./functions/propagation.cpp
	$(COMPILER) -c ./functions/propagation.cpp

applyBoundaryConditions.o: ./functions/applyBoundaryConditions.cpp
	$(COMPILER) -c ./functions/applyBoundaryConditions.cpp

LBAlgorithm.o: ./functions/LBAlgorithm.cpp
	$(COMPILER) -c ./functions/LBAlgorithm.cpp

makematrix.o: ./functions/makematrix.cpp
	$(COMPILER) -c ./functions/makematrix.cpp

commReadData.o: ./functions/commReadData.cpp
	$(COMPILER) -c ./functions/commReadData.cpp

exchangeMask.o: ./functions/exchangeMask.cpp
	$(COMPILER) -c ./functions/exchangeMask.cpp

generateGlobalMask.o: ./functions/generateGlobalMask.cpp
	$(COMPILER) -c ./functions/generateGlobalMask.cpp

exchangePhi.o: ./functions/exchangePhi.cpp
	$(COMPILER) -c ./functions/exchangePhi.cpp

generatePhiGlobal.o: ./functions/generatePhiGlobal.cpp
	$(COMPILER) -c ./functions/generatePhiGlobal.cpp

exchangeDensities.o: ./functions/exchangeDensities.cpp
	$(COMPILER) -c ./functions/exchangeDensities.cpp

distruct.o: ./functions/distruct.cpp
	$(COMPILER) -c ./functions/distruct.cpp

String.o: ./functions/String.cc
	$(COMPILER) -c ./functions/String.cc

writeZPlanXVelocityFile.o: ./functions/writeZPlanXVelocityFile.cpp
	$(COMPILER) -c ./functions/writeZPlanXVelocityFile.cpp

duplicateDrop1.o: ./functions/duplicateDrop1.cpp
	$(COMPILER) -c ./functions/duplicateDrop1.cpp

duplicateDrop2.o: ./functions/duplicateDrop2.cpp
	$(COMPILER) -c ./functions/duplicateDrop2.cpp

generateNGlobal.o: ./functions/generateNGlobal.cpp
	$(COMPILER) -c ./functions/generateNGlobal.cpp

generateGlobal.o: ./functions/generateGlobal.cpp
	$(COMPILER) -c ./functions/generateGlobal.cpp

duplicateArray.o: ./functions/duplicateArray.cpp
	$(COMPILER) -c ./functions/duplicateArray.cpp

duplicateArray_int.o: ./functions/duplicateArray_int.cpp
	$(COMPILER) -c ./functions/duplicateArray_int.cpp

computeContactArea.o: ./functions/computeContactArea.cpp
	$(COMPILER) -c ./functions/computeContactArea.cpp

writeInfoFile.o: ./functions/writeInfoFile.cpp
	$(COMPILER) -c ./functions/writeInfoFile.cpp

writeJumpFile.o: ./functions/writeJumpFile.cpp
	$(COMPILER) -c ./functions/writeJumpFile.cpp

initialiseSurfaceLater.o: ./functions/initialiseSurfaceLater.cpp
	$(COMPILER) -c ./functions/initialiseSurfaceLater.cpp

computeDissipation.o: ./functions/computeDissipation.cpp
	$(COMPILER) -c ./functions/computeDissipation.cpp

writeDissipationFile.o: ./functions/writeDissipationFile.cpp
	$(COMPILER) -c ./functions/writeDissipationFile.cpp

exchangeVelocities.o: ./functions/exchangeVelocities.cpp
	$(COMPILER) -c ./functions/exchangeVelocities.cpp

exchangeChemPot.o: ./functions/exchangeChemPot.cpp
	$(COMPILER) -c ./functions/exchangeChemPot.cpp

LGConfigRev.o: ./functions/LGConfigRev.cpp
	$(COMPILER) -c ./functions/LGConfigRev.cpp

exchangeDensities_ffgg.o: ./functions/exchangeDensities_ffgg.cpp
	$(COMPILER) -c ./functions/exchangeDensities_ffgg.cpp

Dimensionlessnum.o: ./functions/Dimensionlessnum.cpp
	$(COMPILER) -c ./functions/Dimensionlessnum.cpp

writehandlefile.o: ./functions/writehandlefile.cpp
	$(COMPILER) -c ./functions/writehandlefile.cpp

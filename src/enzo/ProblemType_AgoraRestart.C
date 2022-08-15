/***********************************************************************
/
/  Agora isolated galaxy restart
/
/  written by: Nathan Goldbaum
/  date:       March, 2013
/
/  PURPOSE:
/  https://sites.google.com/site/projectagoraworkspace/metagroup1/group2
/  https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAHZyjrfTz88aG12H0Q_Rqla
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "ProblemType.h"
#include "EventHooks.h"
#include "phys_constants.h"


#define VCIRC_TABLE_LENGTH 1000000

void mt_init(unsigned_int seed);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
inline int nlines(const char* fname);

int nlines(const char* fname) {

  FILE* fptr = fopen(fname, "r");
  int ch, n = 0;

  do
  {
    ch = fgetc(fptr);
    if(ch == '\n')
      n++;
  } while (ch != EOF);

  fclose(fptr);
  if (debug) fprintf(stderr,"Read %"ISYM" lines \n", n);
  return n;
}

class ProblemType_AgoraRestart;

class AgoraRestartGrid : private grid
{
  friend class ProblemType_AgoraRestart;
};

class ProblemType_AgoraRestart : public EnzoProblemType
{
private:
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT CenterPosition[MAX_DIMENSION];
  float Bfield[MAX_DIMENSION];
  FLOAT ScaleLength;
  FLOAT ScaleHeight;
  float DiskMass;
  float GasFraction;
  float DiskTemperature;
  float DiskMetallicity;
  float HaloMass;
  float HaloTemperature;
  float HaloMetallicity;
  FLOAT VCircRadius[VCIRC_TABLE_LENGTH];
  float VCircVelocity[VCIRC_TABLE_LENGTH];
  //added by Shin
  float DiskTheta;
  float DiskPhi;
  float DiskTheta2;
  float DiskPhi2;
  FLOAT CenterPosition2[MAX_DIMENSION];
  FLOAT VCircRadius2[VCIRC_TABLE_LENGTH];
  float VCircVelocity2[VCIRC_TABLE_LENGTH];
  FLOAT GalaxyVelocity[MAX_DIMENSION];
  FLOAT GalaxyVelocity2[MAX_DIMENSION];
  //added by Joohyun Lee
  FLOAT ScaleLength2;
  FLOAT ScaleHeight2;
  float DiskMass2;
  float GasFraction2;

  int RefineAtStart;

public:
  ProblemType_AgoraRestart() : EnzoProblemType()
  {
    if (MyProcessorNumber == 0)
      std::cout << "Creating problem type Agora Restart" << std::endl;
  }

  ~ProblemType_AgoraRestart() {}

  virtual int InitializeFromRestart(
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    return SUCCESS;
  }

  virtual int InitializeSimulation(
    FILE *fptr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    if(debug)
    {
      printf("Entering AgoraRestartInitialize\n");
      fflush(stdout);
    }

    char *DensName = "Density";
    char *TEName   = "TotalEnergy";
    char *GEName   = "GasEnergy";
    char *Vel1Name = "x-velocity";
    char *Vel2Name = "y-velocity";
    char *Vel3Name = "z-velocity";
    char *ElectronName = "Electron_Density";
    char *HIName    = "HI_Density";
    char *HIIName   = "HII_Density";
    char *HeIName   = "HeI_Density";
    char *HeIIName  = "HeII_Density";
    char *HeIIIName = "HeIII_Density";
    char *HMName    = "HM_Density";
    char *H2IName   = "H2I_Density";
    char *H2IIName  = "H2II_Density";
    char *DIName    = "DI_Density";
    char *DIIName   = "DII_Density";
    char *HDIName   = "HDI_Density";
    char *MetalName = "Metal_Density";
    char *MetalSNIaName = "MetalSNIa_Density";
    char *MetalSNIIName = "MetalSNII_Density";
    char *BxName = "Bx";
    char *ByName = "By";
    char *BzName = "Bz";
    char *PhiName = "Phi";
    char *MachNum = "Mach";

    /* local declarations */

    char line[MAX_LINE_LENGTH];
    int  i, ret, level;

    /* make sure it is 3D */

    if (MetaData.TopGridRank != 3)
    {
      printf("Cannot do AgoraRestart in %"ISYM" dimension(s)\n",
	     MetaData.TopGridRank);
      ENZO_FAIL("Agora Restart simulations must be 3D!");
    }

    for (i=0; i < MAX_DIMENSION; i++)
    {
      this->CenterPosition[i] = 0.5;
      this->CenterPosition2[i] = 0.5;//By Shin
      this->Bfield[i] = 0.;
    }

    // These come from Oscar's sample output.  The units are:
    // Velocity: km/s
    // Mass: 10^9 Msun
    // Length: kpc
    // Temperature: K
    this->ScaleLength         = .0343218;
    this->ScaleLength2        = .0343218;
    this->ScaleHeight         = .00343218;
    this->ScaleHeight2        = .00343218;
    this->DiskMass            = 42.9661;
    this->DiskMass2           = 42.9661;
    this->GasFraction         = 0.2;
    this->GasFraction2        = 0.2;
    this->DiskTemperature     = 1e4;
    this->DiskMetallicity     = 0.0;
    this->HaloMass            = 0.10000;
    this->HaloTemperature     = this->DiskTemperature;
    this->HaloMetallicity     = 0.0;
    this->RefineAtStart       = TRUE;

    // set this from global data (kind of a hack)
    TestProblemData.MultiSpecies = MultiSpecies;

    /* read input from file */
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret = 0;
      ret += sscanf(line, "AgoraRestartCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		    CenterPosition, CenterPosition+1, CenterPosition+2);
      ret += sscanf(line, "AgoraRestartCenterPosition2 = %"PSYM" %"PSYM" %"PSYM,
                    CenterPosition2, CenterPosition2+1, CenterPosition2+2);

      ret += sscanf(line, "AgoraRestartDiskTheta = %"PSYM,
		    &DiskTheta);
      ret += sscanf(line, "AgoraRestartDiskPhi = %"PSYM,
                    &DiskPhi);
      ret += sscanf(line, "AgoraRestartDiskTheta2 = %"PSYM,
                    &DiskTheta2);
      ret += sscanf(line, "AgoraRestartDiskPhi2 = %"PSYM,
                    &DiskPhi2);
      ret += sscanf(line, "AgoraRestartScaleLength = %"PSYM, &ScaleLength);
      ret += sscanf(line, "AgoraRestartScaleLength2 = %"PSYM, &ScaleLength2);
      ret += sscanf(line, "AgoraRestartScaleHeight = %"PSYM, &ScaleHeight);
      ret += sscanf(line, "AgoraRestartScaleHeight2 = %"PSYM, &ScaleHeight2);
      ret += sscanf(line, "AgoraRestartDiskMass = %"FSYM, &DiskMass);
      ret += sscanf(line, "AgoraRestartDiskMass2 = %"FSYM, &DiskMass2);
      ret += sscanf(line, "AgoraRestartGasFraction = %"FSYM, &GasFraction);
      ret += sscanf(line, "AgoraRestartGasFraction2 = %"FSYM, &GasFraction2);
      ret += sscanf(line, "AgoraRestartDiskTemperature = %"FSYM,
		    &DiskTemperature);
      ret += sscanf(line, "AgoraRestartDiskMetallicity = %"FSYM,
		    &DiskMetallicity);
      ret += sscanf(line, "AgoraRestartHaloMass = %"FSYM, &HaloMass);
      ret += sscanf(line, "AgoraRestartHaloTemperature = %"FSYM,
		    &HaloTemperature);
      ret += sscanf(line, "AgoraRestartHaloMetallicity = %"FSYM,
                    &HaloMetallicity);

      ret += sscanf(line, "AgoraRestartMagneticField = %"FSYM" %"FSYM" %"FSYM,
		    Bfield, Bfield+1, Bfield+2);

      ret += sscanf(line, "AgoraRestartGalaxyVelocity = %"PSYM" %"PSYM" %"PSYM,
                    GalaxyVelocity, GalaxyVelocity+1, GalaxyVelocity+2);
      ret += sscanf(line, "AgoraRestartGalaxyVelocity2 = %"PSYM" %"PSYM" %"PSYM,
                    GalaxyVelocity2, GalaxyVelocity2+1, GalaxyVelocity2+2);

      ret += sscanf(line, "AgoraRestartRefineAtStart = %"ISYM,
		    &RefineAtStart);
      ret += sscanf(line, "AgoraRestartHydrogenFractionByMass = %"FSYM,
		    &TestProblemData.HydrogenFractionByMass);
      ret += sscanf(line, "AgoraRestartHeliumFractionByMass = %"FSYM,
		    &TestProblemData.HeliumFractionByMass);
      ret += sscanf(line, "AgoraRestartMetalFractionByMass = %"FSYM,
		    &TestProblemData.MetalFractionByMass);
      ret += sscanf(line, "AgoraRestartDeuteriumToHydrogenRatio = %"FSYM,
		    &TestProblemData.DeuteriumToHydrogenRatio);
      ret += sscanf(line, "AgoraRestartInitialHIFraction  = %"FSYM,
		    &TestProblemData.HI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHIIFraction  = %"FSYM,
		    &TestProblemData.HII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIFraction  = %"FSYM,
		    &TestProblemData.HeI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIIFraction  = %"FSYM,
		    &TestProblemData.HeII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIIIFraction  = %"FSYM,
		    &TestProblemData.HeIII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHMFraction  = %"FSYM,
		    &TestProblemData.HM_Fraction);
      ret += sscanf(line, "AgoraRestartInitialH2IFraction  = %"FSYM,
		    &TestProblemData.H2I_Fraction);
      ret += sscanf(line, "AgoraRestartInitialH2IIFraction  = %"FSYM,
		    &TestProblemData.H2II_Fraction);
      ret += sscanf(line, "AgoraRestartInitialDIFraction  = %"FSYM,
		    &TestProblemData.DI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialDIIFraction  = %"FSYM,
		    &TestProblemData.DII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHDIFraction  = %"FSYM,
		    &TestProblemData.HDI_Fraction);
      ret += sscanf(line, "AgoraRestartUseMetallicityField  = %"ISYM,
		    &TestProblemData.UseMetallicityField);


      if (ret == 0 && strstr(line, "=") &&
	  (strstr(line, "AgoraRestart") || strstr(line, "TestProblem")) &&
	  line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr,
		"*** warning: the following parameter line from AgoraRestart was not interpreted:\n%s\n",
		line);

    } // end input from parameter file


    // Read in circular velocity table

    this->ReadInVcircData();
    this->ReadInVcircData2();


    /* set up top grid */

    float dummy_density = 1.0;
    float dummy_gas_energy = 1.0; // Only used if DualEnergyFormalism is True
    float dummy_total_energy = 1.0;
    float dummy_velocity[3] = {0.0, 0.0, 0.0};
    float dummy_b_field[3] = {1e-20, 1e-20, 1e-20}; // Only set if HydroMethod = mhd_rk

    if (this->InitializeUniformGrid(
	  TopGrid.GridData, dummy_density, dummy_total_energy,
	  dummy_gas_energy, dummy_velocity, dummy_b_field) == FAIL)
    {
      ENZO_FAIL("Error in InitializeUniformGrid");
    }

    this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData);

    this->InitializeParticles(TopGrid.GridData, TopGrid, MetaData);

    /* Convert minimum initial overdensity for refinement to mass
       (unless MinimumMass itself was actually set). */
    //printf("###############\n");
    //printf("1.MinimumOverDensityForRefinement=%e\n",MinimumOverDensityForRefinement[0]);
    //printf("2.TopGridRank=%d\n", MetaData.TopGridRank);


    if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
      for (int dim = 0; dim < MetaData.TopGridRank; dim++)
	{	MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
	  //printf("delta-DomainRightEdge[%d]=%lf\n",dim,DomainRightEdge[dim]-DomainLeftEdge[dim]);
	  //printf("MetaData.TopGridDims[%d]=%d\n",dim,MetaData.TopGridDims[dim]);
	    }
    }
    //printf("3.MinimumMassForRefinement=%e\n",MinimumMassForRefinement[0]);

    /* If requested, refine the grid to the desired level. */

    if (RefineAtStart)
    {
      /* Declare, initialize, and fill out the first level of the LevelArray. */
      LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	LevelArray[level] = NULL;
      AddLevel(LevelArray, &TopGrid, 0);

      /* Add levels to the maximum depth or until no new levels are created,
	 and re-initialize the level after it is created. */
      for (level = 0; level < MaximumRefinementLevel; level++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	  fprintf(stderr, "Error in RebuildHierarchy.\n");
	  return FAIL;
	}
	if (LevelArray[level+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level+1];
	while (Temp != NULL) {
	  if (this->InitializeGrid(Temp->GridData, TopGrid, MetaData) == FAIL)
	  {
	    ENZO_FAIL("Error in AgoraRestart->InitializeGrid");
	  }
	  Temp = Temp->NextGridThisLevel;
	} // end: loop over grids on this level
      } // end: loop over levels
    }

    /* set up field names and units */
    int count = 0;
    DataLabel[count++] = DensName;
    DataLabel[count++] = Vel1Name;
    if(MetaData.TopGridRank > 1)
      DataLabel[count++] = Vel2Name;
    if(MetaData.TopGridRank > 2)
      DataLabel[count++] = Vel3Name;
    DataLabel[count++] = TEName;
    if (DualEnergyFormalism)
      DataLabel[count++] = GEName;

    if (HydroMethod == MHD_RK) {
      DataLabel[count++] = (char*) BxName;
      DataLabel[count++] = (char*) ByName;
      DataLabel[count++] = (char*) BzName;
      DataLabel[count++] = (char*) PhiName;
    }

    if (MultiSpecies)
    {
      DataLabel[count++] = ElectronName;
      DataLabel[count++] = HIName;
      DataLabel[count++] = HIIName;
      DataLabel[count++] = HeIName;
      DataLabel[count++] = HeIIName;
      DataLabel[count++] = HeIIIName;
      if (MultiSpecies > 1)
      {
	DataLabel[count++] = HMName;
	DataLabel[count++] = H2IName;
	DataLabel[count++] = H2IIName;
      }
      if (MultiSpecies > 2)
      {
	DataLabel[count++] = DIName;
	DataLabel[count++] = DIIName;
	DataLabel[count++] = HDIName;
      }
    }
    if (TestProblemData.UseMetallicityField)
      DataLabel[count++] = MetalName;
    if (StarMakerTypeIaSNe)
        DataLabel[count++] = MetalSNIaName;
    if (StarMakerTypeIISNeMetalField)
        DataLabel[count++] = MetalSNIIName;
    for (i = 0; i < count; i++)
      DataUnits[i] = NULL;
    

    if (MyProcessorNumber == ROOT_PROCESSOR)
    {
      fprintf(Outfptr, "AgoraRestartCenterPosition          = %"
	      PSYM" %"PSYM" %"PSYM"\n",
	      CenterPosition[0], CenterPosition[1], CenterPosition[2]);
      fprintf(Outfptr, "AgoraRestartCenterPosition2          = %"
              PSYM" %"PSYM" %"PSYM"\n",
              CenterPosition2[0], CenterPosition2[1], CenterPosition2[2]);
      fprintf(Outfptr, "AgoraRestartGalaxyVelocity               = %"
              PSYM" %"PSYM" %"PSYM"\n",
              GalaxyVelocity[0], GalaxyVelocity[1], GalaxyVelocity[2]);
      fprintf(Outfptr, "AgoraRestartGalaxyVelocity2                = %"
              PSYM" %"PSYM" %"PSYM"\n",
	      GalaxyVelocity2[0], GalaxyVelocity2[1], GalaxyVelocity2[2]);
      fprintf(Outfptr, "AgoraRestartMagneticField           = %"FSYM" %"FSYM" %"FSYM,
		    Bfield[0], Bfield[1], Bfield[2]);
      fprintf(Outfptr, "AgoraRestartScaleLength             = %"PSYM"\n",
	      ScaleLength);
      fprintf(Outfptr, "AgoraRestartScaleLength2             = %"PSYM"\n",
	      ScaleLength);
      fprintf(Outfptr, "AgoraRestartScaleHeight             = %"PSYM"\n",
	      ScaleHeight);
      fprintf(Outfptr, "AgoraRestartScaleHeight2             = %"PSYM"\n",
	      ScaleHeight);
      fprintf(Outfptr, "AgoraRestartDiskMass                = %"FSYM"\n",
	      DiskMass);
      fprintf(Outfptr, "AgoraRestartDiskMass2                = %"FSYM"\n",
	      DiskMass);
      fprintf(Outfptr, "AgoraRestartGasFraction             = %"FSYM"\n",
	      GasFraction);
      fprintf(Outfptr, "AgoraRestartGasFraction2             = %"FSYM"\n",
	      GasFraction);
      fprintf(Outfptr, "AgoraRestartDiskTemperature         = %"FSYM"\n",
	      DiskTemperature);
      fprintf(Outfptr, "AgoraRestartHaloMass                = %"FSYM"\n",
	      HaloMass);
      fprintf(Outfptr, "AgoraRestartHaloTemperature         = %"FSYM"\n",
	      HaloTemperature);
      fprintf(Outfptr, "AgoraRestartRefineAtStart           = %"ISYM"\n",
	      RefineAtStart);
      fprintf(Outfptr, "AgoraRestartHydrogenFractionByMass = %"FSYM"\n",
	      TestProblemData.HydrogenFractionByMass);
      fprintf(Outfptr, "AgoraRestartHeliumFractionByMass = %"FSYM"\n",
	      TestProblemData.HeliumFractionByMass);
      fprintf(Outfptr, "AgoraRestartMetalFractionByMass = %"FSYM"\n",
	      TestProblemData.MetalFractionByMass);
      fprintf(Outfptr, "AgoraRestartInitialHIFraction  = %"FSYM"\n",
	      TestProblemData.HI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHIIFraction  = %"FSYM"\n",
	      TestProblemData.HII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIFraction  = %"FSYM"\n",
	      TestProblemData.HeI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIIFraction  = %"FSYM"\n",
	      TestProblemData.HeII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIIIIFraction  = %"FSYM"\n",
	      TestProblemData.HeIII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHMFraction  = %"FSYM"\n",
	      TestProblemData.HM_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialH2IFraction  = %"FSYM"\n",
	      TestProblemData.H2I_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialH2IIFraction  = %"FSYM"\n",
	      TestProblemData.H2II_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialDIFraction  = %"FSYM"\n",
	      TestProblemData.DI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialDIIFraction  = %"FSYM"\n",
	      TestProblemData.DII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHDIFraction  = %"FSYM"\n",
	      TestProblemData.HDI_Fraction);
      fprintf(Outfptr, "AgoraRestartUseMetallicityField  = %"ISYM"\n",
	      TestProblemData.UseMetallicityField);
    }
 
    return SUCCESS;

  } // InitializeSimulation

  int InitializeGrid(grid *thisgrid_orig, HierarchyEntry &TopGrid,
		     TopGridData &MetaData)
  {

    if(debug)
      printf("Entering AgoraRestart InitializeGrid\n");

    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    /* Identify physical quantities */
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, MetalNum;

    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					     Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    if (TestProblemData.MultiSpecies)
      if (thisgrid->IdentifySpeciesFields(
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
	    HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
	ENZO_FAIL("Error in grid->IdentifySpeciesFields.");

    int MetallicityField = FALSE;
    if ((MetalNum = FindField(
	   Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)
	  ) != -1)
      MetallicityField = TRUE;
    else
      MetalNum = 0;

    int dim, i, j, k, size, index=0;
    float RhoZero, RhoZero2, DiskGasEnergy, DiskDensity, HaloGasEnergy, HaloDensity, BoxVolume, vcirc, vcirc2, mu;
    FLOAT cellwidth, x,y,z,radius,xy_radius,x2,y2,z2,radius2,xy_radius2,theta,theta2,phi,phi2,Rx,Ry,Rz,Rx2,Ry2,Rz2;

    /* Compute size of this grid */
    size = 1;
    for (dim = 0; dim < thisgrid->GridRank; dim++)
      size *= thisgrid->GridDimension[dim];
    cellwidth = thisgrid->CellWidth[0][0];

    /* Compute the size of the box */
    BoxVolume = 1.;
    for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
      BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

    /* Find the mean molecular weight */

    if (TestProblemData.MultiSpecies == FALSE)
      mu = Mu;
    else
    {
      // Atomic hydrogen
      mu = TestProblemData.HydrogenFractionByMass *
	(TestProblemData.HI_Fraction + 2.0*TestProblemData.HII_Fraction);

      // Helium
      mu += TestProblemData.HeliumFractionByMass / 4.0 *
	(TestProblemData.HeI_Fraction + 2.0*TestProblemData.HeII_Fraction +
	 3.0*TestProblemData.HeIII_Fraction);

      // Molecular hydrogen, ignore Deuterium
      if (TestProblemData.MultiSpecies > 1)
	mu += TestProblemData.HydrogenFractionByMass / 2.0 *
	  (TestProblemData.H2I_Fraction + 2.0*TestProblemData.H2II_Fraction);

      // Metals
      if (TestProblemData.UseMetallicityField)
	mu += TestProblemData.MetalFractionByMass / 16.0;

      mu = POW(mu, -1);

    }

    /* Find global physical properties */
    //RhoZero = this->DiskMass * this->GasFraction / (4.*pi) / (POW((this->ScaleLength),2)*(this->ScaleHeight));
    //RhoZero2 = this->DiskMass2 * this->GasFraction2 / (4.*pi) / (POW((this->ScaleLength2),2)*(this->ScaleHeight2));
      
    // Dec 2021 JL isothermal initialization, assuming 5*ScaleLength is the boundary of the gaseous sphere
    // 5 = Rmax/Rc
    RhoZero = this->DiskMass * this->GasFraction / (5.*pi) / (5 - atan(5)) / POW((this->ScaleLength),3);
    RhoZero2 = this->DiskMass2 * this->GasFraction2 / (5.*pi) / (5 - atan(5)) / POW((this->ScaleLength2),3);


    HaloGasEnergy = this->HaloTemperature / mu / (Gamma - 1) /
      TemperatureUnits;

    HaloDensity = this->HaloMass / BoxVolume;

    DiskGasEnergy = this->DiskTemperature / mu / (Gamma - 1) /
      TemperatureUnits;



    /* Loop over the mesh. */

    for (k = 0; k < thisgrid->GridDimension[2]; k++)
    {
      for (j = 0; j < thisgrid->GridDimension[1]; j++)
      {
	for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
	{

	  thisgrid->BaryonField[Vel1Num][index] = 0;
	  thisgrid->BaryonField[Vel2Num][index] = 0;
	  thisgrid->BaryonField[Vel3Num][index] = 0;

	  /* Compute position */

	  x = (thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i]) *
	    LengthUnits;
	  y = (thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j]) *
	    LengthUnits;
	  z = (thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]) *
	    LengthUnits;
	  x2 = (thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i]) *
            LengthUnits;
      	  y2 = (thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j]) *
            LengthUnits;
      	  z2 = (thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]) *
            LengthUnits;


	  x -= this->CenterPosition[0]*LengthUnits;
	  y -= this->CenterPosition[1]*LengthUnits;
	  z -= this->CenterPosition[2]*LengthUnits;

	  x2 -= this->CenterPosition2[0]*LengthUnits;
      	  y2 -= this->CenterPosition2[1]*LengthUnits;
      	  z2 -= this->CenterPosition2[2]*LengthUnits;
      
      // converting into radian unit
	  theta=this->DiskTheta*pi/180.0;
	  theta2=this->DiskTheta2*pi/180.0;
	  phi=this->DiskPhi*pi/180.0; 
	  phi2=this->DiskPhi2*pi/180.0;
	  
      // ZXZ euler angles rotation matrix
	  Rx=x*cos(theta)+(-y*cos(theta)+z*sin(theta))*sin(phi);
	  Ry=x*sin(theta)+(y*cos(theta)-z*sin(theta))*cos(phi);
	  Rz=y*sin(theta)+z*cos(theta);

	  Rx2=x2*cos(theta2)+(-y2*cos(theta2)+z2*sin(theta2))*sin(phi2);
	  Ry2=x2*sin(theta2)+(y2*cos(theta2)-z2*sin(theta2))*cos(phi2);
	  Rz2=y2*sin(theta2)+z2*cos(theta2);

	  radius = sqrt(POW(x, 2) +
			POW(y, 2) +
			POW(z, 2) );

	  xy_radius = sqrt(POW(Rx, 2) +
			   POW(Ry, 2) );

	  radius2 = sqrt(POW(x2, 2) +
            POW(y2, 2) +
            POW(z2, 2) );

      xy_radius2 = sqrt(POW(Rx2, 2) +
               POW(Ry2, 2) );

	  /* Find disk density, halo density and internal energy */

	  DiskDensity = gauss_mass(RhoZero, Rx/LengthUnits, Ry/LengthUnits,
				   Rz/LengthUnits, cellwidth) / POW(cellwidth, 3)
	    +gauss_mass(RhoZero2, Rx2/LengthUnits, Ry2/LengthUnits,
			Rz2/LengthUnits, cellwidth) / POW(cellwidth, 3);

	  float maxlen=5.0;
	  //if ( HaloDensity*HaloTemperature >= DiskDensity*DiskTemperature ) // In halo
	  if (radius>(ScaleLength*LengthUnits*maxlen) && radius2>(ScaleLength*LengthUnits*maxlen))
	  {
	    thisgrid->BaryonField[DensNum][index] = HaloDensity;
        thisgrid->BaryonField[TENum][index] = HaloGasEnergy;
        if (DualEnergyFormalism)
            thisgrid->BaryonField[GENum][index] = HaloGasEnergy;

        thisgrid->BaryonField[Vel1Num][index] = 0;
        thisgrid->BaryonField[Vel2Num][index] = 0;
        thisgrid->BaryonField[Vel3Num][index] = 0;

	    if (TestProblemData.UseMetallicityField)
              thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
                TestProblemData.MetalFractionByMass * HaloMetallicity;
        
        // galaxy halo velocity initialization
        if (radius<=ScaleLength*LengthUnits*maxlen*3) // Galaxy 1
	    {
	      thisgrid->BaryonField[Vel1Num][index] = GalaxyVelocity[0]*1e5/VelocityUnits;
	      thisgrid->BaryonField[Vel2Num][index] = GalaxyVelocity[1]*1e5/VelocityUnits;
	      thisgrid->BaryonField[Vel3Num][index] = GalaxyVelocity[2]*1e5/VelocityUnits;
	    }
        if (radius2<=ScaleLength*LengthUnits*maxlen*3) // Galaxy 2, overlapping region should be considered carefully (not for halo)
	    {
	      thisgrid->BaryonField[Vel1Num][index] = GalaxyVelocity2[0]*1e5/VelocityUnits;
	      thisgrid->BaryonField[Vel2Num][index] = GalaxyVelocity2[1]*1e5/VelocityUnits;
	      thisgrid->BaryonField[Vel3Num][index] = GalaxyVelocity2[2]*1e5/VelocityUnits;
	    }
	  }
      
	  else // Ok, we're in the disk
	  {
	    thisgrid->BaryonField[DensNum][index] = DiskDensity; 
	    thisgrid->BaryonField[TENum][index] = DiskGasEnergy;
        if (DualEnergyFormalism) {
            thisgrid->BaryonField[GENum][index] = DiskGasEnergy;
        }
        
	
        if (radius<=radius2) {
            vcirc = this->InterpolateVcircTable(xy_radius);
            
            // Dec 2021 JL no-rotation isothermal initialization
            //vcirc = 0.0;
		
            thisgrid->BaryonField[Vel1Num][index] = (-vcirc*y/xy_radius/VelocityUnits)*cos(theta)-(vcirc*x/xy_radius/VelocityUnits)*cos(theta)*sin(phi)
                                                    +GalaxyVelocity[0]*1e5/VelocityUnits;
            thisgrid->BaryonField[Vel2Num][index] = (-vcirc*y/xy_radius/VelocityUnits)*sin(theta)+(vcirc*x/xy_radius/VelocityUnits)*cos(theta)*cos(phi)
                                                    +GalaxyVelocity[1]*1e5/VelocityUnits;
            thisgrid->BaryonField[Vel3Num][index] = (vcirc*x/xy_radius/VelocityUnits)*sin(theta)
                                                    +GalaxyVelocity[2]*1e5/VelocityUnits;
		}//Galaxy1                                                                                           
	    
        else {
            vcirc2 = this->InterpolateVcircTable2(xy_radius2);
            
            // Dec 2021 JL isothermal initialization
            vcirc2 = 0.0;
            
            thisgrid->BaryonField[Vel1Num][index] = -vcirc2*y2/xy_radius2/VelocityUnits*cos(theta2)-(vcirc2*x2/xy_radius2/VelocityUnits)*cos(theta2)*sin(phi2)
                                                    +GalaxyVelocity2[0]*1e5/VelocityUnits;
            thisgrid->BaryonField[Vel2Num][index] = -vcirc2*y2/xy_radius2/VelocityUnits*sin(theta2)+(vcirc2*x2/xy_radius2/VelocityUnits)*cos(theta2)*cos(phi2)
                                                    +GalaxyVelocity2[1]*1e5/VelocityUnits;
            thisgrid->BaryonField[Vel3Num][index] = (vcirc2*x2/xy_radius2/VelocityUnits)*sin(theta2)
                                                    +GalaxyVelocity2[2]*1e5/VelocityUnits;
        }//Galaxy2
          
        // Metallicity assignment
        if (TestProblemData.UseMetallicityField) {
            thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
                                                     TestProblemData.MetalFractionByMass * DiskMetallicity;
        }
      }
      
      // Energy assignment
      if(HydroMethod != Zeus_Hydro) {
        thisgrid->BaryonField[TENum][index] +=  0.5 *
        (POW(thisgrid->BaryonField[Vel1Num][index],2) +
        POW(thisgrid->BaryonField[Vel2Num][index],2) +
        POW(thisgrid->BaryonField[Vel3Num][index],2));
      }
      
      if (StarMakerTypeIaSNe) {
          int SNIaNum = FindField(MetalSNIaDensity , thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
          if(SNIaNum != -1) {
              thisgrid->BaryonField[SNIaNum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
          }
      }
      if (StarMakerTypeIISNeMetalField) {
          int SNIINum = FindField(MetalSNIIDensity , thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
          if(SNIINum != -1) {
              thisgrid->BaryonField[SNIINum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
          }
          else {
              ENZO_FAIL("Thought we would find a SNII field but did not.");
          }
      }

	  if(TestProblemData.MultiSpecies)
	  {
	    thisgrid->BaryonField[HINum][index] = TestProblemData.HI_Fraction *
	      TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HIINum][index] = TestProblemData.HII_Fraction *
	      TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeINum][index] = TestProblemData.HeI_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeIINum][index] = TestProblemData.HeII_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeIIINum][index] = TestProblemData.HeIII_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    if(TestProblemData.MultiSpecies > 1){
	      thisgrid->BaryonField[HMNum][index] = TestProblemData.HM_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	      thisgrid->BaryonField[H2INum][index] = 2 * TestProblemData.H2I_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	      thisgrid->BaryonField[H2IINum][index] = 2 * TestProblemData.H2II_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];
	    }

	    if (TestProblemData.MultiSpecies > 1)
	      thisgrid->BaryonField[HIINum][index] -=
		(thisgrid->BaryonField[HMNum][index] + thisgrid->BaryonField[H2IINum][index]
		 + thisgrid->BaryonField[H2INum][index]);

	    // Electron "density" (remember, this is a factor of m_p/m_e scaled
	    // from the 'normal' density for convenience) is calculated by
	    // summing up all of the ionized species.  The factors of 0.25 and
	    // 0.5 in front of HeII and HeIII are to fix the fact that we're
	    // calculating mass density, not number density (because the
	    // thisgrid->BaryonField values are 4x as heavy for helium for a single
	    // electron)
	    thisgrid->BaryonField[DeNum][index] = thisgrid->BaryonField[HIINum][index] +
	      0.25*thisgrid->BaryonField[HeIINum][index] +
	      0.5*thisgrid->BaryonField[HeIIINum][index];

	    if (TestProblemData.MultiSpecies > 1)
	      thisgrid->BaryonField[DeNum][index] += 0.5*thisgrid->BaryonField[H2IINum][index] -
		thisgrid->BaryonField[HMNum][index];

	    // Set deuterium species (assumed to be a negligible fraction of the
	    // total, so not counted in the conservation)
	    if(TestProblemData.MultiSpecies > 2){
	      thisgrid->BaryonField[DINum ][index] =
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HINum][index];
	      thisgrid->BaryonField[DIINum][index] =
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HIINum][index];
	      thisgrid->BaryonField[HDINum][index] = 0.75 *
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[H2INum][index];
	    }
	  } // if(TestProblemData.MultiSpecies)

	    if (HydroMethod == MHD_RK)
	      {
		thisgrid->BaryonField[B1Num][index] = Bfield[0];
		thisgrid->BaryonField[B2Num][index] = Bfield[1];
		thisgrid->BaryonField[B3Num][index] = Bfield[2];

		thisgrid->BaryonField[TENum][index] += 
		  0.5*(POW(thisgrid->BaryonField[B1Num][index], 2) +
		       POW(thisgrid->BaryonField[B2Num][index], 2) + 
		       POW(thisgrid->BaryonField[B3Num][index], 2))/thisgrid->BaryonField[DensNum][index];
	      }



	} // i
      } // j
    } // k

    return SUCCESS;

  } // InitializeGrid

  void InitializeParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid,
			  TopGridData &MetaData)
  {
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    mt_init(thisgrid->ID);

    if(debug)
      printf("Entering AgoraRestart InitializeParticles\n");

    // Determine the number of particles of each type
    int nBulge, nDisk, nHalo, nParticles;
    nBulge = 2*nlines("bulge.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Bulge Particles %"ISYM"\n", nBulge);
    nDisk = 2*nlines("disk.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Disk Particles %"ISYM"\n", nDisk);
    nHalo = 2*nlines("halo.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Halo Particles %"ISYM"\n", nHalo);
    nParticles = nBulge + nDisk + nHalo;
    if(debug) fprintf(stderr, "InitializeParticles: Total Number of Particles %"ISYM"\n", nParticles);


    // Initialize particle arrays
    PINT *Number = new PINT[nParticles];
    int *Type = new int[nParticles];
    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION];
    for (int i = 0; i < thisgrid->GridRank; i++)
    {
      Position[i] = new FLOAT[nParticles];
      Velocity[i] = new float[nParticles];
    }
    float *Mass = new float[nParticles];
    float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    for (int i = 0; i < NumberOfParticleAttributes; i++)
    {
      Attribute[i] = new float[nParticles];
      for (int j = 0; j < nParticles; j++)
	Attribute[i][j] = FLOAT_UNDEFINED;
    }

    FLOAT dx = thisgrid->CellWidth[0][0];

    // Read them in and assign them as we go
    int count = 0;
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "bulge.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "disk.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx);

    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    MetaData.NumberOfParticles = count;
    if(debug) fprintf(stderr, "InitializeParticles: Set Number of Particles %"ISYM"\n", count);

  }

  float gauss_mass(
    float RhoZero, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT cellwidth)
  {
    // Computes the total mass in a given cell by integrating the density
    // profile using 5-point Gaussian quadrature.
    // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
    FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
    FLOAT xResult [5];
    FLOAT yResult [5];
    FLOAT r, z;
    float Mass = 0;
    int i,j,k;

    for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      for (j=0;j<5;j++)
      {
	yResult[j] = 0.0;
	for (k=0;k<5;k++)
	{
    // Dec 2021 JL isothermal initialization
      r = sqrt((POW(xpos+EvaluationPoints[i]*cellwidth/2.0, 2.0) + POW(ypos+EvaluationPoints[j]*cellwidth/2.0, 2.0) + POW(zpos+EvaluationPoints[k]*cellwidth/2.0, 2.0) ) );
	  z = fabs(zpos+EvaluationPoints[k]*cellwidth/2.0);
	  yResult[j] +=
	    cellwidth/2.0 * Weights[k] * RhoZero / (1 + POW((r/this->ScaleLength), 2.0));
	}
	xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
      }
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }
    return Mass;
  }

  void ReadInVcircData(void)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int i=0, ret;
    float vcirc;
    FLOAT rad;

    fptr = fopen("vcirc.dat" , "r");

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret += sscanf(line, "%"PSYM" %"FSYM, &rad, &vcirc);
      this->VCircRadius[i] = rad*kpc_cm; // 3.08567758e21 = kpc/cm
      this->VCircVelocity[i] = vcirc*1e5; // 1e5 = (km/s)/(cm/s)
      i += 1;
    }

    fclose(fptr);
  } // ReadInVcircData
  
  void ReadInVcircData2(void)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int i=0, ret;
    float vcirc;
    FLOAT rad;

    fptr = fopen("vcirc1.dat" , "r");

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret += sscanf(line, "%"PSYM" %"FSYM, &rad, &vcirc);
      this->VCircRadius2[i] = rad*kpc_cm; // 3.08567758e21 = kpc/cm     
      this->VCircVelocity2[i] = vcirc*1e5; // 1e5 = (km/s)/(cm/s)
      i += 1;
    }

    fclose(fptr);
  } // ReadInVcircData2
  
  float InterpolateVcircTable(FLOAT radius)
  {
    int i;

    for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
      if (radius < this->VCircRadius[i])
	break;

    if (i == 0)
      return (VCircVelocity[i]) * (radius - VCircRadius[0]) / VCircRadius[0];
    else if (i == VCIRC_TABLE_LENGTH)
      ENZO_FAIL("Fell off the circular velocity interpolation table 1");

    // we know the radius is between i and i-1
    return VCircVelocity[i-1] +
      (VCircVelocity[i] - VCircVelocity[i-1]) *
      (radius - VCircRadius[i-1])  /
      (VCircRadius[i] - VCircRadius[i-1]);
  }

  float InterpolateVcircTable2(FLOAT radius2)
  {
    int i;

    for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
      if (radius2 < this->VCircRadius2[i])
        break;
    
    // printf("%f",radius2)
    // printf("%f",VcircRadius[i])
    
    if (i == 0)
      return (VCircVelocity2[i]) * (radius2 - VCircRadius2[0]) / VCircRadius2[0];
    else if (i == VCIRC_TABLE_LENGTH)
      ENZO_FAIL("Fell off the circular velocity interpolation table 2");

    // we know the radius is between i and i-1                                                                                  
    return VCircVelocity2[i-1] +
      (VCircVelocity2[i] - VCircVelocity2[i-1]) *
      (radius2 - VCircRadius2[i-1])  /
      (VCircRadius2[i] - VCircRadius2[i-1]);
  }

  int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int ret;
    FLOAT x, y, z, x1, y1, z1, x2, y2, z3, x_original;
    float vx, vy, vz, vx1, vy1, vz1, vx2, vy2, vz2;
    double mass;
    float theta,phi,Rx,Ry,Rz,Rvx,Rvy,Rvz, theta2,phi2,Rx2,Ry2,Rz2,Rvx2,Rvy2,Rvz2;
    
    theta=this->DiskTheta*pi/180.0;
	theta2=this->DiskTheta2*pi/180.0;
	phi=this->DiskPhi*pi/180.0; 
	phi2=this->DiskPhi2*pi/180.0;
      
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1, TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    fptr = fopen(fname, "r");

    while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret +=
	sscanf(line,
	       "%"PSYM" %"PSYM" %"PSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	       &x, &y, &z, &vx, &vy, &vz, &mass);
      
      x_original=x;
      
      if (x_original>=0){
          x-=(this->CenterPosition[0]-0.5)*LengthUnits/kpc_cm;
          y-=(this->CenterPosition[1]-0.5)*LengthUnits/kpc_cm;
          z-=(this->CenterPosition[2]-0.5)*LengthUnits/kpc_cm;
          vx-=GalaxyVelocity[0];
          vy-=GalaxyVelocity[1];
          vz-=GalaxyVelocity[2];
          
          Rx=x*cos(theta)+(-y*cos(theta)+z*sin(theta))*sin(phi);
          Ry=x*sin(theta)+(y*cos(theta)-z*sin(theta))*cos(phi);
          Rz=y*sin(theta)+z*cos(theta);

          Rvx=vx*cos(theta)+(-vy*cos(theta)+vz*sin(theta))*sin(phi);
          Rvy=vx*sin(theta)+(vy*cos(theta)-vz*sin(theta))*cos(phi);
          Rvz=vy*sin(theta)+vz*cos(theta);
      
          Position[0][c] = Rx * kpc_cm / LengthUnits + this->CenterPosition[0];
          Position[1][c] = Ry * kpc_cm / LengthUnits + this->CenterPosition[1];
          Position[2][c] = Rz * kpc_cm / LengthUnits + this->CenterPosition[2];

          Velocity[0][c] = Rvx * km_cm / VelocityUnits + GalaxyVelocity[0]*km_cm/VelocityUnits;
          Velocity[1][c] = Rvy * km_cm / VelocityUnits + GalaxyVelocity[1]*km_cm/VelocityUnits;
          Velocity[2][c] = Rvz * km_cm / VelocityUnits + GalaxyVelocity[2]*km_cm/VelocityUnits;
      } // galaxy 1
      
      else if (x_original<0){
          x-=(this->CenterPosition2[0]-0.5)*LengthUnits/kpc_cm;
          y-=(this->CenterPosition2[1]-0.5)*LengthUnits/kpc_cm;
          z-=(this->CenterPosition2[2]-0.5)*LengthUnits/kpc_cm;
          vx-=GalaxyVelocity2[0];
          vy-=GalaxyVelocity2[1];
          vz-=GalaxyVelocity2[2];
          
          Rx2=x*cos(theta2)+(-y*cos(theta2)+z*sin(theta2))*sin(phi2);
          Ry2=x*sin(theta2)+(y*cos(theta2)-z*sin(theta2))*cos(phi2);
          Rz2=y*sin(theta2)+z*cos(theta2);

          Rvx2=vx*cos(theta2)+(-vy*cos(theta2)+vz*sin(theta2))*sin(phi2);
          Rvy2=vx*sin(theta2)+(vy*cos(theta2)-vz*sin(theta2))*cos(phi2);
          Rvz2=vy*sin(theta2)+vz*cos(theta2);

          Position[0][c] = Rx2 * kpc_cm / LengthUnits + this->CenterPosition2[0];
          Position[1][c] = Ry2 * kpc_cm / LengthUnits + this->CenterPosition2[1];
          Position[2][c] = Rz2 * kpc_cm / LengthUnits + this->CenterPosition2[2];

          Velocity[0][c] = Rvx2 * km_cm / VelocityUnits + GalaxyVelocity2[0]*km_cm/VelocityUnits;
          Velocity[1][c] = Rvy2 * km_cm / VelocityUnits + GalaxyVelocity2[1]*km_cm/VelocityUnits;
          Velocity[2][c] = Rvz2 * km_cm / VelocityUnits + GalaxyVelocity2[2]*km_cm/VelocityUnits;
      } // galaxy 2

      // Particle masses are actually densities.                                                                                                                                                                                
      Mass[c] = mass * 1e9 * SolarMass / MassUnits / dx / dx / dx;
      Type[c] = particle_type;
      Number[c] = c++;
      }
    fclose(fptr);
    // printf("after2:c=%d,\n",c);
    return c;
  } // ReadParticlesFromFile

}; // class declaration


//.. register:
namespace {
    EnzoProblemType_creator_concrete<ProblemType_AgoraRestart>
        agora_restart("AgoraRestart");
}




#endif

## Modeling growth and metabolic exchange of multispecies community in mixed batch culture

This protocol demonstrates the capacity of COMETS to simulate multispecies community dynamics and metabolic exchange using the MATLAB toolbox.

The metabolic models used in this simulation modelsCommunity.mat and the growth medium mediumCommunity.mat can be downloaded from the [COMETS GitHub repository](https://github.com/segrelab/COMETS_Protocols/tree/master/COMETS_protocols/COMETS_example_MultispeciesGrowthExchange). In this example, the growth phenotypes of *Bacillus subtilis* [1], *Escherichia coli* [2], *Klebsiella pneumoniae* [3], *Lactococcus lactis* [4], *Methylobacterium extorquens* [5], *Pseudomonas aeruginosa* [6], *Porphyromonas gingivalis* [7], *Rhodobacter sphaeroides* [8], *Shigella boydii* [9], *Saccharomyces cerevisiae* [10], *Salmonella enterica* [11], *Shewanella oneidensis* [12], *Synechocystis* sp. PCC6803 [13], and *Zymomonas mobilis* [14] are simulated. These organisms, whose growth and metabolic exchange phenotypes were previously analyzed in silico [15], were selected to represent a broad cross-section of nutrient utilization capabilities. Here, the growth profile of a multispecies combination of these organisms in the presence of D-glucose and L-alanine over 12 hours will be analyzed.

### Import the models and medium conditions, and initialize the COMETS layout

Unzip the ‘modelsCommunity.mat.zip’ file. If this is unzipped as a directory, it will be necessary to move the 'modelsCommunity.mat' file back into the top directory.
In MATLAB, import the metabolic model MATLAB structure.

```matlab
>> load modelsCommunity.mat
```

This will result in a structure models being loaded into the workspace.

Initialize a COMETS layout object using the CometsLayout class and add the models using the addModel function.

```matlab
>> layout = CometsLayout();
modelNames = fieldnames(models);
for m = 1:length(modelNames)
layout = addModel(layout,models.(modelNames{m}));
end
```

This step results in a CometsLayout object layout, with the attributes below, being added to the workspace.
layout = 

```matlab
  CometsLayout with properties:

                  models: {1×14 cell}
                    xdim: 1
                    ydim: 1
                    mets: {725 cell}
               media_amt: [725×1 double]
                  params: [1×1 CometsParams]
     diffusion_constants: [725×2 double]
    global_media_refresh: [725×1 double]
           media_refresh: [725×1 double]
     global_static_media: [725×2 double]
            static_media: [725×1 double]
           initial_media: 0
                 barrier: 0
             initial_pop: 0
           external_rxns: [0×0 table]
       external_rxn_mets: [0×0 table]
```

Here, the object property models is an array with the individual metabolic models.

Load the medium file.

```matlab
>> load mediumCommunity.mat
```

The mediumCommunity file contains four objects that will be added to the workspace:

minMed: A cell array of 32 molecules (in BIGG format) that are contained in the minimal medium to be added at nonlimiting concentrations:

```matlab
minMed =

  32×1 cell array

    '4abz[e]'
    'btn[e]'
    'ca2[e]'
    'cbl1[e]'
    'chol[e]'
    'cl[e]'
    'cobalt2[e]'
    'cu2[e]'
    'fe2[e]'
    'fe3[e]'
    'fol[e]'
    'h2[e]'
    'h2o[e]'
    'k[e]'
    'lipoate[e]'
    'mg2[e]'
    'mn2[e]'
    'mobd[e]'
    'na1[e]'
    'ncam[e]'
    'nh4[e]'
    'ni2[e]'
    'no3[e]'
    'o2[e]'
    'pi[e]'
    'pnto-R[e]'
    'pydx[e]'
    'ribflv[e]'
    'slnt[e]'
    'so4[e]'
    'thm[e]'
    'zn2[e]'
```

Nutrients: A cell array of five carbon sources (in BIGG format) that will be added at limiting concentrations:

```matlab
nutrients =

  1×2 cell array

    'glc-D[e]'    'ala-L[e]'
```

nutrientNames: A cell array containing human-readable names of the five carbon sources:

```matlab
nutrientNames =

  1×5 cell array

    'D-glucose'    ‘L-alanine’
```

numCarbons: A numerical vector containing the number of carbon atoms contained in each of the five carbon sources:

```matlab
numCarbons =

     6     3
```

Add the minimal media components at nonlimiting concentrations using the setMedia function.

```matlab
>> for mm = 1:length(minMed)
       layout = layout.setMedia(minMed{mm},1000);
   end
```

The quantity 1000 is used to denote 1000mmol, an amount of molecule that is effectively non limiting for this simulation. If a continuous culture environment is desired, the setStaticMedia function can be used within the for loop to fix the amount of each minimal medium molecule at 1000mmol.

```matlab
>> layout = setStaticMedia(layout,1,1,minMed{mm},1000);
```

Add the carbon sources at limiting amounts (5e-4 mmol) using the setMedia function.

```matlab
>> for n = 1:length(nutrients)
       layout = layout.setMedia(nutrients{n},5e-4);
   end
```

Alternatively, each nutrient can be added at equal carbon ratios within the for loop:

```matlab
>> layout = layout.setMedia(nutrients{n},5e-4/numCarbons(n)/length(nutrients));
```

If a continuous culture environment is desired, the global_media_refresh layout parameter can be used within the for loop in addition to the setMedia function to refresh the amount of each carbon source.

```matlab
>> layout.global_media_refresh(find(ismember(layout.mets,nutrients(n)))) = 50/(12/0.1)
```

The global_media_refresh  parameter defines the amount of additional nutrient to be added at each time step, analogous to a constant influx of fresh nutrients in a continuous culture device.

Define the initial organism populations.

```matlab
>> layout.initial_pop = ones(length(modelNames),1).*1.e-7;
```

Here, each organism is added to the environment at an abundance of 1.0·10-7  grams dry weight. By default, the layout is initialized with dimensions of 1 x 1 cells to simulate a well-mixed environment.

### Set the COMETS parameters and run the simulation

Define the COMETS working directory, log file names, and simulation parameters.

```matlab
>> cometsDirectory = 'CometsRunDir';

layout.params.writeBiomassLog = true;
layout.params.biomassLogRate = 1;
layout.params.biomassLogName = 'biomassLog';
layout.params.biomassLogFormat = 'MATLAB';
layout.params.writeMediaLog = true;
layout.params.mediaLogRate = 1;
layout.params.mediaLogName = 'mediaLog';
layout.params.mediaLogFormat = 'MATLAB';
layout.params.writeFluxLog = true;
layout.params.fluxLogRate = 1;
layout.params.fluxLogName = 'fluxLog.m';
layout.params.fluxLogFormat = 'MATLAB';

layout.params.maxSpaceBiomass = 1e3;
layout.params.timeStep = 0.01;  
layout.params.maxCycles = 1200;
layout.params.deathRate = 0.1;
```

Here, the simulation records the biomass of each organism, the amounts of medium components, and the fluxes of all organisms at each time step. The time step is set to 0.01 1/h and the maximum number of cycles is set to 2400 for a total simulation time of 24 hours. The death rate is set to a fixed value of 0.1, under which 10% of the population of each organism will be eliminated from the simulation at each time step. This death rate can also be set to approximate an organism dilution rate in a continuous culture environment.

Prepare the metabolic models for the COMETS simulation. Here, the lower bounds of certain reactions within the models are being altered in order to allow uptake of medium components. Specifically, exchange reactions associated with uptake of the minimal medium components are being unconstrained (lower bound set to -1000), while the exchange reactions associated with uptake of the limiting carbon sources are being opened (lower bound set to -10). The reactions to be unconstrained are located on a model-by-model basis by first identifying the index of a specific external metabolite and matching the appropriate exchange reaction to it via the S matrix.

```matlab
for m = 1:length(modelNames)
    modelCurr = models.(modelNames{m});
    minMedMets = find(ismember(modelCurr.mets,minMed));
    for i = 1:length(minMedMets) 
modelCurr.lb(intersect(find(findExcRxns(modelCurr)),find(modelCurr.S(minMedMets(i),:)))) = -1000; % Allow unlimited uptake of nonlimiting nutrients
    end
    limitingMets = find(ismember(modelCurr.mets,nutrients));
    for i = 1:length(limitingMets) 
modelCurr.lb(intersect(find(findExcRxns(modelCurr)),find(modelCurr.S(limitingMets(i),:)))) = -10; % Allow limited uptake of limiting nutrients
    end
    models.(modelNames{m}) = modelCurr;
end
```

Run the COMETS simulation.

```matlab
>> runComets(layout,cometsDirectory)
```

### Parse the COMETS output logs and visualize data 

Parse the media log using the parseBiomassLog function and format the log into a matrix.

```matlab
>> biomassLogRaw = parseBiomassLog([cometsDirectory '/' layout.params.biomassLogName]);
biomassLog = zeros(size(biomassLogRaw,1)/length(modelNames),length(modelNames));
for i = 1:length(modelNames)
    biomassLog(:,i) = biomassLogRaw.biomass(i:length(modelNames):end);
end
```

Format the names of the organisms for plotting and plot the biomass over time.

```matlab
>> modelNamesFormatted = cell(length(modelNames),1);
for m = 1:length(modelNames)
    s = split(modelNames{m},'_');
    modelNamesFormatted{m} = [s{1} '. ' s{2}];
end
 
close all
figure
plotColors = parula(length(modelNames));

for m = 1:length(modelNames)
plot([1:layout.params.maxCycles+1]*layout.params.timeStep,biomassLog(:,m),'LineWidth',4,’Color’,plotColors(m,:))
end

set(gca,'FontSize',16)
ylabel('Nutrient Amount (mmol)')
xlabel('Time (h)')
legend(modelNamesFormatted)
```

The results of this action are shown in Figure 1a below.

Parse the media log using the parseMediaLog function and format the log into a matrix.

```matlab
>> allMetsFromModels = layout.mets;
COMETSCycles = layout.params.maxCycles;
mediaLogMat = zeros(length(allMetsFromModels),COMETSCycles);

mediaLogRaw = parseMediaLog([cometsDirectory '/' layout.params.mediaLogName]);
mediaLogMetOrder = zeros(length(allMetsFromModels),1);

% Re-order the medium components to match the list in layout.mets
for i = 1:length(allMetsFromModels)
        mediaLogMetOrder(i) = find(ismember(mediaLogRaw.metname(1:length(allMetsFromModels)),allMetsFromModels(i)));
end

for i = 1:COMETSCycles
currentMedia = mediaLogRaw.amt(find(mediaLogRaw.t == i));
mediaLogMat(:,i) = currentMedia(mediaLogMetOrder);
end
```

The resulting mediaLogMat matrix has dimensions M x N, where M is the number of metabolites in the COMETS simulation and N is the number of time steps.

Plot the nutrient abundances.

```matlab
>> % Make a list of indices from allMetsFromModels ordered according to the elements in nutrients
nutrientsToPlot = zeros(length(nutrients),1);
for i = 1:length(nutrients)
    nutrientsToPlot(i) = find(ismember(allMetsFromModels,nutrients{i}));
end

figure
plot([1:layout.params.maxCycles]*layout.params.timeStep,mediaLogMat(nutrientsToPlot,:)', 'LineWidth',4)
set(gca,'FontSize',16)
ylabel('Nutrient Amount (mmol)')
xlabel('Time (h)')
legend(nutrientNames)
```

This action allows visualization of the nutrient abundances, shown in Figure 1b below.

These results, along with directly examining the nutrient abundances in mediaLogMat, allow us to infer that *M. extorquens* and *S. oneidensis* initially grew to high abundances on D-glucose and L-alanine, rapidly outcompeting the remaining organisms. We are also able to observe low growth of multiple other organisms, such as *S. boydii* and *Z. mobilis*, which peak upon exhausting D-glucose at hour 6. When both primary resources are exhausted between hours 6 and 7, the abundance of all organisms begins to decay.

Obtain the metabolites secreted, absorbed, and exchanged by each organism using the getSecAbsExcMets function.

```matlab
>> [secMets,absMets,excTable] = getSecAbsExcMets([cometsDirectory '/' layout.params.fluxLogName],models,layout);
```

This function outputs two matrices, secMets and absMets, each having dimensions N x M where N is the number of metabolites present in the layout and M is the number of metabolic models. For any model-metabolite pair, the quantity present in secMets or absMets denotes the highest flux that the metabolite was secreted or absorbed by the corresponding organism in the entire simulation. 

The function also outputs a matrix excTable, which contains information about each metabolite that was exchanged in the simulation. The first column of excTable contains the index (corresponding to the order in models) of the organism that secreted a metabolite, the second column contains the index of the organism that absorbed that metabolite, and the third column contains the index (corresponding to allMetsFromModels) of the metabolite that was exchanged. Here, a truncated version of excTable is as follows:

```matlab
excTable =

     2     5    17
     2    12    17
     3     5    17
     3    12    17
     6     5    17
     6    12    17
     9     5    17
     9    12    17
    11     5    17
    11    12    17
    14     5    17
    14    12    17
     6    12    36
     1     5    65
     2     5    65
     3     5    65
     5     5    65
```

The first line, for example, denotes that acetate (metabolite 17) is secreted by *E. coli* (organism 2) and absorbed by *M. extorquens* (organism 5).

Use the secMets matrix to identify molecule secretion.

Here, the source organism of each metabolite secreted during the simulation will be identified. First, a list of metabolites with nonzero secretion flux is generated.

```matlab
>>  nonzeroSecMetIndices = find(sum(secMets,2));
```

This action yields a list of indices that correspond to allMetsFromModels, which can be read using the command:

```matlab
>>  allMetsFromModels(nonzeroSecMetIndices);
```

This yields a list of metabolites, a truncated version of which is below:

```matlab
ans =

  18×1 cell array

    '5mtr[e]'
    'ac[e]'
    'ala-L[e]'
    'co2[e]'
    'fe2[e]'
    'for[e]'
    'glyclt[e]'
```

Molecules of particular interest can be selected and analyzed further. Here, secretion of acetate (ac[e]) and formate (for[e]) will be analyzed.

```matlab
>> selectSecMets = {'ac[e]','for[e]'};
selectSecMetIndices = intersect(find(ismember(allMetsFromModels,selectSecMets)),nonzeroSecMetIndices);
```

The abundances of each of these metabolites over time can be plotted:

```matlab
>> figure
plotColors2 = winter(2);
for s = 1:length(selectSecMets)
plot([1:layout.params.maxCycles]*layout.params.timeStep,smoothdata(mediaLogMat(selectSecMetIndices(s),:)'),'LineWidth',4,'Color',plotColors2(s,:))
    hold on
end
set(gca,'FontSize',16)
ylabel('Metabolite Amount (mmol)')
xlabel('Time (h)')
legend({'Acetate','Formate'})
```

This action yields Figure 1c below, showing rapid accumulation and subsequent consumption of acetate and formate.

![Figure 1](https://github.com/segrelab/comets-manual/blob/master/docs/img/multispecies_1.png)
**Figure 1. Growth and metabolic exchange of 14-species microbial community. a)** Biomass production of all 14 organisms over time. **b)** Consumption of limiting carbon sources over time. **c)** Secretion and consumption of metabolic byproducts over time.

### References
[1] Henry, C. S., Zinner, J. F., Cohoon, M. P. & Stevens, R. L. iBsu1103: a new genome-scale metabolic model of Bacillus subtilis based on SEED annotations. Genome Biol. 10, R69 (2009).

[2] Orth, J. D. et al. A comprehensive genome-scale reconstruction of Escherichia coli metabolism--2011. Mol. Syst. Biol. 7, 535 (2011).

[3] Liao, Y.-C. et al. An experimentally validated genome-scale metabolic reconstruction of Klebsiella pneumoniae MGH 78578, iYL1228. J. Bacteriol. 193, 1710–1717 (2011).

[4] Flahaut, N. A. L. et al. Genome-scale metabolic model for Lactococcus lactis MG1363 and its application to the analysis of flavor formation. Appl. Microbiol. Biotechnol. 97, 8729–8739 (2013).

[5] Peyraud, R. et al. Genome-scale reconstruction and system level investigation of the metabolic network of Methylobacterium extorquens AM1. BMC Syst. Biol. 5, 189 (2011).

[6] Oberhardt, M. A., Puchałka, J., Fryer, K. E., Martins dos Santos, V. A. P. & Papin, J. A. Genome-scale metabolic network analysis of the opportunistic pathogen Pseudomonas aeruginosa PAO1. J. Bacteriol. 190, 2790–2803 (2008).

[7] Mazumdar, V., Snitkin, E. S., Amar, S. & Segrè, D. Metabolic network model of a human oral pathogen. J. Bacteriol. 191, 74–90 (2009).

[8] Imam, S. et al. iRsp1095: a genome-scale reconstruction of the Rhodobacter sphaeroides metabolic network. BMC Syst. Biol. 5, 116 (2011).

[9] Monk, J. M. et al. Genome-scale metabolic reconstructions of multiple Escherichia coli strains highlight strain-specific adaptations to nutritional environments. Proc. Natl. Acad. Sci. U. S. A. 110, 20338–20343 (2013).

[10] Zomorrodi, A. R. & Maranas, C. D. Improving the iMM904 S. cerevisiae metabolic model using essentiality and synthetic lethality data. BMC Syst. Biol. 4, 178 (2010).

[11] Thiele, I. et al. A community effort towards a knowledge-base and mathematical model of the human pathogen Salmonella Typhimurium LT2. BMC Syst. Biol. 5, 8 (2011).

[12] Pinchuk, G. E. et al. Constraint-based model of Shewanella oneidensis MR-1 metabolism: a tool for data analysis and hypothesis generation. PLoS Comput. Biol. 6, e1000822 (2010).

[13] Nogales, J., Gudmundsson, S., Knight, E. M., Palsson, B. O. & Thiele, I. Detailing the optimality of photosynthesis in cyanobacteria through systems biology analysis. Proc. Natl. Acad. Sci. U. S. A. 109, 2678–2683 (2012).

[14] Lee, K. Y., Park, J. M., Kim, T. Y., Yun, H. & Lee, S. Y. The genome-scale metabolic network analysis of Zymomonas mobilis ZM4 explains physiological features and suggests ethanol and succinic acid production strategies. Microb. Cell Fact. 9, 94 (2010).

[15] Pacheco, A. R., Moel, M. & Segrè, D. Costless metabolic secretions as drivers of interspecies interactions in microbial ecosystems. Nat. Commun. 10, 103 (2019).

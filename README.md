# Free energy calculations (estimations)

Shiny app that allows quick free energy calculations of metabolisms of interest based on physicochemical parameters. 
- Outputs are tabular and plotted
- Activity coefficient estimates for environments of various Ionic Strengths come from Amend and LaRowe (2019), linked in this repo.

Limitations:
1. This is just an estimation! Although I have found it to be decently accurate. Values should be compared with specific results from software like GWB.
2. Can only calculate reactions with species for which there is thermodynamicdata within CHNOSZ.
3. At the moment, results are reported in KJ/mol and should be manually normalized depending on question on hand.
  - Next version will have automatic normalization per electron transferred.

## Usage

1. After cloning the repo, please run `renv::restore()` in a R session to ensure the required packages and versions are installed.
   
2. Run app: `Rscript app.R`

3. Upload a physicochemical .csv file that follows the template (there is a link for the template within the app). A couple of notes about the file:
  - First 4 column must remain as: Sample, Temperature (C), Pressure (bar), pH.
  - Column 5+ can be changed but a column with *H2O must always be present with value of 1* for every row.
  - Species names must follow [CHNOSZ](https://chnosz.net/vignettes/anintro.html) nomeculature.

4. Select the ionic strength of you environment.

5. Fill out the metabolic reactions you want to calculate (these must be PROPERLY BALANCED) in this format:
> nameOfChoice = speciesInReaction(comma,separated)|speciesCoefficients(comma,separated) <br><br>
> speciesCoefficients rules: negative means the spp is a reactant, positive means the spp is a product <br><br>
> example: FeOx_NitrateRed = Fe+2,NO3-,H+,Fe+3,NO2-,H2O|-2,-1,-2,2,1,1

5. Fill out any minerals that are used in the reactions
  - Leave blank if none

6. Click calculate.
  - Table and plot should appear.
  - Can customize y-axis of plot with the dropdown.
  

Please contact Didi Bojanova at diana.p.bojanova@gmail.com with any questions.

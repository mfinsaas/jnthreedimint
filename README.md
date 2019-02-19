# Visualizing Continuous x Continuous Interactions in 3D
#### M. C. Finsaas & B. L. Goldstein, Department of Psychology, Stony Brook University

#### https://mfinsaas.shinyapps.io/3DInteractions/

Please email Megan Finsaas at megan.finsaas@stonybrook.edu if you have any questions or comments about the program. 

This program allows researchers to visualize interactions between two continuous variables in 3D space along with the 
results of the Johnson-Neyman test.

## Notes about use:
### Browser recommendation
- We highly recommend using Google Chrome as your browser when using this program. This will increase the 
likelihood that you can see the plot fully. 

### Trouble seeing the plot?
- If you do run into trouble viewing the plot (e.g., the text boxes overlap), please...
    1. Double check that you're using Google Chrome.
    2. Expand your browser window fully.
    3. Zoom out in your browser window. This option is usually located under "View." (on Macs: the shortcut is "command" + "-")
    4. Use the Plotly options to zoom, pan, and rotate the plot. These options appear at the top right side of the plot when you hover
    your mouse over the plot. Once you have clicked on zoom, pan, or rotate, click on the plot and hold while moving the mouse side to side
    or up and down to change the view. 
    
If you continue to have problems, please email me at megan.finsaas@stonybrook.edu.

### Saving the plot
- To save a 2D image of your plot, you can either take a screen shot (on Macs, command + shift + 4; on PCs, 
Alt + PrtScn) or you can use the Plotly option "Download Plot as a png" (this option is located at the top of the plot frame as a
camera icon. The former will allow you to also capture the model information, Johnson-Neyman results, and crossover points. 
If you wish to save a rotatable version of the plot, we recommend using taking a video screenshot while rotating the plot manually 
(on Macs: https://support.apple.com/en-us/HT208721; on PCs: https://www.pcmag.com/news/349410/how-to-capture-video-clips-in-windows-10). 

### Accepted file types
- Currently, the program accepts SPSS (.sav), Excel (.xlsx), and .csv files. In some cases, some of these file types 
may appear grayed out/unavailable in the file upload window. To get around this, simply change your setting in the file upload window 
to view all file types.

### Reactivity
1. The program will estimate the model and provide any plots as soon as you enter three unique variables for the 
predictors and outcome. After this point, the program will automatically update the model and plot any time you update a predictor,
covariate, or outcome. You can verify which variables are in your model by looking at the Model Results on the 3D Plot tab.

2. There are several plot features that can be toggled on and off. The scatterplot, regression plane, and 95% CI around the 
predicted values are available for all models, regardless of whether an interaction term is significant. They will appear in the sidebar 
as soon as two unique predictors and an outcome variable are entered. Other plot features are only available if the interaction term is 
significant and if certain values fall within the range of the observed data; these include the crossover points (or lines, as they 
appear in 3D space) and the shaded Johnson-Neyman regions of significance (ROS). Additionally, the option to include a gradient that reflects the 
width of the 95% confidence band around the slope estimate appears if the respective ROS box is checked.

### Missing values
- At the current time, the program only accepts blanks as missing values. Please recode all missing values as blanks 
before inputting your data. If you opt to center or standardize your variables in the program, it will do so using only those cases that 
have values present on all variables in the regression model. 

### Higher order terms
- At the current time, the program does not allow for squared or other higher order terms.

### Abbreviations
    ROS: region of significance resulting from the Johnson-Neyman test
    CI: confidence interval (in this case, around the predicted values)
    CB: confidence band (in this case, around the slope estimate)

The tool was created using R (R Development Core Team, 2016) in the Shiny framework (Chang, Cheng, Allaire, Xie, & McPherson, 2017), 
and the three-dimensional plotting was done using the plotly package (Sievert et al., 2017). To use the tool, researchers must 
upload a data file and specify the variables that act as the predictors and the outcome in their model, as well as any covariates. 
It is not necessary for the user to first run their regression model in a different program and pull out specific values. The primary 
output is the 3D regression plane with shaded regions of significance. Our tool also has several optional features that can be toggled 
on and off: observed raw data as a scatterplot, which can be especially helpful for identifying multivariate outliers; 95% confidence 
intervals around the predicted values as additional planes floating above and below the primary regression plane; slopes corresponding 
to the crossover points for both variables; and 95% confidence bands for the slope estimates in the form of a gradient over the regions 
of significance. In addition to the 3D regression plane, the tool outputs descriptive statistics on raw and standardized/centered data, 
the results of the regression model, the Johnson-Neyman values demarcating the regions of significance, the direction of the slope and 
percentage of cases within these regions and whether these regions fall within the observed data, Johnson-Neyman marginal effects plots, 
and tables containing the slope values and 95% confidence bands at all values of the moderator. 

As a brief overview, the tool produces the plot by dividing the predictors into 100 equally spaced cells and 
adding the Johnson Neyman values to these vectors (if they fall within the range of observed data). The Johnson-Neyman 
values and confidence bands are derived using the jtools package (Long, 2018). These vectors are then entered into the
regression model to produce a matrix of predicted values, with any covariates set to their means. This matrix of predicted 
values makes up the regression plane. The values that demarcate the Johnson-Neyman regions of significance are used to 
identify the areas on the plane to shade. The axes of the plot are determined by the range of the observed data, which, 
as noted above, can also be overlaid on the plane as a 3D scatterplot. See the "app" file in this repository for the source code.






Primary package citations
Long, J. A. (2018). jtools: Analysis and presentation of social scientific data. 
R package version 1.1.0, https://cran.r.project.org/package=jtools.

Sievert, C., Parmer, C., Hocking, T., Chamberlain, S., Ram, K., Corvellec, M., & Despouy, P. (2017). 
plotly: Create Interactive Web Graphics via ‘plotly.js’. R package version 4.7.1, https://CRAN.R-project.org/package=plotly.

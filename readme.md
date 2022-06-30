# GeothermoTool: An open-source software for basic geothermal calculations and plots

# Introduction

GeothermalTool provides thermal conductivity correction, temperature log analysis, the estimation of the potential of geothermal reservoirs, temperature field computing (1D and 2D), and Tectono-Thermal Evolution modeling. Besides, the calculation results can be visualizated by the Matplotlib module and exported to files including MS EXCEL and Comma-Separated Values (CSV). The plots can be saved as bitmaps and vectorgraphs allowing second edit. Moreover, GeothermoTool offers a standalone PC application with a simple interface and is easily to be run. Numpy, Scipy and Numba are used to speed up the calculations.  Users who are capable of programming in Python can join the development team, fix errors, and add more functions.


# Getting started

## 1) GeothermoTool offers a standalone PC application with a simple interface and is easily to be run. 
All packages are bundled in the applicatuon. It's a recommand way to run the software.
## 2) Run the source code. Required modules: Numpy, Scipy, Pandas, Matplotlib, Pwlf, openpyxl，and PysimpleGUI. Anaconda is a good choice if you want to run the source code. 

# Manual and Publication

GeothermoTool has a truly simple and intuitive interface. The user can easily get the results after reading the exampes of the paper.

# Input and output

## Files
GeothermoTool uses MS Excel XLSX and CSV (Comma Separated Values) files formats for the input and output files. The input data file structure is given by the template files. The headings must remain unchanged because GeothermalTool parses the data by the headings. After calculations, GeothermoTool can write the results to XLSX or CSV files. 

## Figures
 The plotting results can be saved in various file formats, including bitmap and vectorgraph, which can be reedited in Adobe Illustrator and CorelDraw. For batch computing, we just need to choose the data folder. GeothermoTool will automatically create a folder named “[name of input folder]+_output” and put the calculated results and plotted figures into the file.



# License

GeothermoTool is distributed under the GNU General Public License, version 3:
http://www.gnu.org/copyleft/gpl.html

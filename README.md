Analysis of Dendrite/Axon Exit Position
===================
by Album Shen

Categorizes the exit position of primary dendrites and axon from the soma into 10 vertical locations.

### Data needed: ###
*	Scaling information used to construct neurites model
*	Scaling information of the image (pixel to um conversion)
*	Lowest and highest slice in which cell body appears

### Parameters: ###
* 	['PA.py', input, x_d, y_d, z_d, x_s, y_s, z_s, s_l, s_h]
* 	input:		.swc file to be read; results will be output to the filename plus the extension '_output.txt'
* 	scale_default:	refers to the scaling used when plotting the neurites in NeuronStudio
* 	scale_actual:	refers to the actual scaling of the images
* 	s_l:		lowest slice in which the cell body appears
* 	s_h:		highest slice in which the cell body appears

### Assumptions about neurites file: ###
*	Soma is indexed with integer label n = '1', and neuronal type T = '1')
*	Axon is labeled with neuronal type T = '2'
*	All dendritic ends are labeled with neuronal type T = '6'
*	There is a single point of origin (one soma with T = '1')

### DO ONE OF THE FOLLOWING BEFORE USING: ###
* Uncomment line below "# Sample run" to allow for manual entry of input filename and pertinent scales directly into the script. Script can be run using "python P10.py" at the command line or simply clicking to execute.

* Uncomment line below "# Running through command line" to run by passing parameters through command line. 
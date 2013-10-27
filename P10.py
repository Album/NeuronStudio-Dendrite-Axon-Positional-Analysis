# Album Shen
# Deans Lab
# 2013 June
#
# Analysis of Dendrite/Axon Exit Position
#
# Data needed:
#	Scaling information used to construct neurites model
#	Scaling information of the image (pixel to um conversion)
#	Lowest and highest slice in which cell body appears
#
# Parameters:
# 	['PA.py', input, x_d, y_d, z_d, x_s, y_s, z_s, s_l, s_h]
# 	input:			.swc file to be read;
#					results will be output to the filename plus the extension '_output.txt'
# 	scale_default:	refers to the scaling used when plotting the neurites in NeuronStudio
# 	scale_actual:	refers to the actual scaling of the images
# 	s_l:			lowest slice in which the cell body appears
# 	s_h:			highest slice in which the cell body appears
#
# Assumptions about neurites file:
#	Soma is indexed with integer label n = '1', and neuronal type T = '1')
#	Axon is labeled with neuronal type T = '2'
#	All dendritic ends are labeled with neuronal type T = '6'
#	There is a single point of origin (one soma with T = '1')

import sys, string, os
from math import sqrt, pow, ceil

class PositionAnalysis:
	# initialize
	def __init__(self, params):
		global input, scale_default, scale_actual, slice_l, slice_h
		global output, output_name, neurites, ends, branches, dendrites, axons, axon_end, endpoints, end_um, end_dist
		global center, coor_swc, coor_image, coor_relative, coor_microns, coor_surface, num_pos, positions, width
		input = open(params[1], 'rb')
		output = open(params[1].replace('.swc', '_output_10.txt'), 'wb')
		output_name = (params[1].replace('.swc', '_output_10.txt'))
		scale_default = map(float, params[2:5])
		scale_actual = map(float, params[5:8])
		slice_l = float(params[8])
		slice_h = float(params[9])
		neurites = []
		ends = []
		branches = []
		dendrites = []
		axons = []
		center = []
		coor_swc = []
		coor_image = []
		coor_relative = []
		coor_microns = []
		coor_surface = []
		positions = []
		endpoints = []
		end_dist = []
		num_pos = 10
		width = 12
		
	# counts dendrites and axons for each branch
	# performs position analysis
	# outputs results to output file
	def work(self):

		import time
		start = time.clock()
		
		# neurites
		# n T x y z R P
		for line in input.readlines():
			if line.startswith('#'):
				continue
			neurites.append(line.split(" "))
		
		# correct for index at 0
		neurites[0] = [0] * 7
		
		# find ends
		# label indices of ends
		for n in neurites:
			if n[1] == '6':
				ends.append(n[0])
			if n[1] == '2':
				axon_end = n[0]
		
		# find primary dendrite branches
		for n in neurites:
			if n[6] == '1\r\n':
				branches.append(n[0])
		dendrites = [0] * len(branches)
		axons = [0] * len(branches)
		# print test
		# print '\nNumber of Ends:\t\t' + str(len(ends))
		# print 'Number of Branches:\t' + str(len(branches))

		for e in ends:
			prev = neurites[int(e)][6]
			while prev != '1\r\n':
				last = prev
				prev = self.traceback(int(prev))

			i = branches.index(str(int(last)))
			dendrites[i] += 1
		while True:
			try:
				prev = neurites[int(axon_end)][6]
				break
			except UnboundLocalError:
				print "The axon end could not be found. Please modify .swc file."
				print "Locate the neurite point of the axon end and set the neurite type (second element) to 2."
				sys.exit(1)
				
		while prev != '1\r\n':
			last = prev
			prev = self.traceback(int(prev))
		i = branches.index(str(int(last)))
		axons[i] += 1
		
		# extracting cell body coordinate from neurites file
		center = map(float, neurites[1][2:5])
		# print '\nCell body (um):\t\t' + str(center)
		midslice = float(slice_h + slice_l) / 2
		z_radius = ((slice_h - slice_l) / 2) * scale_actual[2]
		# cell body image coordinates
		center[0] = center[0] / scale_default[0]
		center[1] = center[1] / scale_default[1]
		center[2] = center[2] / scale_default[2] + 1
		# adjust z-center to midslice
		center[2] = midslice
		
		# print test
		'''
		print 'Cell body (px, px, sl):\t' + str(center)
		print 'Midslice:\t\t' + str(midslice)
		print 'Position height (slices):\t' + str(float(slice_h - slice_l) / num_pos)
		print 'Position height (um):\t' + str((float(slice_h - slice_l) / num_pos) * scale_actual[2])
		print 'Z-radius (um):\t\t' + str(z_radius)
		'''
		
		coor_swc = [0] * len(branches)
		coor_image = [0] * len(branches)
		coor_relative = [0] * len(branches)
		coor_microns = [0] * len(branches)
		coor_surface = [0] * len(branches)
		positions = [0] * len(branches)
		
		for i in xrange(len(branches)):
			# .swc coordinates of each branch
			coor_swc[i] = map(float, neurites[int(branches[i])][2:5])
			# image coordinates of each branch
			coor_image[i] = [coor_swc[i][0] / scale_default[0], coor_swc[i][1] / scale_default[1], coor_swc[i][2] / scale_default[2] + 1]
			# coordinates relative to center
			coor_relative[i] = [coor_image[i][j] - center[j] for j in xrange(3)]
			# coordinates converted to microns
			coor_microns[i] = [coor_relative[i][j] * scale_actual[j] for j in xrange(3)]
			# coordinates adjusted to surface (point of origin)
			multiplier = z_radius / sqrt(pow(coor_microns[i][0], 2) + pow(coor_microns[i][1], 2) + pow(coor_microns[i][2], 2))
			coor_surface[i] = [coor_microns[i][j] * multiplier for j in xrange(3)]
			# position analysis
			positions[i] = ceil(coor_surface[i][2]/((float(slice_h - slice_l) / num_pos) * scale_actual[2]) + (num_pos / 2))
		# corrects positions out of bound
		for p in positions:
			if p < 1:
				p = 1
			if p > num_pos:
				p = num_pos

		# Print Tests
		''' 	
		# print .swc coordinates
		print '\n.swc Coordinates (pre-calculated um) '.ljust(5*width, '-')
		for i in xrange(len(branches)):
			print 'Branch ' + branches[i] + ':\t' + '\t'.join(map(str, coor_swc[i]))
		# print image coordinates
		print '\nImage Coordinates (px, px, slice) '.ljust(5*width, '-')
		for i in xrange(len(branches)):
			print 'Branch ' + branches[i] + ':\t' + '\t'.join(map(str, coor_image[i]))			
		# print relative coordinates
		print '\nCoordinates Relative to Center (px, px, slice) '.ljust(5*width, '-')
		for i in xrange(len(branches)):
			print 'Branch ' + branches[i] + ':\t' + '\t'.join(map(str, coor_relative[i]))			
		# print micron coordinates
		print '\nCoordinates Relative to Center (um) '.ljust(5*width, '-')
		for i in xrange(len(branches)):
			print 'Branch ' + branches[i] + ':\t' + '\t'.join(map(str, coor_microns[i]))			
		# print surface coordinates
		print '\nSurface Coordinates (Points of Origin) (um) '.ljust(5*width, '-')
		for i in xrange(len(branches)):
			print 'Branch ' + branches[i] + ':\t' + '\t'.join(map(str, coor_surface[i]))			
				
		print ('\nBranch #'.ljust(width) + 'Dendrites'.rjust(width) + 'Axons'.rjust(width) + 'Position\n'.rjust(width) + ''.ljust((4*width), '-') + '\n')
		for i in xrange(len(branches)):
			print (('Branch ' + branches[i]).ljust(width) + str(dendrites[i]).rjust(width) + str(axons[i]).rjust(width) + str(int(positions[i])).rjust(width) + '\n')
		
		# information on end points		
		print '\nEndpoint Coordinates '.ljust(5*width, '-')
		
		for e in ends:
			print ('Endpoint ' + e + ':').ljust(width) + '\t' + '\t'.join(map(str, neurites[int(e)][2:5]))		
		'''
		endpoints = [0] * len(ends)
		end_um = [0] * len(ends)
		end_dist = [0] * len(ends)
		
		for i in xrange(len(ends)):
			# .swc coordinates of each point
			endpoints[i] = map(float, neurites[int(ends[i])][2:5])
			# image coordinates of each point
			endpoints[i] = [endpoints[i][0] / scale_default[0], endpoints[i][1] / scale_default[1], endpoints[i][2] / scale_default[2] + 1]
			# coordinates relative to center in microns
			end_um[i] = [((endpoints[i][j] - center[j]) * scale_actual[j]) for j in xrange(3)]
			# distance from center
			end_dist[i] = sqrt(pow(end_um[i][0], 2) + pow(end_um[i][1], 2) + pow(end_um[i][2], 2))
		
		# 	Print Test
		'''
		for i in xrange(len(endpoints)):
			print endpoints[i]
			print end_um[i]
			print end_dist[i]
			print 'minus xy-r: ' + str(end_dist[i] - (float(neurites[1][5])/scale_default[0]) * scale_actual[0])
			print 'minus z-r: ' + str(end_dist[i] - z_radius)
		'''
		
		# completed output written to file
		#
		# output to copy to Excel
		output.write('Copy to Excel: \n')
		for i in xrange(len(branches)):
			output.write(branches[i] + '\t' + str(dendrites[i]) + '\t' + str(axons[i]) + '\t' + str(int(positions[i])) + '\n')

		# output for viewing purposes
		#
		# endpoint coordinates
		output.write('\nEndpoint Coordinates '.ljust(6*width, '-'))
		output.write('\nEndpoint'.ljust(width) + '\t\t' + 'x (px)'.ljust(width) + 'y (px)'.ljust(width) + 'z (sl)\t' + 'um from surface (xy | z - radius)'.ljust(width))
		for i in xrange(len(ends)):		
			output.write(('\nEndpoint ' + ends[i] + ':').ljust(width) + '\t' + '\t'.join(map(str, endpoints[i])) + '\t' + str(end_dist[i] - (float(neurites[1][5])/scale_default[0]) * scale_actual[0]) + '\t' + str(end_dist[i] - z_radius))	
					
		# position analysis table
		output.write('\n\n Branch #'.ljust(width) + 'Dendrites'.rjust(width) + 'Axons'.rjust(width) + 'Position\n'.rjust(width) + ''.ljust((4*width), '-') + '\n')
		for i in xrange(len(branches)):
			output.write(('Branch ' + branches[i]).ljust(width) + str(dendrites[i]).rjust(width) + str(axons[i]).rjust(width) + str(int(positions[i])).rjust(width) + '\n')
		output.write('\nNumber of Ends:\t\t' + str(len(ends)))
		output.write('\nNumber of Branches:\t' + str(len(branches)))
		
		# information on cell body center and end points		
		output.write('\n\nCell body (px, px, sl):'.ljust(3*width) + str(center))
		output.write('\nMidslice:'.ljust(3*width) + str(midslice))
		output.write('\nPos height (slices):'.ljust(3*width) + str(float(slice_h - slice_l) / num_pos))
		output.write('\nPos height (um):'.ljust(3*width) + str((float(slice_h - slice_l) / num_pos) * scale_actual[2]))
		output.write('\nZ-radius (um):'.ljust(3*width) + str(z_radius))
		
		# coordinates
		# .swc coordinates
		output.write('\n\n.swc Coordinates (pre-calculated um) '.ljust(6*width, '-'))
		for i in xrange(len(branches)):
			output.write(('\nBranch ' + branches[i] + ':').ljust(width) + '\t\t' + '\t\t'.join(map(str, coor_swc[i])))
		# image coordinates
		output.write('\n\nImage Coordinates (px, px, slice) '.ljust(6*width, '-'))
		for i in xrange(len(branches)):
			output.write(('\nBranch ' + branches[i] + ':').ljust(width) + '\t\t' + '\t\t'.join(map(str, coor_image[i])))		
		# relative coordinates
		output.write('\n\nCoordinates Relative to Center (px, px, slice) '.ljust(5*width, '-'))
		for i in xrange(len(branches)):
			output.write(('\nBranch ' + branches[i] + ':').ljust(width) + '\t\t' + '\t\t'.join(map(str, coor_relative[i])))	
		# micron coordinates
		output.write('\n\nCoordinates Relative to Center (um) '.ljust(6*width, '-'))
		for i in xrange(len(branches)):
			output.write(('\nBranch ' + branches[i] + ':').ljust(width) + '\t\t' + '\t\t'.join(map(str, coor_microns[i]))	)
		# surface coordinates
		output.write('\n\nSurface Coordinates (Points of Origin) (um) '.ljust(6*width, '-'))
		for i in xrange(len(branches)):
			output.write(('\nBranch ' + branches[i] + ':').ljust(width) + '\t\t' + '\t\t'.join(map(str, coor_surface[i]))	)		

		# open file for viewing in notepad++
		# os.system("start notepad++ " + output_name)
		
		end = time.clock()
		print "%.2gs" % (end-start)
		
	# traces ends back to plotted point of origin
	def traceback(self, n):
		return neurites[n][6]
	
if __name__ == '__main__':
	# read parameters
	# ['PA.py', input, x_d, y_d, z_d, x_s, y_s, z_s, s_l, s_h]
	# input:			.swc file to be read; results will be output to the filename plus the extension '_output.txt'
	# scale_default:	refers to the scaling used when plotting the neurites in NeuronStudio
	# scale_actual:		refers to the actual scaling of the images
	# s_l:				lowest slice in which the cell body appears
	# s_h:				highest slice in which the cell body appears

	# Running through command line
	# ['PA.py', input, x_d, y_d, z_d, x_s, y_s, z_s, s_l, s_h]
	# PA = PositionAnalysis(sys.argv)
	
	# Running with manual input
	# Sample run
	PA = PositionAnalysis(['P10.py', 'sample.swc', 0.142, 0.142, 0.5, 0.142, 0.142, 0.5, 15, 31])
	
	# counts dendrites and axons for each branch
	# performs position analysis
	# outputs results to output file
	PA.work()
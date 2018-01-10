#
#	NESTED SAMPLING - THOM PIJNENBURG
#
#	Solution is spread over this program and some more comments and notes in the hw6.nb mathematica notebook,
#	which I also attached as pdf. The output of the program is saved in the output.txt doc, and figures are attached as svg's
#

from math import *
import pylab as pl
from random import *
import sys
import os
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

# This is the main part of the program
def main():
	global priorXPoints
	global priorYPoints
	global livePoints
	global lowestLikelihoods
	global highestPriorMass
	global M
	global Evidence
	global f
	global posterior

	# Define here the number of samples to be scattered across the parameter space. 
	M=3000;

	# Initialise the lists and variables needed to keep track of certain quantities
	priorXPoints = [];
	priorYPoints = [];
	livePoints = [];
	lowestLikelihoods = [];
	runs = [];
	Evidence=0;
	evidenceList = [];
	lHoodList = [];
	priorMassList = [];
	posterior = []
	
	j=0;
	alpha = 10**(-4)

	# Make a text file to save output
	f = open('output.txt','w')

	# Call the function that generates the Live Points, drawn from the priors
	generateLivePoints()

	plt.figure('Live point distribution. Before nested sampling.')
	xs = [x[0] for x in livePoints]
	ys = [x[1] for x in livePoints]
	plt.plot(xs,ys, marker=',', color= 'g', linestyle='None')
	save("LivePointsBefore", ext="svg", close=True, verbose=True)

	# Calculate the initial prior mass (this will be close to 1)
	highestPriorMass = priorMass(livePoints[0][0],livePoints[0][1])
	Lmax = livePoints[-1][2]

	# This is where we start picking up our point with lowest likelihood,
	# and replace it with a new insertion in our parameter space with higher likelihood
	print('Start nested sampling routine:')
	f.write('Start nested sampling routine: \n\n')

	for i in range(200000):
		lowestPoint = livePoints[0]
		lowestLikelihoods.append(lowestPoint[2])

		# Calculate the contribution of the point with lowest likelihood to the evidence integral
		contribution= computeAddedEvidence(lowestPoint,highestPriorMass) 
		Evidence += contribution[0]
		highestPriorMass= contribution[1]

		# Now replace the point by a new one with higher likelihood
		replacePointLowestLikelihood(lowestPoint)

		# This part is for terminal feedback about the progress of the program.
		# Every 1000 loops it prints the values of the smallest likelihood, highest prior mass,
		# and evidence
		Lmax = livePoints[-1][2]
		j+=1
		if j%1000==0:
			runs.append(j)
			evidenceList.append(Evidence)
			lHoodList.append(lowestPoint[2])
			priorMassList.append(highestPriorMass)

			print("Loop nr %d, smallest likelihood= %.15f, highest prior mass= %.15f, evidence= %.15f" %(j, lowestPoint[2], highestPriorMass, Evidence))
			f.write("Loop nr %d, smallest likelihood= %.15f, highest prior mass= %.15f, evidence= %.15f \n" %(j, lowestPoint[2], highestPriorMass, Evidence))
			
			pass

		# Termination condition: break out of loop when contribution to integral is too small
		if Lmax*highestPriorMass < alpha*Evidence:
			print("Accuracy reached! Computation terminated after %d iterations." % j)
			f.write("Accuracy reached! Computation terminated after %d iterations.\n\n" % j)

			break
		pass

	# Display value of the integral in the terminal
	print("Evidence integral = %.10f" % Evidence)
	f.write("Evidence integral = %.10f \n" % Evidence)

	# Close output file
	f.close()

	# Make figures of evolution of evidence, likelihood and prior mass.
	plt.figure('Evolution evidence, likelihood, priormass')
	plt.subplot(311)
	plt.scatter(runs,evidenceList)
	plt.ylabel('Evidence')
	plt.xlabel('# runs')

	plt.subplot(312)
	plt.scatter(runs, lHoodList,marker='h',c='r')
	plt.ylabel('Smallest Likelihood')
	plt.xlabel('# runs')

	plt.subplot(313)
	plt.scatter(runs, priorMassList,marker='s',c='y')
	plt.ylabel('highest prior mass')
	plt.xlabel('# runs')
	save("Evolution_evidence_likelihood_priormass", ext="svg", close=True, verbose=True)

	plt.figure('Live point distribution. After nested sampling.')
	xs = [x[0] for x in livePoints]
	ys = [x[1] for x in livePoints]
	plt.plot(xs,ys, marker=',', color= 'm', linestyle='None')
	save("LivePointsAfter", ext="svg", close=True, verbose=True)

	posterior.sort(key=lambda x:x[2])
	xp = [x[0] for x in posterior]
	yp = [x[1] for x in posterior]
	zp=[]
	for x in posterior:
		zp.append(x[2]/Evidence*x[3])
		pass

	plt.subplots()
	plt.scatter(xp, yp, c=zp, s=1, edgecolor='')
	plt.colorbar()
	save("Posterior_Density", ext="svg", close=True, verbose=True)
	# Make figures
	plt.show()

# End of program
######################################
# define functions:

# Prior
def priorX(x) :
	return -2.0*(x-0.5)

# Integrated prior
def cdfX(x) :
	return 0.75+x-x**2

# Rewrite CDF value into actual parameter value
def cdfX_to_X(x) :
	return 0.5*(1-2*(1-x)**0.5)

# Prior
def priorY(y) :
	return 2*(y+0.5)

# Integrated prior
def cdfY(y) :
	return 0.25+y+y**2

# Rewrite CDF value into actual parameter value
def cdfY_to_Y(y) :
	return 0.5*(-1+2*(y)**(0.5))

# Calculate prior mass. Used for first point
def priorMass(x , y):
	def integrand(x, y):
	    return priorX(x)*priorY(y)

	ans, err = dblquad(integrand, -fabs(x),fabs(x) ,
		                   lambda x: -fabs(y),
		                   lambda x: fabs(y))
	return ans

# Likelihood function
def likelihoodL(x , y):
	sigmaX= 0.35;
	sigmaY=0.5;
	return 1/(2*pi*sigmaX*sigmaY)*exp(-(x**2/(2*sigmaX**2))-(y**2/(2*sigmaY**2)))

# Generate M live points drawn from priors
def generateLivePoints():

	# First draw points from prior according distribution
	for i in range(M):
		x = cdfX_to_X(uniform(0,1))
		y = cdfY_to_Y(uniform(0,1))
		priorXPoints.append(x)
		priorYPoints.append(y)
		lHood = likelihoodL(x,y)
		livePoints.append([x,y,lHood])
		pass

	# Sort the list on likelihood values
	livePoints.sort(key=lambda x:x[2])
	print("Generated %d live points" % M)
	f.write("Generated %d live points\n\n" % M)
	pass

# Draw a new point from the prior
def generateNewPoint():
	newpoint = []
	xPoint = cdfX_to_X(uniform(0,1))
	yPoint = cdfY_to_Y(uniform(0,1))
	lHood = likelihoodL(xPoint,yPoint)
	newpoint = [xPoint,yPoint,lHood]
	return newpoint

# Compute the contribution of the point with lowest likelihood to the evidence integral
def computeAddedEvidence(item,priorMass):

	# Statistically assign the new lowest prior mass by drawing from t distribution
	t = (uniform(0.,1.))**(1./(M))
	highestPriorMass=t*priorMass
	deltaX=(1-t)*priorMass

	# Gather point with likelihood value and delta X to construct the posterior density later
	posteriorItem = item
	posteriorItem.append(deltaX)
	posterior.append(posteriorItem)
	return [item[2]*deltaX,highestPriorMass]

# Replace the point by a new one drawn from the prior
def replacePointLowestLikelihood(item):
	replaced = False
	count = 0
	# Keep trying until a new point has been found
	while replaced == False:
		newPoint = generateNewPoint()

		# Keep new point if the likelihood is higher than the old point
		if newPoint[2] > item[2]:
			livePoints.remove(item)
			livePoints.append(newPoint)
			livePoints.sort(key=lambda x:x[2])
			replaced = True
			pass
		count +=1
		if count >1000000:
			print('No new point found in %d tries' % count)
			break
		pass
	pass

# Function to save figures to local directory
def save(path, ext='png', close=True, verbose=True):
	# Extract the directory and filename from the given path
	directory = os.path.split(path)[0]
	filename = "%s.%s" % (os.path.split(path)[1], ext)
	if directory == '':
		directory = '.'
	# If the directory does not exist, create it
	if not os.path.exists(directory):
		os.makedirs(directory)
	# The final path to save to
	savepath = os.path.join(directory, filename)
	if verbose:
		print("Saving figure to '%s'..." % savepath),
		# Actually save the figure
	plt.savefig(savepath)
		# Close it
	if close:
		plt.close()
	if verbose:
		print("Done")

if __name__ == '__main__':
    sys.exit(main())

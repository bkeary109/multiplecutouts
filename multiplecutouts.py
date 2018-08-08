import numpy as np
import scipy
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import sys

### Error catching ###
if len(sys.argv) != 3:
    sys.exit("Failure. Run with start and end file numbers.")

# Fit funtion
def exp1(x, a, b, c):
    return (a*np.exp(-x/b) + c)

####### file name here ######
inputfileroot = "20180808_"
startfile = int(sys.argv[1])
endfile = int(sys.argv[2])
backgroundfile = "background.asc"
time_offset = 0. # ns

# cutouts
cutout1 = 735 #pixelnum 585 nm
cutout2 = 1039 #pixelnum 606 nm
cutout3 = 1593 #1359pixelnum 628 nm

#plotpoint (pixel number)
plotpoint = 1024


for q in xrange (startfile, endfile+1):
    print q

    # File names
    qstr = str(q)
    inputfile = inputfileroot + qstr + ".asc"
    outputfile = inputfileroot + qstr + "_lifetimesback.txt"
    outputspecfile = inputfileroot + qstr + "_cleanback.txt"
    cutoutfile = inputfileroot + qstr + "_cutout.txt"

    dark = 0.

    # Find num_steps and stepsize
    fullfile = open(inputfile, 'r')
    all_lines = fullfile.readlines()
    steps_line = all_lines[-22]
    num_steps = int("".join([steps_line[i] for i in xrange(-6,-1)]))
    size_line = all_lines[-10]
    stepsize = int("".join([size_line[i] for i in xrange(-6,-1)]))
    fullfile.close()

    rawspectra = np.genfromtxt(inputfile, delimiter=',', skip_footer = 28)
    background = np.genfromtxt(backgroundfile, delimiter=',', skip_footer = 28)

    #Initialise arrays
    times = []
    points = [0 for j in xrange(num_steps)]
    testpoints = [0 for j in xrange(num_steps)]
    wavelengths = []
    lifetimes = [0 for j in xrange(2047)]
    res_vals = [0 for j in xrange(2047)]
    backgrounddata = [0 for j in xrange(2047)]
    spectra = [[0 for j in xrange(num_steps)] for i in xrange(2047)]
    roughspec = [[0 for j in xrange(num_steps)] for i in xrange(2047)]

    #Output arrays
    lifedata = np.zeros(shape=(2047,3))
    outputspectra = np.zeros(shape=(2047,num_steps))
    cutoutpoints1 = np.zeros(shape=(num_steps, 1))
    cutoutpoints2 = np.zeros(shape=(num_steps, 1))
    cutoutpoints3 = np.zeros(shape=(num_steps, 1))
    cutoutdata = np.zeros(shape=(num_steps, 4))

    for i in range (1, num_steps+1):
        times.append(int(stepsize*(i-1) + time_offset))

    times = np.array(times)

    for i in range(0, 2047):
        wavelengths.append(rawspectra[i,0])
        backgrounddata[i] = background[i,1]

    # Remove background and smooth
    for i in xrange(1, num_steps+1):
        for j in xrange(0, 2047):
            newval = rawspectra[j,i] - (backgrounddata[j])
            if newval < 0:
                newval = 1

            roughspec[j][i-1] = newval

        for j in xrange(3, 2043):
            avgspec = [0 for a in range(6)]
            avgspec[0] = roughspec[j-3][i-1]
            avgspec[1] = roughspec[j-2][i-1]
            avgspec[2] = roughspec[j-1][i-1]
            avgspec[3] = roughspec[j+1][i-1]
            avgspec[4] = roughspec[j+2][i-1]
            avgspec[5] = roughspec[j+3][i-1]

            averg = np.mean(avgspec)

            if roughspec[j][i-1] > (averg*1.01):
                roughspec[j][i-1] = averg
            elif roughspec[j][i-1] < (averg/1.01):
                roughspec[j][i-1] = averg

            spectra[j][i-1] = roughspec[j][i-1]


    # Test plot used to determine fitting limits
    for t in xrange(0, num_steps):
        testpoints[t] = spectra[plotpoint][t]

    plt.semilogy(times, testpoints, 'o')
    plt.xlabel(q)
    plt.show()
    print "Enter fitting limits (ns): "
    sys.stdout.flush()
    str_start, str_end = raw_input().split()
    fit_start = int(np.floor(int(str_start)/stepsize))
    fit_end = int(np.floor(int(str_end)/stepsize))

    # Determine lifetimes
    for l in xrange(350, 1700):
        for n in xrange(0, num_steps):
            points[n] = spectra[l][n]

        if l == cutout1:
            for p in range(0, num_steps):
                cutoutpoints1[p] = points[p]
        elif l == cutout2:
            for t in range(0, num_steps):
                cutoutpoints2[t] = points[t]
        elif l == cutout3:
            for r in range(0, num_steps):
                cutoutpoints3[r] = points[r]

        popt, pcov = curve_fit(exp1, times[fit_start:fit_end], points[fit_start:fit_end], p0=(10000, 100, 1000), maxfev=5000)
        lifet = popt[1]
        lifetimes[l] = lifet

        resid = points[fit_start:fit_end] - exp1(times[fit_start:fit_end], *popt)
        ss_res = np.sum(resid**2)
        ss_tot = np.sum((points[fit_start:fit_end] - np.mean(points[fit_start:fit_end]))**2)
        r_squared = 1 - (ss_res/ss_tot)
        res_vals[l] = r_squared

	# Print lifetime at l = cutout1
	if l == cutout1:
	    print("\n\nPixel "+str(cutout1))
	    print("Cutout1 lifetime = "+str(lifet))
	    print("r^2 = "+str(r_squared))

	# Print lifetime at l = cutout2
	if l == cutout2:
	    print("\n\nPixel "+str(cutout2))
	    print("Cutout2 lifetime = "+str(lifet))
	    print("r^2 = "+str(r_squared))

        # Print lifetime at l = cutout3
        if l == cutout3:
	    print("\n\nPixel "+str(cutout3))
	    print("Cutout3 lifetime = "+str(lifet))
	    print("r^2 = "+str(r_squared))

	# Plot sample at l = plotpoint
	if (l == plotpoint):
            x_new = np.linspace(times[fit_start], times[fit_end], 500)
            y_new = exp1(x_new, *popt)
            plt.plot(times, points, 'o')
            plt.plot(x_new, y_new)
	    plt.title(r'$\tau = ' + str(int(lifet)) + ' ns , r^{2} = $' + str(r_squared))
	    plt.xlabel(q)
            plt.show()


    # Smooth lifetimes
    for a in xrange(340,1690):
	avglif = [0 for i in range (6)]
	avglif[0] = lifetimes[a-3]
	avglif[1] = lifetimes[a-2]
	avglif[2] = lifetimes[a-1]
	avglif[3] = lifetimes[a+1]
	avglif[4] = lifetimes[a+2]
	avglif[5] = lifetimes[a+3]
	lifnow = lifetimes[a]


	averglif = np.mean(avglif)

	if lifnow > (averglif*1.05):
	    lifnow = averglif
	elif lifnow < (averglif/1.05):
	    lifnow = averglif

	lifetimes[a] = lifnow

    # Create arrays for output
    for i in xrange(0,2047):
        lifedata[i][0] = wavelengths[i]
        lifedata[i][1] = lifetimes[i]
        lifedata[i][2] = res_vals[i]

        for j in xrange(0, num_steps-1):
            outputspectra[i][0] = wavelengths[i]
            outputspectra[i][j+1] = spectra[i][j]

    plt.plot(wavelengths, lifetimes)
    plt.xlabel(q)
    plt.show()

    for k in xrange(0, num_steps):
        cutoutdata[k][0] = times[k]
        cutoutdata[k][1] = cutoutpoints1[k]
        cutoutdata[k][2] = cutoutpoints2[k]
        cutoutdata[k][3] = cutoutpoints3[k]


    print("Save data? (Y/N)")
    sys.stdout.flush()
    save = raw_input()

    if save == "Y" or save == "y":
        np.savetxt(outputspecfile, outputspectra, delimiter=',')
        np.savetxt(outputfile, lifedata, delimiter=',')
        np.savetxt(cutoutfile, cutoutdata, delimiter=',')
	print("\n\nLifetime and time-profile data saved.\n\n")
    else:
	print("\n\nLifetime and time-profile data not saved.\n\n")


import utility
import argparse
import numpy as np
from scipy.signal import correlate

parser = argparse.ArgumentParser()
parser.add_argument("filename",help="Path to the file to open")
parser.add_argument("-first_sigma",help="Pulse width of primary pulse",default=None,type=float)
parser.add_argument("-second_sigma",help="Pulse width of secondary pulse",default=None,type=float)
parser.add_argument("-the_delay",help="Time delay of the secondary pulse",default=None,type=float)
parser.add_argument("-the_amplitude",help="Amplitude of the secondary pulse",default=None,type=float)
args = parser.parse_args()

#Data: 16 channels, 512 samples per channel, 1.5 GHz sampling rate
fs = 1.5
trace_length = 512
n_channels = 7
data = utility.load_csv_data(open(args.filename,newline=''),trace_length,n_channels)
csw = utility.get_csw(data,2,[0,1,2,3,4,5,6])
times = np.arange(0,trace_length)/fs
f0 = 0.19
gamma = 0.02
secondary_amplitude = 0.01
plot_delay = 106
model = np.zeros(trace_length)

if(args.first_sigma and args.second_sigma and args.the_delay and args.the_amplitude):
	for i in range(trace_length):
		model[i] = utility.math_conv(times[i]-plot_delay,1,1,args.first_sigma,f0,gamma)
		model[i] += utility.math_conv(times[i]-args.the_delay-plot_delay,1,1,args.second_sigma,f0,gamma)*args.the_amplitude
	model = model/np.sqrt(np.inner(model,model))
	print(np.max(np.abs(correlate(model,csw))))
	for i in range(trace_length):
		print(times[i],csw[i],model[i])
else:
	rho_max = 0
	sigmas = np.arange(0.3,3.0,0.1)
	best_sigma_main = 0
	for sigma_1 in sigmas:
		for i in range(trace_length):
			model[i] = utility.math_conv(times[i],1,1,sigma_1,f0,gamma)
		model = model/np.sqrt(np.inner(model,model))
		rho = np.max(np.abs(correlate(model,csw)))
		if(rho>rho_max):
			best_sigma_main=sigma_1
			rho_max = rho
	print(rho_max,best_sigma_main)
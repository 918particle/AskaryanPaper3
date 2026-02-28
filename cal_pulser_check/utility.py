import csv
import numpy
import argparse
import numpy as np
from scipy.signal import convolve,hilbert,butter,filtfilt
from scipy.special import dawsn as D
from scipy.special import wofz as w
from scipy.integrate import quad
from scipy.signal import correlate,correlation_lags

def shifted_laplace_transform(f,s,t0,t_max=1000):
    real_integral = quad(lambda t: np.exp(-s.real*t)*np.cos(s.imag*t)*f(t-t0),0,t_max)[0]
    imag_integral = quad(lambda t: np.exp(-s.real*t)*np.sin(s.imag*t)*f(t-t0),0,t_max)[0]
    return complex(real_integral,-imag_integral)
def dwdq(q):
	J = complex(0,1)
	return 2*J/np.sqrt(np.pi)-2*q*w(q)
def s(t,E0,sigma_t):
	return -E0*t*np.exp(-0.5*t*t/sigma_t/sigma_t)
def r(t,R0,f0,gamma):
	if(t>=0):
		return R0*np.cos(2*np.pi*f0*t)*np.exp(-2*np.pi*gamma*t)
	else:
		return 0
def math_conv(t,E0,R0,sigma_t,f0,gamma):
    J = complex(0,1)
    units = R0*E0*sigma_t*sigma_t
    x = t/(sigma_t*np.sqrt(2))
    z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
    q = -J*(x+0.5*z)
    Re_wq = np.real(w(q))
    Re_dwdx = np.real(-J*dwdq(q))
    result = -np.sqrt(np.pi)*units*(x*np.exp(-x*x)*Re_wq-0.5*np.exp(-x*x)*Re_dwdx)
    if(np.isinf(result) or np.isnan(result)):
        return 0.0
    else:
        return result
def math_env(t,E0,R0,sigma_t,f0,gamma):
	units = R0*E0*sigma_t*sigma_t
	J = complex(0,1)
	x = t/(sigma_t*np.sqrt(2))
	z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
	k=-z
	q = -J*(x+z/2)
	first_term = -np.sqrt(np.pi)*units*(x*np.exp(-x*x)*w(q)+0.5*J*np.exp(-x*x)*dwdq(q))
	second_term = 2/np.sqrt(np.pi)*units*(D(x) + k*shifted_laplace_transform(D,k,x))
	if(np.isinf(first_term) or np.isinf(second_term) or np.isnan(first_term) or np.isnan(second_term)):
		return 0.0
	else:
		result = first_term+J*second_term
	return 0.5*np.abs(result)
def get_csw(trace_set,index,good_channels):
	ref = trace_set[index]
	n = len(ref)
	sum_trace = np.copy(ref)
	lag_array = correlation_lags(n,n,mode='full')
	for i, trace in enumerate(trace_set):
		if i in good_channels:
			trace -= np.mean(trace)
			corr = correlate(ref,trace,mode='full')/n
			best_lag = lag_array[np.argmax(corr)]
			aligned_trace = np.roll(trace,best_lag)
			sum_trace += aligned_trace
	norm = np.sqrt(np.inner(sum_trace,sum_trace))
	sum_trace/=norm
	return sum_trace
def lowpass_butter(x,cutoff,fs,order=4):
    nyq = 0.5*fs
    b,a = butter(order,cutoff/nyq,btype='low')
    return filtfilt(b,a,x)
def load_csv_data(csvfile,trace_length,n_channels):
	data = np.zeros((n_channels,trace_length))
	reader = csv.reader(csvfile,delimiter=',')
	i=0
	for row in reader:
		data_row = list(map(float,row))
		data_row = lowpass_butter(data_row,0.35,1.5)
		data[i] = data_row/np.sqrt(np.inner(data_row,data_row))
		i+=1
		if(i==n_channels):
			break
	return data
import numpy as np
from scipy import signal

class Lockin:
	def __init__(self, sigfunc, fs, BW=100, freq=1e3, order=4):
		self.sigfunc = sigfunc # function that create the signal
		self.freq		 = freq # Lockin frequency
		self.LP_BW		= BW # LP filter bandwidth
		self.LP_order = order	# Order of the low pass-filter (ie. 6dB/oct per order)
		self.fs			 = fs # Sampling rate

	def getSignal(self):
		return self.sigfunc(self.freq)

	def getXY(self, freq=None):
		if freq is None: freq = self.freq
		sig = self.getSignal()
		self.sig = sig
		t=np.linspace(0,len(sig)/self.fs,len(sig))
		self.t = t
		nyq_rate = self.fs / 2.0
		cutoff = self.LP_BW/nyq_rate
		b, a = signal.butter(self.LP_order, cutoff, 'low', analog=False)
		sX = 2*signal.lfilter(b,a,sig*np.sin(2*np.pi*t*freq))
		sY = 2*signal.lfilter(b,a,sig*np.cos(2*np.pi*t*freq))
		return sX,sY

	def getRP(self, freq=None):
		if freq is None: freq = self.freq
		X,Y = self.getXY(freq)
		R = np.sqrt(X**2+Y**2)
		P = np.arctan2(Y,X)
		return R,P

	def sweep(self, _from, _to=None, step=None, N=1000, TC=300):
		if _to is None:
			assert type(_from)==list or type(_from)==tuple or type(_from)==np.ndarray
			freqs = _from
		elif step in None:
			freqs = np.linspace(_from, _to, N)
		else:
			freq = np.arange(_from, _to, step)
		Dat = {'R':[],'P':[],'X':[],'Y':[]}
		for f in freqs:
			X,Y = self.getXY(f)
			R = np.sqrt(X**2+Y**2)
			P = np.arctan2(Y,X)
			i = np.argmin(abs(self.t-TC*self.LP_BW/self.fs))
			Dat['X'].append(X[i])
			Dat['Y'].append(Y[i])
			Dat['R'].append(R[i])
			Dat['P'].append(P[i])
		for i in 'XYRP': Dat[i]=np.array(Dat[i])
		return Dat

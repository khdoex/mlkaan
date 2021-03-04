# **Paper Section for Featurization**
### Contents
#### 1- Abstract
#### 2- Continous Wavelet Transform
#### 3- Implementation on Seismic Data
#### 4- Experiments
#### 5- Results 


### 1-Abstract

Purpose of Featurization experiments is developing the EQTransformer model by adding some features to it. For this purpose, we focused on the CWT(continuous Wavelet Transform) for an analysis of the seismic trace by getting the time-frequency representation of the trace.
After we implied the CWT using Obspy, we separated the frequencies from each other and added them to the original stream as an external channel. After that, we used the new stream as input in the earthquake detection algorithms Eqtransformer and the classic LTA/STA.  
 We believe that with the usage of CWT, we can get more information about the seismic trace and therefore detect earthquakes in a more efficient way. 


### 2- Continous Wavelet Transform
Continuous Wavelet Transform is a tool for getting the time-frequency based reconstruction of the original signal by using a chosen mother wavelet. CWT reconstructs the signal by letting the scale and translation parameters of the wavelets vary continuously.

#### 2.1 Mathematically
Wavelet Transform of $$x$$$$\big(t\big)$$ at position $$u$$ and scale $$s$$ is equal to:   
##### $$X$$ $$_{w}$$ $$\big(u,s\big)$$ =     $$\int_{-\infty}^\infty$$ $$\frac{1}{|s|^{1/2}}$$ $$x$$ $$\big(t\big)$$ $$\bar{\psi}$$ $$\big(\frac{t-u}{s})$$ $$dt.$$

Where $$\psi$$ is the wavelet function, and $$\bar{\psi}$$ is complex conjugate of $$\psi$$. 
We worked on the Morlet Wavelet for analyzing the seismic signal, which is:
##### $$\psi$$  $$\big(t\big)$$ =  $$\pi^{-1/4}$$ $$e^{i w_{0}  t}$$ $$e^{-t^{2}/2 }$$
here $$t$$ is time and  $$w_{0}$$ is wavenumber of the wavelet.

#### 2.2 Code
For implying this transformation in an computationally efficient way, it is used with Fast Fourier Transformation.
from the [Obspy Cwt function](<https://docs.obspy.org/_modules/obspy/signal/tf_misfit.html>)


    `def cwt(st, dt, w0, fmin, fmax, nf=100, wl='morlet'):
    Continuous Wavelet Transformation in the Frequency Domain.

        .. seealso:: [Kristekova2006]_, eq. (4)

        :param st: time dependent signal.
        :param dt: time step between two samples in st (in seconds)
        :param w0: parameter for the wavelet, tradeoff between time and frequency
        resolution
        :param fmin: minimum frequency (in Hz)
        :param fmax: maximum frequency (in Hz)
        :param nf: number of logarithmically spaced frequencies between fmin and
        fmax
        :param wl: wavelet to use, for now only 'morlet' is implemented

        :return: time frequency representation of st, type numpy.ndarray of complex
        values, shape = (nf, len(st)).
    
        npts = len(st) * 2
        tmax = (npts - 1) * dt
        t = np.linspace(0., tmax, npts)
        f = np.logspace(np.log10(fmin), np.log10(fmax), nf)

        cwt = np.zeros((npts // 2, nf), dtype=np.complex)

        if wl == 'morlet':

            def psi(t):
                return np.pi ** (-.25) * np.exp(1j * w0 * t) * \
                np.exp(-t ** 2 / 2.)

            def scale(f):
                return w0 / (2 * np.pi * f)
        else:
            raise ValueError('wavelet type "' + wl + '" not defined!')

        nfft = util.next_pow_2(npts) * 2
        sf = np.fft.fft(st, n=nfft)

        # Ignore underflows.
        with np.errstate(under="ignore"):
            for n, _f in enumerate(f):
                a = scale(_f)
                # time shift necessary, because wavelet is defined around t = 0
                psih = psi(-1 * (t - t[-1] / 2.) / a).conjugate() / np.abs(a) ** .5
                psihf = np.fft.fft(psih, n=nfft)
                tminin = int(t[-1] / 2. / (t[1] - t[0]))
                cwt[:, n] = np.fft.ifft(psihf * sf)[tminin:tminin + npts // 2] * \
                (t[1] - t[0])

        return cwt.T`

### 3- Implementation on Seismic Data

"seismic signal image"

We take each channel, one by one for the cwt.
then implying the cwt with chosen parameters 
dt (time step between two samples in seconds) = 1/100                                

below parameters are chosen according to experiments, this is only one example.

w0 (wavenumber) =5
fmin (minimum freq. in Hz) = 10 
fmax (maximum freq. in Hz) = 40 
nf (number of logarithmically spaced frequencies between fmin and fmax) = 4 


After we imply the cwt on the seismic trace, it will compute an array of shape (nf, len(trace).

Each new channel will be look like these

"cwt output image"




# **Paper Section for Featurization**
### Contents
#### 1- Abstract
#### 2- Continous Wavelet Transform     
#####  2.1-Mathematically                     
#####   2.2-Code                
#### 3- Implementation on Seismic Data
##### 3.1 CWT On Noise
##### 3.2 CWT On Earthquake
#### 4- Experiments and Results
##### 4.1- STA/LTA
##### 4.2- EQTransformer
###### 4.2.1- Dataset
###### 4.2.2- Results
#
#
### 1-Abstract

Purpose of Featurization experiments is developing the EQTransformer model by adding some features to it. For this purpose, we focused on the CWT(continuous Wavelet Transform) for an analysis of the seismic trace by getting the time-frequency representation of the trace.
After we implied the CWT using Obspy, we separated the frequencies from each other and added them to the original stream as an external channel. After that, we used the new stream as input in the earthquake detection algorithms Eqtransformer and the classic STA/LTA.  
 We believe that with the usage of CWT, we can get more information about the seismic trace and therefore detect earthquakes in a more efficient way. 


### 2- Continous Wavelet Transform
Continuous Wavelet Transform is a tool for getting the time-frequency based reconstruction of the original signal by using a chosen mother wavelet. CWT reconstructs the signal by letting the scale and translation parameters of the wavelets vary continuously.

#### 2.1 Mathematically
Wavelet Transform of x at position u and scale s is equal to:

![form 1_b.png](https://www.dropbox.com/s/9hwiedspq1ofplk/form%201_b.png?dl=0&raw=1)
<!--- $$X$$ $$_{w}$$ $$\big(u,s\big)$$ =     $$\int_{-\infty}^\infty$$ $$\frac{1}{|s|^{1/2}}$$ $$x$$ $$\big(t\big)$$ $$\bar{\psi}$$ $$\big(\frac{t-u}{s})$$ $$dt.$$ --->

Where psi is the wavelet function, and psi bar is complex conjugate of psi. 
We worked on the Morlet Wavelet for analyzing the seismic signal, which is:

![form2 24.png](https://www.dropbox.com/s/n7xnmmqm9p10dda/form2%2024.png?dl=0&raw=1)
<!--- $$\psi$$  $$\big(t\big)$$ =  $$\pi^{-1/4}$$ $$e^{i w_{0}  t}$$ $$e^{-t^{2}/2 }$$ --->
here t is time and  w0 is wavenumber of the wavelet.

#### 2.2 Code
For applying this transformation in an computationally efficient way, it is used with Fast Fourier Transformation.
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
#### 3.1- CWT On Noise
As an example, we have a trace without an earthquake. We will apply CWT on this data with chosen parameters. 
![Normal_Trace1.png](https://www.dropbox.com/s/thesnhr3tkinzuq/Normal_Trace1.png?dl=0&raw=1)

We choose parameters for CWT as follows:
`dt` (time step between two samples in seconds) = 1/100                                
below parameters are chosen according to experiments, this is only one example.
`w0` (wavenumber) =5
`fmin` (minimum freq. in Hz) = 10 
`fmax` (maximum freq. in Hz) = 40 
`nf` (number of logarithmically spaced frequencies between fmin and fmax) = 2 

Here the nf is the number of frequency channels, If we choose nf=1, then the only frequency channel has frequency='fmin'. If we choose nf=2, then we choose nf=2, then we have two frequency channel with frequencies 'fmin' and 'fmax'. Furthermore, if we choose nf>2, then the  frequencies between 'fmin' and 'fmax' will be calculated spaced evenly on a log scale (as an example, with fmin=10 and fmax=40, nf=3 will give [10,20,40] as frequencies).


After we apply the cwt on the seismic trace, it will compute an array of shape `(nf, len(trace))`.

Each new channel will be look like these:

![CWT_0_normal1.png](https://www.dropbox.com/s/g5t1n0t5s584vnb/CWT_0_normal1.png?dl=0&raw=1)![CWT_1_normal1.png](https://www.dropbox.com/s/slwu8x9wxbw66fk/CWT_1_normal1.png?dl=0&raw=1)

#### 3.2- CWT on Earthquake

![ANormal_Trace.png](https://www.dropbox.com/s/5xzot49vn2mqdak/ANormal_Trace.png?dl=0&raw=1)

We choose parameters for CWT same with the noise case:
`dt` (time step between two samples in seconds) = 1/100                                
below parameters are chosen according to experiments, this is only one example.
`w0` (wavenumber) =5
`fmin` (minimum freq. in Hz) = 10 
`fmax` (maximum freq. in Hz) = 40 
`nf` (number of logarithmically spaced frequencies between fmin and fmax) = 2 
Each new channel will be look like these:

![CWT_0_Anormal.png](https://www.dropbox.com/s/ep995i7ar23zdzd/CWT_0_Anormal.png?dl=0&raw=1)![CWT_1_Anormal.png](https://www.dropbox.com/s/931j4mi1k8gtnpq/CWT_1_Anormal.png?dl=0&raw=1)

After we apply the cwt, we can observe that the signal can be analyzed using cwt.

### 4- Experiments
#### 4.1- STA/LTA

Experimenting STA/LTA algorithm on the transformed data is also can be useful.
We used Obspy's Classic STA/LTA function for this purpose.
After we apply the STA/LTA algorithm on both transformed and original signal and plot them with the original signal as follows:
#####  STA/LTA on original signal and transformed signal
We Plot the whole procedure we applied on the raw noise data. First plot is the raw data, second plot is the regular sta/lta, third and fourth plots are the sta/lta applied on CWT data. 
![feature_4pl.jpg](https://www.dropbox.com/s/rlhodvzdr5vgqug/feature_4pl.jpg?dl=0&raw=1)
Red line on the STA/LTA plots shows how big is the LTA window. 
Same setting of plot with the earthquake data.
![feature_4pl_e.jpg](https://www.dropbox.com/s/siwgeusj4v6qq0c/feature_4pl_e.jpg?dl=0&raw=1)
Red line on the STA/LTA plots shows how big is the LTA window. 

#### 4.2- EQTransformer
EQTransformer is a multi-task deep neural network for simultaneous earthquake detection and phase picking with a hierarchical attentive model. It mainly consists of one very deep encoder and three separate decoders (detector, P-picker, and S-picker branches) with an attention mechanism.  

EQTransformer model emulates this through two levels of attention mechanism in a hierarchical structure. one at the global level for identifying an earthquake signal in the input time series, and one at the local level for identifying different seismic phases within that earthquake signal. Two levels of self-attention (global and local) help the model capture and exploit dependencies between local (individual phases) and global (full-waveform) features within an earthquake signal.

![eqt_model.png](https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41467-020-17591-w/MediaObjects/41467_2020_17591_Fig1_HTML.png?as=webp)

We modified EQTransformer to take more than 3 channels. Our extra channels contains CWT version of the original data. Number of extra channels differs according to the `nf` parameter of CWT. 

The aim of the Featurization Experiments was the developing the EQTransformer model using CWT. We added extra channels to the original signal and gave that signal to the model.  Model become much more heavy with CWT. This is the disadvantage of using CWT. 

However we are still experimenting the parameters, by changing them systematically.

##### 4.2.1 Dataset
The Dataset that we are using is from STEAD dataset, Stead-Mini and Stead Micro. Stead-Mini consists 23542 noise and 103023 earthquakes. Stead-Micro consists 2354 noise and 10302 earthquakes.

##### 4.2.2 Results
The Results of these experiments can be found in the table below. EQT is the original model of EQTransformer.
If not specified in the name, Obspy CWT function is used.

| | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|model_name|Dataset|Frequency Range|Number Of Frequency | wavelet number|det_recall|det_precision|d_tp|d_fp|d_tn|d_fn|p_recall|p_precision|p_mae|p_rmse|p_tp|p_fp|p_tn|p_fn|s_recall|s_precision|s_mae|s_rmse|s_tp|s_fp|s_tn|s_fn|#events|#noise|
|EQT|mini|0|0|0|0.999416966281217|0.997139948616026|20570|59|4672|12|0.998299484986882|0.997717781878217|1.856441648|14.7431762|20547|47|4684|35|0.997959381984258|0.997426309911135|3.507730266|60.46299362|20540|53|4678|42|20582|4731|
|CWT|micro|10-40|5|5|0.999076863278593|0.996848943184022|20563|65|4666|19|0.997424934408707|0.997909780283881|1.85705513|14.74175762|20529|43|4688|53|0.994461179671558|0.997271487039563|3.422096949|45.70946891|20468|56|4675|114|20582|4731|
|CWT|micro|10-40|5|5|0.99758337361044|0.996138996138996|2064|8|455|5|0.993716771387144|0.998058252427184|2.737354322|25.24684706|2056|4|459|13|0.982600289995167|0.99656862745098|3.885558616|70.24723164|2033|7|456|36|2069|463|
|CWT|mini|30-35|3|5|0.999611310854144|0.996078431372549|20574|81|4650|8|0.998882518705665|0.997331910352188|1.928166708|15.97408395|20559|55|4676|23|0.993829559809542|0.99668664425279|3.293230533|30.32674047|20455|68|4663|127|20582|4731|
|CWT|mini|10-40|4|5|0.99975706928384|0.995789779326365|20577|87|4644|5|0.998445243416577|0.997088791848617|2.025456656|15.77075064|20550|60|4671|32|0.997473520551939|0.996118389131489|3.421787582|46.46852805|20530|80|4651|52|20582|4731|
|CWT|mini|10-40|5|5|0.99951413856768|0.99309678976587|20572|143|4588|10|0.996890486833155|0.995632763975155|2.04719615|19.0691514|20518|90|4641|64|0.99003984063745|0.994922122943216|3.371328266|33.91297117|20377|104|4627|205|20582|4731|
|CWT|micro|10-35|2|5|0.998547215496368|0.991346153846154|2062|18|449|3|0.989346246973366|0.996099463676256|3.054173889|30.01596897|2043|8|459|22|0.990314769975787|0.991274842462433|4.937576802|131.5476544|2045|18|449|20|2065|467|
|CWT|micro|10-35|3|5|0.998572108519752|0.990557129367327|2098|20|411|3|0.993336506425512|0.995231282784931|2.736538867|26.9977962|2087|10|421|14|0.986673012851023|0.992340832934418|3.84651915|60.2700139|2073|16|415|28|2101|431|
|CWT|micro|10-40|3|3|0.998548621190131|0.988979396262578|2064|23|442|3|0.99225931301403|0.993220338983051|2.87214452|21.92395412|2051|14|451|16|0.995645863570392|0.990375360923965|3.91592073|64.1083559|2058|20|445|9|2067|465|
|CWT|micro|30-35|3|5|0.998538723818802|0.988904968644477|2050|23|456|3|0.993667803214808|0.990772219524041|2.808254932|23.40365176|2040|19|460|13|0.987335606429615|0.990229604298974|3.978056581|63.96795571|2027|20|459|26|2053|479|
|CWT|micro|10-40|2|8|0.998044009779951|0.988856589147287|2041|23|464|4|0.990220048899756|0.990220048899756|2.683512162|14.21422561|2025|20|467|20|0.990220048899756|0.991189427312775|3.773251469|41.68546877|2025|18|469|20|2045|487|
|CWT|micro|1-15|5|5|0.998524348253812|0.988796882610813|2030|23|476|3|0.994589276930644|0.993611793611793|2.724206426|15.87625451|2022|13|486|11|0.988686669945893|0.990147783251232|4.062582389|112.6021531|2010|20|479|23|2033|499|
|CWT|micro|30-50|3|5|0.997582205029014|0.98755385351843|2063|26|438|5|0.995164410058027|0.991807228915662|2.801344464|28.40949359|2058|17|447|10|0.987911025145067|0.988867376573088|4.255854982|94.78550226|2043|23|441|25|2068|464|
|ScipyMorletCWT|micro|1-50|1|2|0.997575169738118|0.987518002880461|2057|26|444|5|0.989330746847721|0.992217898832685|3.321146405|28.48948607|2040|16|454|22|0.99078564500485|0.989825581395349|4.0647004|53.84479108|2043|21|449|19|2062|470|
|ScipyMorletCWT|micro|10-35|5|5|0.995145631067961|0.98700048146365|2050|27|445|10|0.992718446601942|0.989356555394291|2.597813626|15.34745072|2045|22|450|15|0.978155339805825|0.990658800393314|3.670178567|32.76312167|2015|19|453|45|2060|472|
|PywtMexhCWT|micro|10-40|3|not specified|0.997138769670958|0.985855728429986|2091|30|405|6|0.985693848354793|0.988995215311005|2.993208283|24.12684651|2067|23|412|30|0.988555078683834|0.989026717557252|4.056493081|64.6020487|2073|23|412|24|2097|435|
|ScipyMorletCWT|micro|1-50|1|5|0.996593673965937|0.985563041385948|2048|30|447|7|0.992700729927007|0.987415295256534|2.85737971|15.96826142|2040|26|451|15|0.983941605839416|0.990205680705191|4.19492386|90.26341892|2022|20|457|33|2055|477|
|CWT|micro|30-35|5|5|0.998042094958394|0.984548527281506|2039|32|457|4|0.986784140969163|0.988235294117647|2.658979598|21.27885924|2016|24|465|27|0.978952520802741|0.989119683481701|3.866437103|36.73935038|2000|22|467|43|2043|489|
|ScipyRickerCWT|micro|1-50|3|not specified|0.995611896635787|0.983622350674374|2042|34|447|9|0.977571916138469|0.991592482690406|3.104106242|23.03609123|2005|17|464|46|0.984397854705022|0.986803519061584|3.747948251|35.9180189|2019|27|454|32|2051|481|
|CWT|micro|1-15|3|5|0.999521759923482|0.982604607428303|2090|37|404|1|0.995695839311334|0.988603988603989|3.077628151|19.53554585|2082|24|417|9|0.99091343854615|0.986197049024274|3.915714788|36.68606498|2072|29|412|19|2091|441|
|PywtGaussianDerivativeCWT|micro|1-50|3|not specified|0.995633187772926|0.982288176160842|2052|37|434|9|0.991266375545852|0.989346246973366|3.430691635|27.18983548|2043|22|449|18|0.983017952450267|0.98636806231743|4.047528974|73.7800269|2026|28|443|35|2061|471|
|PywtGaussianDerivativeCWT|micro|1-50|3|not specified|0.995658465991317|0.975886524822695|2064|51|408|9|0.986975397973951|0.985074626865672|3.53121218|31.37822568|2046|31|428|27|0.977809937288953|0.980648282535075|4.561624592|77.53995919|2027|40|419|46|2073|459|
|ScipyMorletCWT|micro|1-50|1|10|0.997577519379845|0.974443918599148|2059|54|414|5|0.994186046511628|0.982758620689655|3.195907807|36.36210722|2052|36|432|12|0.986434108527132|0.977436389822372|4.648558201|114.7453662|2036|47|421|28|2064|468|
|ScipyMorletCWT|micro|10-30|3|5|0.999513381995134|0.973921289710763|2054|55|422|1|0.994160583941606|0.977979894686453|3.064786555|21.77954126|2043|46|431|12|0.992700729927007|0.976543800861656|4.246188365|91.26571409|2040|49|428|15|2055|477|
|CWT|micro|30-35|3|5|0.999520383693046|0.97337692666978|2084|57|390|1|0.997601918465228|0.980669495520981|3.370963269|19.20225126|2080|41|406|5|0.989928057553957|0.981921979067555|4.807955523|134.033205|2064|38|409|21|2085|447|
|ScipyMorletCWT|micro|10-35|6|5|0.996122152205526|0.969339622641509|2055|65|404|8|0.989335918565196|0.979836773883821|3.825506045|39.5589133|2041|42|427|22|0.983519146873485|0.974075852136342|4.912571272|117.7313624|2029|54|415|34|2063|469|
|PywtMorletCWT|micro|10-35|3|5|0.99803536345776|0.96900333810205|2032|65|431|4|0.985756385068762|0.982859941234084|3.742294885|35.06960267|2007|35|461|29|0.984282907662082|0.971871968962173|5.087535143|118.8629969|2004|58|438|32|2036|496|



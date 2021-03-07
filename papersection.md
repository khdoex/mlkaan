# **Paper Section for Featurization**
### Contents
#### 1- Abstract
#### 2- Continous Wavelet Transform
#### 3- Implementation on Seismic Data
#### 4- Experiments and Results


### 1-Abstract

Purpose of Featurization experiments is developing the EQTransformer model by adding some features to it. For this purpose, we focused on the CWT(continuous Wavelet Transform) for an analysis of the seismic trace by getting the time-frequency representation of the trace.
After we implied the CWT using Obspy, we separated the frequencies from each other and added them to the original stream as an external channel. After that, we used the new stream as input in the earthquake detection algorithms Eqtransformer and the classic STA/LTA.  
 We believe that with the usage of CWT, we can get more information about the seismic trace and therefore detect earthquakes in a more efficient way. 


### 2- Continous Wavelet Transform
Continuous Wavelet Transform is a tool for getting the time-frequency based reconstruction of the original signal by using a chosen mother wavelet. CWT reconstructs the signal by letting the scale and translation parameters of the wavelets vary continuously.

#### 2.1 Mathematically
Wavelet Transform of x at position u and scale s is equal to:

![form 1_b.png](https://www.dropbox.com/s/9hwiedspq1ofplk/form%201_b.png?dl=0&raw=1)

Where psi is the wavelet function, and psi bar is complex conjugate of psi. 
We worked on the Morlet Wavelet for analyzing the seismic signal, which is:

![form2 24.png](https://www.dropbox.com/s/n7xnmmqm9p10dda/form2%2024.png?dl=0&raw=1)

here t is time and  w0 is wavenumber of the wavelet.

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
#### 3.1- CWT On Noise
![Normal_Trace.png](https://www.dropbox.com/s/954khku48ql11t7/Normal_Trace.png?dl=0&raw=1)

We take each channel, one by one for the cwt.
then implying the cwt with chosen parameters 
`dt` (time step between two samples in seconds) = 1/100                                

below parameters are chosen according to experiments, this is only one example.

`w0` (wavenumber) =5
`fmin` (minimum freq. in Hz) = 10 
`fmax` (maximum freq. in Hz) = 40 
`nf` (number of logarithmically spaced frequencies between fmin and fmax) = 4 


After we imply the cwt on the seismic trace, it will compute an array of shape (nf, len(trace).

Each new channel will be look like these:

![CWT_0_normal.png](https://www.dropbox.com/s/bvtn2kow9a00pe5/CWT_0_normal.png?dl=0&raw=1)![CWT_1_normal.png](https://www.dropbox.com/s/uvl078gk53fiunv/CWT_1_normal.png?dl=0&raw=1)

#### 3.2- CWT on Earthquake

![ANormal_Trace.png](https://www.dropbox.com/s/5xzot49vn2mqdak/ANormal_Trace.png?dl=0&raw=1)

Simply, we imply the same procedure with same parameters.
Each new channel will be look like these:

![CWT_0_Anormal.png](https://www.dropbox.com/s/ep995i7ar23zdzd/CWT_0_Anormal.png?dl=0&raw=1)![CWT_1_Anormal.png](https://www.dropbox.com/s/931j4mi1k8gtnpq/CWT_1_Anormal.png?dl=0&raw=1)

After we imply the cwt, we can observe that the signal can be analyzed better using cwt.

### 4- Experiments
#### 4.1- STA/LTA

Experimenting STA/LTA algorithm on the transformed data is also can be useful.
We used Obspy's Classic STA/LTA function for this purpose.

After we imply the STA/LTA algorithm on both transformed and original signal and plot them as follows:
#####  STA/LTA on original signal
On Noise:

![stalta_n.png](https://www.dropbox.com/s/escqenb5vpnio0q/stalta_n.png?dl=0&raw=1)

On Earthquake:

![stalta_an.png](https://www.dropbox.com/s/abjffxp5vlrqt4k/stalta_an.png?dl=0&raw=1)

#### STA/LTA on transformed signal

On Noise:

![staltacw_n.png](https://www.dropbox.com/s/sd1iozcnzi6wy4c/staltacw_n.png?dl=0&raw=1) ![staltacw1_n.png](https://www.dropbox.com/s/fn7ltps6tflunc0/staltacw1_n.png?dl=0&raw=1)

On Earthquake:

![staltacw_an.png](https://www.dropbox.com/s/i70ux1tlb1b8rz4/staltacw_an.png?dl=0&raw=1) ![staltacw1_an.png](https://www.dropbox.com/s/5h1fl2qpx22etpj/staltacw1_an.png?dl=0&raw=1)
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
The Results of these experiments can be found in the table below. Starting with EQT_CWT_micro, till the end are the experiments of CWT.
If not specified, every model(that uses cwt) is trained with micro dataset. Models that are trained using Mini dataset are have mini extension on its name.

Here the nf is number of logarithmically spaced frequencies between fmin and fmax. Values after fr(frequency range) are fmin and fmax. And w is wavenumber of the wavelet.


| | | | | | | | | | | | | | | | | | | | | | | | | |
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|model_name|det_recall|det_precision|d_tp|d_fp|d_tn|d_fn|p_recall|p_precision|p_mae|p_rmse|p_tp|p_fp|p_tn|p_fn|s_recall|s_precision|s_mae|s_rmse|s_tp|s_fp|s_tn|s_fn|#events|#noise|
|EQT|0.999416966281216|0.997139948616026|20570|59|4672|12|0.998299484986881|0.997717781878216|1.856441648|14.7431762|20547|47|4684|35|0.997959381984258|0.997426309911134|3.507730266|60.46299362|20540|53|4678|42|20582|4731|
|Salvation|0.999659896997376|0.993769319938176|20575|129|4602|7|0.996307453114371|0.996259048729534|1.889314985|15.67784302|20506|77|4654|76|0.99203187250996|0.995611468695143|3.364538025|45.34527868|20418|90|4641|164|20582|4731|
|Genesis|0.999659896997376|0.99636803874092|20575|75|4656|7|0.998153726557186|0.997475237910273|2.173377654|17.85295549|20544|52|4679|38|0.995238557963268|0.997127975466095|3.283805775|29.63703477|20484|59|4672|98|20582|4731|
|Vanilla|0.99927120785152|0.996994522274468|20567|62|4669|15|0.991157321931785|0.998678220002937|2.216988|15.73635624|20400|27|4704|182|0.985618501603342|0.997541306058222|4.021425353|96.14685711|20286|50|4681|296|20582|4731|
|GRU|0.999854241570304|0.994202618483984|20579|120|4611|3|0.999368380137984|0.996029247978306|1.987080846|21.58908472|20569|82|4649|13|0.99577300553882|0.995096135171878|3.29897954|37.35907761|20495|101|4630|87|20582|4731|
|4x4|0.999465552424448|0.993096456502848|20571|143|4588|11|0.998736760275969|0.994966118102613|2.138280463|19.70681076|20556|104|4627|26|0.993586629093382|0.993876360808709|3.344848962|38.35297944|20450|126|4605|132|20582|4731|
|2LSTM|0.999805655427072|0.996609841146842|20578|70|4661|4|0.999028277135361|0.997380675203725|1.846317987|13.60885379|20562|54|4677|20|0.997424934408706|0.997134252962891|3.346229979|40.51779752|20529|59|4672|53|20582|4731|
|EQUtils|0.998547215496368|0.978178368121442|2062|46|421|3|0.995641646489104|0.987986544930322|2.905260656|21.40837119|2056|25|442|9|0.994673123486682|0.979961832061068|4.355389912|83.45000417|2054|42|425|11|2065|467|
|Kernel|0.991279069767441|0.987928536938677|2046|25|443|18|0.984980620155038|0.99267578125|2.985078649|15.70033522|2033|15|453|31|0.979651162790697|0.991662579695929|4.109833951|47.85803984|2022|17|451|42|2064|468|
|NB_Filter|0.998547918683446|0.991350312349831|2063|18|448|3|0.994675701839303|0.992274263640753|2.572405846|15.34952154|2055|16|450|11|0.994675701839303|0.992753623188405|4.656873053|146.2052428|2055|15|451|11|2066|466|
|EQT Yunus|0.991185112634672|0.98252427184466|2024|36|454|18|0.972086190009794|0.99002493765586|2.670646439|30.29544371|1985|20|470|57|0.980901077375122|0.984759095378564|4.143172387|80.25946589|2003|31|459|39|2042|490|
|BC|0.992410827219833|0.933380918391625|3923|280|4476|30|0.971414115861371|0.961201501877346|3.503028596|25.09710017|3840|155|4601|113|0.935745003794586|0.942420382165605|6.118806287|136.3850341|3699|226|4530|254|3953|4756|
|BCLOS|0.993992490613266|0.974717722140402|3971|103|4611|24|0.989486858573216|0.982355864811133|3.035391847|21.89737163|3953|71|4643|42|0.98072590738423|0.97754491017964|3.811012655|39.66749113|3918|90|4624|77|3995|4714|
|BO|0.990950226244344|0.926657263751763|3942|312|4419|36|0.975615887380593|0.939937030758052|3.259230738|21.32297176|3881|248|4483|97|0.97611865258924|0.940421409542262|5.077923143|125.0135624|3883|246|4485|95|3978|4731|
|EQT Ecem|0.988585607940446|0.954480114997604|3984|190|4489|46|0.977667493796526|0.960038986354775|3.760692786|34.08360059|3940|164|4515|90|0.978908188585608|0.970479704797048|5.646353251|155.121097|3945|120|4559|85|4030|4679|
|Partial Vanilla Ecem|0.998230088495575|0.905920146856356|3948|410|4344|7|0.98811630847029|0.9260663507109|2.654157541|26.30038175|3908|312|4442|47|0.882932996207332|0.939214631522324|4.394066201|110.1605632|3492|226|4528|463|3955|4754|
|EQT_CWT_micro|0.999520383693045|0.97337692666978|2084|57|390|1|0.997601918465227|0.98066949552098|3.370963269|19.20225126|2080|41|406|5|0.989928057553956|0.981921979067554|4.807955523|134.033205|2064|38|409|21|2085|447|
|EQT_CWT|0.999611310854144|0.996078431372549|20574|81|4650|8|0.998882518705665|0.997331910352187|1.928166708|15.97408395|20559|55|4676|23|0.993829559809542|0.996686644252789|3.293230533|30.32674047|20455|68|4663|127|20582|4731|
|EQT_CWT_NF4|0.99975706928384|0.995789779326364|20577|87|4644|5|0.998445243416577|0.997088791848617|2.025456656|15.77075064|20550|60|4671|32|0.997473520551938|0.996118389131489|3.421787582|46.46852805|20530|80|4651|52|20582|4731|
|EQT_CWT_NF5|0.999076863278592|0.996848943184021|20563|65|4666|19|0.997424934408706|0.997909780283881|1.85705513|14.74175762|20529|43|4688|53|0.994461179671557|0.997271487039563|3.422096949|45.70946891|20468|56|4675|114|20582|4731|
|nf-3-fr30-35-obspycwt w=5|0.998538723818801|0.988904968644476|2050|23|456|3|0.993667803214807|0.99077221952404|2.808254932|23.40365176|2040|19|460|13|0.987335606429615|0.990229604298974|3.978056581|63.96795571|2027|20|459|26|2053|479|
|nf-1-fr-30-35-obspycwt w=5|0.998042094958394|0.984548527281506|2039|32|457|4|0.986784140969163|0.988235294117647|2.658979598|21.27885924|2016|24|465|27|0.978952520802741|0.989119683481701|3.866437103|36.73935038|2000|22|467|43|2043|489|
|nf-5-fr-10-40-obspycwt w=5|0.997583373610439|0.996138996138996|2064|8|455|5|0.993716771387143|0.998058252427184|2.737354322|25.24684706|2056|4|459|13|0.982600289995166|0.99656862745098|3.885558616|70.24723164|2033|7|456|36|2069|463|
|nf-5-fr-1-15-obspycwt w=5|0.998524348253812|0.988796882610813|2030|23|476|3|0.994589276930644|0.993611793611793|2.724206426|15.87625451|2022|13|486|11|0.988686669945892|0.990147783251231|4.062582389|112.6021531|2010|20|479|23|2033|499|
|nf-3-fr-1-15-obspycwt w=5|0.999521759923481|0.982604607428302|2090|37|404|1|0.995695839311334|0.988603988603988|3.077628151|19.53554585|2082|24|417|9|0.99091343854615|0.986197049024274|3.915714788|36.68606498|2072|29|412|19|2091|441|
|nf-3-fr-30-50-obspycwt w=5|0.997582205029013|0.987553853518429|2063|26|438|5|0.995164410058027|0.991807228915662|2.801344464|28.40949359|2058|17|447|10|0.987911025145067|0.988867376573088|4.255854982|94.78550226|2043|23|441|25|2068|464|
|nf-3-fr-1-50-scipymorlet w=5|0.996593673965936|0.985563041385948|2048|30|447|7|0.992700729927007|0.987415295256534|2.85737971|15.96826142|2040|26|451|15|0.983941605839416|0.990205680705191|4.19492386|90.26341892|2022|20|457|33|2055|477|
|nf-3-fr-10-30-scipymorlet w=5 |0.999513381995133|0.973921289710763|2054|55|422|1|0.994160583941606|0.977979894686452|3.064786555|21.77954126|2043|46|431|12|0.992700729927007|0.976543800861656|4.246188365|91.26571409|2040|49|428|15|2055|477|
|nf-3-fr-1-50-scipyricker|0.995611896635787|0.983622350674373|2042|34|447|9|0.977571916138469|0.991592482690405|3.104106242|23.03609123|2005|17|464|46|0.984397854705022|0.986803519061583|3.747948251|35.9180189|2019|27|454|32|2051|481|
|nf-1-fr-1-50-scipymorlet w=2|0.997575169738118|0.98751800288046|2057|26|444|5|0.98933074684772|0.992217898832684|3.321146405|28.48948607|2040|16|454|22|0.990785645004849|0.989825581395348|4.0647004|53.84479108|2043|21|449|19|2062|470|
|nf-1-fr-1-50-scipymorlet w=10|0.997577519379845|0.974443918599148|2059|54|414|5|0.994186046511628|0.982758620689655|3.195907807|36.36210722|2052|36|432|12|0.986434108527131|0.977436389822371|4.648558201|114.7453662|2036|47|421|28|2064|468|
|nf-5-fr-10-35-scipymorlet w=5|0.995145631067961|0.987000481463649|2050|27|445|10|0.992718446601941|0.989356555394291|2.597813626|15.34745072|2045|22|450|15|0.978155339805825|0.990658800393313|3.670178567|32.76312167|2015|19|453|45|2060|472|
|nf-6-fr-10-35-scipymorlet w=5|0.996122152205526|0.969339622641509|2055|65|404|8|0.989335918565196|0.979836773883821|3.825506045|39.5589133|2041|42|427|22|0.983519146873485|0.974075852136341|4.912571272|117.7313624|2029|54|415|34|2063|469|
|nf-3-fr-10-35-obspymorlet w=5|0.998572108519752|0.990557129367327|2098|20|411|3|0.993336506425511|0.99523128278493|2.736538867|26.9977962|2087|10|421|14|0.986673012851023|0.992340832934418|3.84651915|60.2700139|2073|16|415|28|2101|431|
|nf-2-fr-10-35-obspymorlet w=5|0.998547215496368|0.991346153846154|2062|18|449|3|0.989346246973365|0.996099463676255|3.054173889|30.01596897|2043|8|459|22|0.990314769975786|0.991274842462433|4.937576802|131.5476544|2045|18|449|20|2065|467|
|nf-3-fr-1-50-gaussianderivative2 |0.995658465991316|0.975886524822695|2064|51|408|9|0.98697539797395|0.985074626865671|3.53121218|31.37822568|2046|31|428|27|0.977809937288953|0.980648282535075|4.561624592|77.53995919|2027|40|419|46|2073|459|
|nf-3-fr-1-50-gaussianderivative1 |0.995633187772925|0.982288176160842|2052|37|434|9|0.991266375545851|0.989346246973365|3.430691635|27.18983548|2043|22|449|18|0.983017952450266|0.986368062317429|4.047528974|73.7800269|2026|28|443|35|2061|471|
|nf-3-fr-10-35-pywtmorlet w=5|0.99803536345776|0.96900333810205|2032|65|431|4|0.985756385068762|0.982859941234084|3.742294885|35.06960267|2007|35|461|29|0.984282907662082|0.971871968962172|5.087535143|118.8629969|2004|58|438|32|2036|496|
|nf-3 fr-10-40- pywtmexh|0.997138769670958|0.985855728429985|2091|30|405|6|0.985693848354792|0.988995215311004|2.993208283|24.12684651|2067|23|412|30|0.988555078683834|0.989026717557252|4.056493081|64.6020487|2073|23|412|24|2097|435|
|nf-2-fr-10-40-obspycwt w=8|0.998044009779951|0.988856589147286|2041|23|464|4|0.990220048899755|0.990220048899755|2.683512162|14.21422561|2025|20|467|20|0.990220048899755|0.991189427312775|3.773251469|41.68546877|2025|18|469|20|2045|487|
|nf-3-fr-10-40-obspycwt w=3|0.99854862119013|0.988979396262578|2064|23|442|3|0.99225931301403|0.99322033898305|2.87214452|21.92395412|2051|14|451|16|0.995645863570392|0.990375360923965|3.91592073|64.1083559|2058|20|445|9|2067|465|
|nf-2-fr-10-40 obspycwt w=5 mini|0.99951413856768|0.99309678976587|20572|143|4588|10|0.996890486833155|0.995632763975155|2.04719615|19.0691514|20518|90|4641|64|0.99003984063745|0.994922122943215|3.371328266|33.91297117|20377|104|4627|205|20582|4731|

EQT_CWT_Micro = nf-3-fr-30-35-obspycwt-w=5 

EQT_CWT =  nf-3-fr-30-35-obspycwt-w=5 Mini version

EQT_CWT_NF4 = nf-4-fr-10-40-obspycwt-w=5 Mini version

EQT_CWT_NF5 = nf-5-fr-10-40-obspycwt-w=5 Mini version



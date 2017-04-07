# ASN
ambient seismic noise

noisez.f90 is a simple z-component noise simulation program, which can only deal with laterally homogeneous medias. Thus the only input file is a Rayleigh-wave phase velocity file, rayleighpv.dat. All sources are set to be Rick wavelet, and the default certral frequency is 10Hz. 24 sensors are arranged near the origin, and the interval is 5m. The noise source distribution default setting is even, inline and outline. The source strength can be changed with the source number or taper windows. The output files are in SAC format per sensor. 

asncor.sh is a bash script for cross-correlation computing. The program will call SAC (Seimic Analysis Code) and use cutsac.sh.

asncoh.sh is a script for cross-coherence computation. The program will use cutsac.sh and asncoh.f90. The output coherence possesses both lags, the positive time and the negative time one.

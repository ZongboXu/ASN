# ASN simulation and cross-correlation
ambient seismic noise

noisez.f90 is a Z-component ambient seismic data simulation program. This program is only for laterally homogeneous medium and only requires a Rayleigh-wave phase velocity file as the input, rayleighpv.dat. All seismic sources are set to be Rick wavelet, and the default certral frequency is 10Hz. 24 sensors are arranged near the origin, and the interval is 5m. The noise source distributions can be even, inline or outline. Users can change the source strengths by chaging the source number or taper windows. The output files are in SAC format per sensor. 

asncor.sh is a bash script for cross-correlation computing. The program will call SAC (Seimic Analysis Code) in cutsac.sh to cut seismic recordings.

asncoh.sh is a script for cross-coherence computation. The program will use cutsac.sh and asncoh.f90. The output coherence possesses both lags, the positive time and the negative time ones.

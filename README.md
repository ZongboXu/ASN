# Ambient seismic noise (ASN) simulation and cross-correlation
This is a gather of ASN codes which I wrote and used a lot. 

noisez.f90 is a Z-component ambient seismic data simulation program. This program is only for laterally homogeneous medium and only requires a Rayleigh-wave phase velocity file as the input, rayleighpv.dat. All seismic sources are set to be Rick wavelet, and the default certral frequency is 10Hz. 24 sensors are arranged near the origin, and the interval is 5m. The noise source distributions can be even, inline or outline. Users can change the source strengths by chaging the source number or taper windows. The output files for sensors are in SAC format. 

asncor.sh is a bash script for computing cross-correlations. The program calls SAC (Seimic Analysis Code) in cutsac.sh to cut seismic recordings and runs cross-correlation.

asncoh.sh is a bash script for computing cross-coherence (i.e. spectral normalized cross-correlation). The program calls cutsac.sh and asncoh.f90. The output coherence possesses both lags, the positive- and negative-time ones.

###References:
Xu, Z., Xia, J., Luo, Y., Cheng, F., & Pan, Y. (2016). [Potential misidentification of love-wave phase velocity based on three-component ambient seismic noise](https://doi.org/10.1007/s00024-015-1160-4). Pure and Applied Geophysics, 173(4), 1115-1124.

Xu, Z., Dylan Mikesell, T., Xia, J., & Cheng, F. (2017). [A comprehensive comparison between the refraction microtremor and seismic interferometry methods for phase-velocity estimation](https://doi.org/10.1190/geo2016-0654.1). Geophysics, 82(6), EN99-EN108.

Xu, Z., & Mikesell, T. D. (2017). [On the reliability of direct Rayleigh-wave estimation from multicomponent cross-correlations](https://doi.org/10.1093/gji/ggx228). Geophysical Journal International, 210(3), 1388-1393.

import pyRing.waveform as wf

m1=10
m2=20
chi1=0
chi2=0
phases={(2,2): 0}
r=8
iota=0
phi=0
modes=[(2,2)]
TGR_params={}


TEOB = wf.TEOBPM(0,m1,m2,chi1,chi2,phases,r,iota,phi,modes,TGR_params)

fit_coefficients = TEOB.fit_coefficients[(2,2)]

sigma_real = fit_coefficients['alpha1']
a1         = fit_coefficients['a1']
a2         = fit_coefficients['a2']
a3         = fit_coefficients['a3']
a4         = fit_coefficients['a4']

real_amp = TEOB.TEOBPM_Amplitude(0, 0, a1,a2,a3,a4)

print(real_amp)
print(TEOB.Mf)
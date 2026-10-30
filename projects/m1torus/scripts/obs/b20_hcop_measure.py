#!/usr/bin/env python3
"""B20: HCO+ minus CO centroid shift in the receding half of the line (NGC 613), the
measurement pre-registered in LOG 2026-10-30 (commit 181abbdc7); errors from noise
realisations at the HCO+ rms."""
import sys
sys.path[:0]=['/home/etrevol/phantom/projects/m1torus/scripts/obs','/home/etrevol/phantom/projects/m1torus/scripts/vis']
import numpy as np, mockobs as mo
from dataclasses import replace
from scipy.ndimage import gaussian_filter
from fits_min import read_fits
C=299792.458; VSYS=1481.0; EA,NA=0.0100,0.0025
def load(path, nu0):
    h,d=read_fits(path); c=d[0]; n=h['NAXIS3']; f=h['CRVAL3']+(np.arange(n)+1-h['CRPIX3'])*h['CDELT3']
    v=C*(nu0/f-1)-VSYS; o=np.argsort(v); pix=abs(h['CDELT2'])*3600
    east=-(np.arange(h['NAXIS1'])+1-h['CRPIX1'])*pix; north=(np.arange(h['NAXIS2'])+1-h['CRPIX2'])*pix
    return c[o],v[o],east,north,pix,np.sqrt(h['BMAJ']*h['BMIN'])*3600
def red_m0(c,v,east,north,pix,beam,lo,hi,free):
    dv=float(np.median(np.diff(v))); sb=beam/mo.FWHM/pix
    rms_s=float(np.std(gaussian_filter(c,(1.5,sb,sb),mode='wrap')[free]))
    t=replace(mo.NGC613,beam_arcsec=beam,dv=dv,pix_arcsec=pix,clip_mjy=0.0)
    k=(v>lo)&(v<hi)
    m0=mo.moments(dict(cube=c[k],vs=v[k],xs=north),t,rms=rms_s)[0]
    return np.where(np.isfinite(m0),m0,0)
def centroid(m0,east,north,R=0.1):
    E,N=np.meshgrid(east,north); w=np.where(np.hypot(E-EA,N-NA)<R,m0,0)
    return (w*E).sum()/w.sum()-EA,(w*N).sum()/w.sum()-NA, w.sum()
pcas=83.39
hc,hv,he,hn,hpix,hbeam=load('/home/etrevol/phantom/materials/data/alma_ngc613/ngc613_hcop43_spw21_cut035.fits',356.7342e9)
cc,cv,ce,cn,cpix,cbeam=load('/home/etrevol/phantom/materials/data/alma_ngc613/ngc613_co32_spw25_cut035.fits',345.79599e9)
print(f"HCO+ cube: v {hv.min():.0f}..{hv.max():.0f} km/s, dv {np.median(np.diff(hv)):.1f}, beam {hbeam:.4f}; rms (v>600) {np.std(hc[hv>600])*1e3:.3f} mJy")
for lo,hi in ((20,300),(20,230)):
    mh=red_m0(hc,hv,he,hn,hpix,hbeam,lo,hi,hv>600); mc=red_m0(cc,cv,ce,cn,cpix,cbeam,lo,hi,np.abs(cv)>450)
    ch=centroid(mh,he,hn); cco=centroid(mc,ce,cn)
    dE,dN=ch[0]-cco[0],ch[1]-cco[1]
    rng=np.random.default_rng(7); sb=hbeam/mo.FWHM/hpix; rms=np.std(hc[hv>600]); shifts=[]
    for i in range(40):
        nz=gaussian_filter(rng.standard_normal(hc.shape),(0,sb,sb),mode='wrap'); nz*=rms/nz.std()
        m=red_m0(hc+nz,hv,he,hn,hpix,hbeam,lo,hi,hv>600); c2=centroid(m,he,hn); shifts.append((c2[0]-ch[0],c2[1]-ch[1]))
    sE,sN=np.std(np.array(shifts),axis=0); sig=np.hypot(sE,sN)*pcas
    print(f"v {lo}..{hi}: CO red-half centroid {np.hypot(*cco[:2])*pcas:.2f} pc at PA {np.degrees(np.arctan2(*cco[:2]))%360:.0f}; "
          f"HCO+ {np.hypot(*ch[:2])*pcas:.2f} pc at PA {np.degrees(np.arctan2(*ch[:2]))%360:.0f}; "
          f"shift HCO+ - CO {np.hypot(dE,dN)*pcas:.2f} pc toward PA {np.degrees(np.arctan2(dE,dN))%360:.0f}, noise error {sig:.2f} pc ({np.hypot(dE,dN)*pcas/sig:.1f} sigma)")

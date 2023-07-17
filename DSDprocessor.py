'''
需要的环境：Fortran
    可参照https://zhuanlan.zhihu.com/p/76613134
使用的库：https://github.com/jleinonen/pytmatrix
Leinonen, J., High-level interface to T-matrix scattering calculations: architecture, capabilities and limitations, Opt. Express, vol. 22, issue 2, 1655-1660 (2014), doi: 10.1364/OE.22.001655.
'''

import numpy as np
import os
import pandas as pd
import datetime
import warnings
warnings.filterwarnings('ignore')


import pytmatrix
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive

diameter=np.array([0.062, 0.187, 0.312, 0.437, 0.562,
                      0.687, 0.812, 0.937, 1.062, 1.187, 
                      1.375, 1.625, 1.875, 2.125, 2.375, 
                      2.750, 3.250, 3.750, 4.250, 4.750,
                      5.500, 6.500, 7.500, 8.500, 9.500,
                      11.00, 13.00, 15.00, 17.00, 19.00, 
                      21.50, 24.50])# mm
velocity=9.65-10.3*np.exp(diameter*-0.6)# m/s
spread=np.array([0.125, 0.125, 0.125, 0.125, 0.125,
                    0.125, 0.125, 0.125, 0.125, 0.125, 
                    0.25, 0.25, 0.25, 0.25, 0.25,
                    0.5, 0.5, 0.5, 0.5, 0.5,
                    1, 1, 1, 1, 1, 
                    2, 2, 2, 2, 2,
                    3, 3])# mm

bins=np.hstack((0,diameter+np.array(spread)/2))
    
def thurai(D_eq):
  if D_eq < 0.7:
    return 1
  elif D_eq>=0.7 and D_eq <=1.5:
    return 1.173-0.5165*D_eq+0.4698*D_eq**(2)-0.1317*D_eq**3-8.5*10**(-3)*D_eq**4
  else:
    return 1.065-6.25*10**(-2)*D_eq-3.99*10**(-3)*D_eq**2+7.66*10**(-4)*D_eq**3-4.095*10**(-5)*D_eq**4

def prv_cal(nds, band):
    z=[]
    zdr=[]
    cc=[]
    kdp=[]
    ah=[]
    for i in range(len(nds)):
    # for i in range(1):
        nd=nds[i]
        # 设置：设置波段和雨滴温度
        scatterer=Scatterer(wavelength=tmatrix_aux.wl_X, m=refractive.m_w_20C[tmatrix_aux.wl_X])
       
        scatterer.psd_integrator=PSDIntegrator()
        # 设置2：雨滴形状模型，这里选择Thurai（2007）的雨滴形状模型
        scatterer.psd_integrator.axis_ratio_func=lambda D: 1.0/thurai(D)
        # 设置3：Dmax
        scatterer.psd_integrator.D_max=8
        scatterer.psd_integrator.geometries=(tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)
        # scatterer.or_pdf = orientation.gaussian_pdf(20.0)
        # scatterer.orient = orientation.orient_averaged_fixed
        scatterer.psd_integrator.init_scatter_table(scatterer)
        # scatterer.psd = GammaPSD(D0=2.0, Nw=1e3, mu=4)
        
        BinnedDSD=pytmatrix.psd.BinnedPSD(bins, nd)
        scatterer.psd=BinnedDSD
        scatterer.set_geometry(tmatrix_aux.geom_horiz_back)
        
        z.append(radar.refl(scatterer))
        zdr.append(radar.Zdr(scatterer))
        cc.append(radar.rho_hv(scatterer))
        scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)
        kdp.append(radar.Kdp(scatterer))
        ah.append(radar.Ai(scatterer))

        
   
    
    z=np.array(z)
    zdr=np.array(zdr)
    cc=np.array(cc)
    kdp=np.array(kdp)
    ah=np.array(ah)
    
    mask=np.where(z<=0) # log(z),z>0
    z[mask]=10**(-3.3)
    zdr[mask]=-33
    cc[mask]=-33
    kdp[mask]=-33
    ah[mask]=-33
    z=10*np.log10(z)
    
    
    prv=np.array([z,zdr,cc,kdp,ah]).transpose()
    df = pd.DataFrame(prv, columns = ['zh','zdr','cc','kdp','ah'])
    return df


def main(path):
    print(path)
    ls = os.listdir(path)
    path_save = os.path.join(path, 'prv')
    if os.path.exists(path_save) == False:
        os.mkdir(path_save)
        print('now making '+path_save)
    
    for f in ls[:]:
        if f.endswith('xlsx'):
            fpath = os.path.join(path, f)
            fpath_save = os.path.join(path_save, f)
            if os.path.exists(fpath_save) == False:
                
            
                df = pd.read_excel(fpath, index_col = 0)
                nds = df.iloc[:, -32:].values
                prvs = prv_cal(nds[:])
                dfnew = pd.concat([df, prvs], axis=1)
                dfnew.to_excel(fpath_save)
    print(path+' done')
    return 1
    # nds=np.genfromtxt("nds.csv",delimiter=",")[:10]
    # df=prv_cal(nds)
    # a=1
'''
test
'''
if __name__ == '__main__':
    t1 = datetime.datetime.now()
    
    # ----工作目录
    path = r'C:\Users\HP\OneDrive\temp\待处理的DSD\53614'
    ls_path = os.listdir(path)
    file = ls_path[-1]
    fpath = os.path.join(path, file)
    df = pd.read_excel(fpath, index_col = 0)
    # df1 = df.values
    nds = df.loc[:, 0:31].values
    prv = prv_cal(nds)
    
    # # ----多文件处理
    # import multiprocessing as mp
    # pool = mp.Pool(5)
    # results = [pool.apply_async(main, args=(os.path.join(path_path, a_path, 'ND_and_RR'),))
    #             for a_path in ls_path]
    # results1 = [p.get() for p in results]
    
    # ----多条目处理
    # import multiprocessing as mp
    # pool = mp.Pool(4)
    # results = [pool.apply_async(main, args=(os.path.join(path_path, a_path, 'ND_and_RR'),))
    #             for a_path in ls_path]
    # results1 = [p.get() for p in results]
 
    
    
    # nds=np.genfromtxt("nds.csv",delimiter=",")[:10]
    # df=prv_cal(nds)
    
    # result = main('G:\\BJ_DSD\\20170725-20170727\\ND_and_RR\\')
    
    t2 = datetime.datetime.now()
    print(t2-t1)
    # np.save("prvs.npy",prvs)
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
import pytmatrix
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive


path="20170822SY/dsd/"
codes=["54399","54416","54419","54431","54433","54499"]
for site_code in codes:
    # ----------[1] N(D)和R计算----------
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
    
    def readfromtxt(fid):
        # 统计文件中所有数据的个数
        f=open(fid)
        total=len(np.fromfile(f,count=-1,sep=" "))
        f.close()
        
        f=open(fid)
        # 第一行
        # 站号，纬度，经度，型号
        info=np.fromfile(f,count=4,sep=" ")
        p=3
        data={}
        while p < total:
            # 第二行
            # 时间（UTC）
            minute=str(np.fromfile(f,count=1,sep=" "))[1:-2]
            # 粒子种数 num_types
            num_types=int(np.fromfile(f,count=1,sep=" "))
            p+=2
            data[minute]={}
            data[minute]["num_types"]=num_types
            
            # 第三行
            # 粒子种类 types_d
            types=np.fromfile(f,count=num_types,sep=" ")
            p+=num_types
            types_d=(types%32-1).astype(np.int)
            types_d[types_d<0]=31
            data[minute]["types_d"]=types_d
            
            # 第四行
            # 各个种类对应的数量 num_per_type
            num_per_type=np.zeros(shape=32)
            temp=np.fromfile(f,count=num_types,sep=" ").astype(np.int)
            for i in range(num_types):  
                num_per_type[types_d[i]]+=temp[i]
            p+=num_types
            data[minute]["num_per_type"]=num_per_type
            
            
            # 根据类型，换算对应的：
            # 粒径 diameter
            dia=np.zeros(shape=32)
            dia[types_d]=diameter[types_d]
            data[minute]["diameter"]=dia
            # 速度 velocity
            vlc=np.zeros(shape=32)
            vlc[types_d]=velocity[types_d]
            data[minute]["velocity"]=vlc
            # 范围 spread
            sprd=np.zeros(shape=32)
            sprd[types_d]=spread[types_d]
            data[minute]["spread"]=sprd
            
            # 计算N(D) ND
            nd=num_per_type/0.0054/60/vlc/sprd
            nd[np.isnan(nd)]=0
            data[minute]["ND"]=nd
            # 计算R rain_rate
            r=6*np.pi/10**4/1*sum(vlc*dia**3*sprd*nd)
            data[minute]["rain_rate"]=r
            
            if p == total-1:
                break
        f.close()
        
        key_del=[]
        for key in data.keys():
            # 删除粒子总数Td小于10
            if np.sum(data[key]["num_per_type"])<10 and key not in key_del:
                key_del.append(key)
            # 删除RR小于0.1mm/h
            if data[key]["rain_rate"]<0.1 and key not in key_del:
                key_del.append(key)
            # 删除Dmax>8mm的数据，避免冰雹等的干扰
            if np.max(data[key]["diameter"])>8 and key not in key_del:
                key_del.append(key)
        for key in key_del:
            data.pop(key)
        
        return data, info
    
    nds={}
    rrs=[]
    for f in os.listdir(path):
        if site_code in f:
            print(f)
            # 读取并计算N(D)和R
            rawData,sitaInfo=readfromtxt(path+f)
            # 汇总每个时刻的N(D)和rain rate
            date=f[10:18]
            hour=int(f[18:20])
            for t in range(hour*100,hour*100+60):
                if str(t) in rawData.keys():# 有数据记录
                    nds[date+str(t)]=rawData[str(t)]["ND"]
                    rrs.append(rawData[str(t)]["rain_rate"])
                else:# 无数据记录
                    nds[date+str(t)]=np.zeros(shape=32)
                    rrs.append(0)
    
    # ----------[2] PRV计算----------
    bins=np.hstack((0,diameter+np.array(spread)/2))
    
    def thurai(D_eq):
      if D_eq < 0.7:
        return 1
      elif D_eq>=0.7 and D_eq <=1.5:
        return 1.173-0.5165*D_eq+0.4698*D_eq**(2)-0.1317*D_eq**3-8.5*10**(-3)*D_eq**4
      else:
        return 1.065-6.25*10**(-2)*D_eq-3.99*10**(-3)*D_eq**2+7.66*10**(-4)*D_eq**3-4.095*10**(-5)*D_eq**4
    
    def prv_cal(nds):
        z=[]
        kdp=[]
        ah=[]
        for key in nds.keys():
            nd=nds[key]
            # X波段，20℃作为雨滴温度
            scatterer=Scatterer(wavelength=tmatrix_aux.wl_X, m=refractive.m_w_20C[tmatrix_aux.wl_X])
            scatterer.psd_integrator=PSDIntegrator()
            # Thurai（2007）作为雨滴形状模型
            scatterer.psd_integrator.axis_ratio_func=lambda D: 1.0/thurai(D)
            # 删除Dmax>8mm的数据，避免冰雹等的干扰
            scatterer.psd_integrator.D_max=8
            scatterer.psd_integrator.geometries=(tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)
            #scatterer.or_pdf=orientation.gaussian_pdf(20.0)
            #scatterer.orient=orientation.orient_averaged_fixed
            scatterer.psd_integrator.init_scatter_table(scatterer)
            
            BinnedDSD=pytmatrix.psd.BinnedPSD(bins, nd)
            scatterer.psd=BinnedDSD
            scatterer.set_geometry(tmatrix_aux.geom_horiz_back)
            
            z.append(radar.refl(scatterer))
            scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)
            kdp.append(radar.Kdp(scatterer))
            ah.append(radar.Ai(scatterer))
        
        z=np.array(z)
        kdp=np.array(kdp)
        ah=np.array(ah)
        
        mask=np.where(z<1)
        z[mask]=1
        kdp[mask]=0
        ah[mask]=0
        z=10*np.log10(z)
        
        mask=np.where(z>55)
        z[mask]=0
        kdp[mask]=0
        ah[mask]=0
        
        prv=np.array([z,kdp,ah]).transpose()
        return prv
    
    lenth=len(nds)
    prv=prv_cal(nds)
    prv_3min=prv[np.arange(0,lenth,3)]
    
    rrs=np.array(rrs)
    loc=np.where(prv[:,0]==0)
    rrs[loc]=0
    rrs_3min=rrs[np.arange(0,lenth,3)]
    
    # 重构、拼接、储存
    rrs=rrs.reshape(lenth,1)
    data_per_min=np.hstack((prv,rrs))
    np.savetxt(path[:-4]+site_code+"_min.csv",data_per_min,delimiter=",")
    rrs_3min=rrs_3min.reshape(int(lenth/3),1)
    data_per_3min=np.hstack((prv_3min,rrs_3min))
    np.savetxt(path[:-4]+site_code+"_3min.csv",data_per_3min,delimiter=",")
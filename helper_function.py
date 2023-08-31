#!/usr/bin/env python
# coding: utf-8
# Writed in August 31, 2023
# author : Yukui Min

import ee
import collections
import os
import geemap
import subprocess
import math



def get_s2_sr_cld_col(aoi, start_date, end_date,MGRS_TILE,SENSING_ORBIT_NUMBER,CLOUDY_PIXEL_PERCENTAGE):
    # Import and filter S2 SR.
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')\
        .filterBounds(aoi)\
        .filterDate(start_date, end_date)\
        .filter(ee.Filter.eq('MGRS_TILE',MGRS_TILE))\
        .filter(ee.Filter.eq('SENSING_ORBIT_NUMBER',SENSING_ORBIT_NUMBER))\
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUDY_PIXEL_PERCENTAGE)))

    # Import and filter s2cloudless.
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')\
        .filterBounds(aoi)\
        .filterDate(start_date, end_date))\
        

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))

def add_cloud_bands(img):
    # Get s2cloudless image, subset the probability band.
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.lte(60).rename('cloud')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def cloud_mask_L8(img):
    qa = img.select('cloud')
    cloudShadowBitMask = 1<<3
    cloudsBitMask = 1<<5
    cloud = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))
    cloud = cloud.rename('cloud')
    return img.addBands(cloud,None,True)


def addtime(img):
    date = ee.Date(img.date())
    doy = ee.Number.parse(date.format('D'))
    img = img.set('DOY',doy)
    return img.addBands(ee.Image(doy).rename('doy').multiply(1000).toInt().divide(1000),None,True)


def addNDWI(img):
    ndwi = img.normalizedDifference(['green','nir']).rename('NDWI')
    return img.addBands(ndwi)

def addNDVI(img): 
    ndvi = img.normalizedDifference(['nir','red']).rename('NDVI')   
    return img.addBands(ndvi)

def applyScaleFactors(img):
    opticalBands = img.select('SR_B.').multiply(0.275).add(-2000).toInt()
    return img.addBands(opticalBands,None,True)
    return img

def toint(img):
    return img.multiply(1000).toInt().divide(1000)

def Spatial_register(img,refer_img):
    
    referRedBand = refer_img.select('red')
    imageRedBand = img.select('red')
    
    displacement = imageRedBand.displacement(
      referenceImage = referRedBand,
      maxOffset = 50.0,
      patchWidth = 100.0
    )

    registered = img.displace(displacement).copyProperties(img)
    registered = ee.Image(registered.copyProperties(img))
    
    return registered

def resample_L8(img,refer_img):

    crs = refer_img.select('red').projection().crs()
    resampled = img.resample('bicubic').reproject(crs = crs,scale = 10) 
    resampled = addtime(resampled)
    
    return resampled.copyProperties(img)

def band_adjust(img,SlopeL8,offsetL8):
    imgL8SR_bandadj = img.multiply(SlopeL8).add(offsetL8).multiply(1000).toInt().divide(1000)
    return imgL8SR_bandadj


#  Filling cloud-pixels by linear interpolation 
def linear_interpolate(img,img_before,img_after):
    
    #generate cloud mask
    mask = (img.select('cloud')).lte(0)
    unmask = (img.select('cloud')).gt(0)
  
    doy_b = img_before.select('doy')
    doy_a = img_after.select('doy')
    doy = img.select('doy')
    
    k = (img_after.subtract(img_before)).divide(doy_a.subtract(doy_b))
    x = doy.subtract(doy_b)
    ch = (x.multiply(k))
    
    linear = (img_before.add(ch))
    linear = linear.updateMask(mask)
    remain = (img.updateMask(unmask))
    
    # combine
    f = ee.ImageCollection([linear,remain])#image
    
    return f.mosaic().multiply(1000).toInt().divide(1000)



#name_No = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
# n = datesat_QY.size().getInfo()
# name_No = [k for k in range(0,n)]
# name_No = [str(k) for k in name_No]

def get_weight(datesat,i_obj,tide_flag,k):
    
    bands = ['doy','NDVI','NDWI']
    
    i = i_obj # i represent No of reconstructed Image
    n = k  # k represent sum of reconstruct Image
    # tide_flag  :if tide_flag == 1,mean neap tide region ,else mean spring tide region
    
    num = datesat.size().getInfo()
    No = [j for j in range(0,num)]
    name_No = [str(j) for j in No]

    datesat = datesat.select(bands)
    datesat_list = datesat.toList(datesat.size())
    
    datesat_array = datesat.toArray()
    datesat_doy = (datesat_array.arraySlice(1,0,1).arrayProject([0])).arrayFlatten([name_No])
    datesat_NDVI = (datesat_array.arraySlice(1,1,2).arrayProject([0])).arrayFlatten([name_No])
    T = (datesat_doy.select(name_No[i+n])).subtract(datesat_doy.select(name_No[i-n]))
    t0 = datesat_doy.select(name_No[i])
    
    #   create quality information for reconstruct Images by using NDWI index  #
    
    u = 1
    datesat_NDWI =  (datesat_array.arraySlice(1,2,3).arrayProject([0])).arrayFlatten([name_No])
    datesat_NDWI_weight = datesat_NDWI.where(datesat_NDWI.gt(0),0)
    datesat_NDWI_weight = datesat_NDWI.where(datesat_NDWI.lt(0),1)
    
    #####################################
                                                                  
    # weit_name = ['w1','w2','w3','w4','w5','w6']
    weit_name = [k for k in range(1,2*n+1)]
    weit_name = ['w'+str(k) for k in weit_name]
                                                                  
    weit_C = ee.Image.constant(ee.List.repeat(1,2*n)).float().rename(weit_name)
    
    if tide_flag == 1:
        for j in range(1,n+1):
            #  t = (t - t0)/(T - S)
            t_b = ((datesat_doy.select(name_No[i-j]).subtract(t0)).divide(T)).float()
            t_a = ((datesat_doy.select(name_No[i+j]).subtract(t0)).divide(T)).float()
            # NDWI 归一化变量（0 /1）
            ndwi_b = (datesat_NDWI_weight.select(name_No[i-j]))
            ndwi_a = (datesat_NDWI_weight.select(name_No[i+j]))
            
            wa = weit_C.select('w'+str(j+n)).multiply(t_a)
            wa = (((wa.divide(wa.abs())).subtract((wa.pow(2)).divide(0.1))).exp()).multiply(1000).toInt().divide(1000)
            wa = wa.multiply(ndwi_a)
            
            wb = weit_C.select('w'+str(n+1-j)).multiply(t_b)
            wb = (((wb.divide(wb.abs())).subtract((wb.pow(2)).divide(0.1))).exp()).multiply(1000).toInt().divide(1000)
            wb = wa.multiply(ndwi_b)
            
            wb = wb.rename('w'+str(n+1-j))
            wa = wa.rename('w'+str(j+n))
            weit_C = weit_C.addBands(wa,None,True)\
                           .addBands(wb,None,True)     
    else:
        for j in range(1,n+1):
            
            t_b = ((datesat_doy.select(name_No[i-j]).subtract(t0)).divide(T)).float()
            t_a = ((datesat_doy.select(name_No[i+j]).subtract(t0)).divide(T)).float()
            # NDWI 归一化变量（0 /1）
            ndwi_b = (datesat_NDWI_weight.select(name_No[i-j]))
            ndwi_a = (datesat_NDWI_weight.select(name_No[i+j]))
            
            wa = weit_C.select('w'+str(j+n)).multiply(t_a)
            wa = (((wa.pow(2)).divide(0.1).multiply(-1)).exp()).multiply(1000).toInt().divide(1000)  
            wa = wa.multiply(ndwi_a)
            
            wb = weit_C.select('w'+str(n+1-j)).multiply(t_b)
            wb = (((wb.pow(2)).divide(0.1).multiply(-1)).exp()).multiply(1000).toInt().divide(1000)
            wb = wa.multiply(ndwi_b)
            
            wb = wb.rename('w'+str(n+1-j))
            wa = wa.rename('w'+str(j+n))
            weit_C = weit_C.addBands(wa,None,True)\
                           .addBands(wb,None,True)
            
    def weit_Normalization(weit_img,n):
        
        sum_weit = (weit_img.toArray().toArray(1)).arrayAccum(0).arrayGet([2*n-1,0])
        Weight = (weit_img.divide(sum_weit)).multiply(1000).toInt().divide(1000) 
    
        return Weight
    
    Weit_img = weit_Normalization(weit_C,n)    
                      
    return Weit_img


def tide_gap_Filling(datesat,No,Bands,k):
    
    # datesat ---> input datesat
    # weight_img ---> weight image 
    # No ---> number of reconstructed Image in datesat                                                              
    # Bands ---> the bands you want gap filling
    # k ---> number of reconstruct Image for gap filling (default = 3)
    
    i = No
        
    No_B = [y for y in range(1,k+1)]
    No_B = [str(y) for y in No_B] #    No_B = ['1','2','3']
    No_A = [y for y in range(k+1,2*k+1)]
    No_A = [str(y) for y in No_A] #    No_A = ['4','5','6']


    NDVI_datesat = (datesat.select(Bands)).toArray().toArray(1)
    NDVI_Before = NDVI_datesat.arraySlice(0,i-k,i).arrayProject([0]).arrayFlatten([No_B])
    NDVI_After = NDVI_datesat.arraySlice(0,i+1,i+k+1).arrayProject([0]).arrayFlatten([No_A])
    NDVI_Collection = NDVI_Before.addBands(NDVI_After)
    
    datesat_list = datesat.toList(datesat.size())
    target_img = ee.Image(datesat_list.get(No))
    tide_region = target_img.updateMask(target_img.select('NDWI').gt(0))
   
    
    # Divide the image by NDWI into the following three areas
    
    # No tidal inundated area (NDWI =< 0)
    
    non_tide_region = target_img.updateMask(target_img.select('NDWI').lte(0))
    non_tide_region = non_tide_region.select(Bands).multiply(1000).toInt().divide(1000).rename(Bands+str(i))
    
    #tidal inundated ( 0 < NDWI < tide_thred)

    tide_thred = 0.4
    NDVI_Collection_Ntide = NDVI_Collection.updateMask(tide_region.select('NDWI').lt(tide_thred))
    weight_img_Ntide = (get_weight(datesat,No,0,k)).updateMask(tide_region.select('NDWI').lt(tide_thred))
    weighted_NDVI_Collection_Ntide = NDVI_Collection_Ntide.multiply(weight_img_Ntide).multiply(1000).toInt().divide(1000)
    composite_NDVI_Ntide = (weighted_NDVI_Collection_Ntide.toArray().toArray(1)).arrayAccum(0).arrayGet([2*k-1,0]).rename(Bands+str(i))  
    
    # tidal inundated (NDWI >= tide_thred)
   
    NDVI_Collection_Stide = NDVI_Collection.updateMask(tide_region.select('NDWI').gte(tide_thred))
    weight_img_Stide = (get_weight(datesat,No,1,k)).updateMask(tide_region.select('NDWI').gte(tide_thred))
    weighted_NDVI_Collection_Stide = NDVI_Collection_Stide.multiply(weight_img_Stide).multiply(1000).toInt().divide(1000)
    composite_NDVI_Stide = (weighted_NDVI_Collection_Stide.toArray().toArray(1)).arrayAccum(0).arrayGet([2*k-1,0]).rename(Bands+str(i))
    
    merge = ee.ImageCollection([composite_NDVI_Stide,composite_NDVI_Ntide,non_tide_region])
    merge = merge.mosaic().multiply(1000).toInt().divide(1000).rename('NDVI')

    return merge


def apply_TGF_to_timeseries(datesat,Bands):
    
    # apply Tide_Gap_filling to timeseries
    
    List = ee.List([])
    img_list = datesat.toList(datesat.size())
    n = datesat.size().getInfo()
    
    img0 = ee.Image(img_list.get(0)).select(['doy','NDVI'])
    img0 = img0.multiply(1000).toInt().divide(1000)
    img1 = ee.Image(img_list.get(1)).select(['doy']).addBands(tide_gap_Filling(datesat,1,Bands,1).rename('NDVI'))
    img2 = ee.Image(img_list.get(2)).select(['doy']).addBands(tide_gap_Filling(datesat,2,Bands,2).rename('NDVI'))
    List = ee.List([img0,img1,img2])
    
    for i in range(3,n-3):
        img =  ee.Image(img_list.get(i)).select(['doy']).addBands(tide_gap_Filling(datesat,i,Bands,3).rename('NDVI'))
        List = List.add(img)   
    
    img22 = ee.Image(img_list.get(n-3)).select(['doy']).addBands(tide_gap_Filling(datesat,n-3,Bands,2).rename('NDVI'))
    img11 = ee.Image(img_list.get(n-2)).select(['doy']).addBands(tide_gap_Filling(datesat,n-2,Bands,1).rename('NDVI'))
    img00 = ee.Image(img_list.get(n-1)).select(['doy','NDVI'])
    List = List.add(img22).add(img11).add(img00)
    
    collection = ee.ImageCollection.fromImages(List)
    
    return collection


def count_VI_dif(datesat):
    
    # input :datesat_dif
    
    n = datesat.size().getInfo()
    img_list = datesat.toList(n)
    list_dif = ee.List([])
    
    print(n)
    def VI_dif(img_before,img_after):
    
    #Calculate the difference between the NDVI and the corresponding DOY.

        img_before_ndvi = img_before.select('NDVI')
        img_after_ndvi = img_after.select('NDVI')
    
        NDVI_DIF = img_before_ndvi.subtract(img_after_ndvi)
        doy_mid = ((img_before.select('doy').add(img_after.select('doy'))).divide(2)).toInt()
        doy_after = img_after.select('doy')
        doy_before = img_before.select('doy')
        doy_dif = doy_after.subtract(doy_before).toInt()
        
        return NDVI_DIF.addBands(doy_mid.rename('H_doy'))\
                        .addBands(doy_dif.rename('dif_d'))\
                        .addBands(doy_after.rename('doy_after'))\
                        .addBands(doy_before.rename('doy_before'))

    for i in range(0,n-1):
        img1 = ee.Image(img_list.get(i))
        img2 = ee.Image(img_list.get(i+1))
        ndvi = VI_dif(img1,img2)
        list_dif = list_dif.add(ndvi)

    dif_datesat = ee.ImageCollection.fromImages(list_dif)
    
    return dif_datesat


# In[9]:


def NDVI_local_maximum(datesat):
    
    n = datesat.size().getInfo()
    datesat_list = datesat.toList(n)
    
    if n%3 == 0:
        num = int(n/3)
    else:
        num = int(n/3)+1

    mosaic_list = ee.List([])
    
    for i in range(0,num):
        
        if i < num-1:
            s = i*3
            e = s+3
        else:
            s = (num-1)*3 
            e = n      
        block_list = datesat_list.slice(s,e,1)
        collection = ee.ImageCollection.fromImages(block_list)
        NDVI_mosaic = collection.qualityMosaic('NDVI')
        mosaic_list = mosaic_list.add(NDVI_mosaic)
        
    composite_datesat = ee.ImageCollection.fromImages(mosaic_list)
    
    return composite_datesat


# In[10]:


def NDVI_lm_sort(datesat_QY_TGF_2021,dif_NDVI_composite_datesat):
    
    median_doy = (dif_NDVI_composite_datesat.qualityMosaic('NDVI')).select('H_doy')

    def count_dif_doy(image):
        doy = median_doy
        dif_doy = (image.select('doy')).subtract(doy).abs()
        return image.addBands(dif_doy.rename('dif_doy'))
    
    #  datesat_QY_TGF_2021 add dif_doy band  （ doy NDVI dif_doy）
    
    datesat_QY_TGF_2021_1 = datesat_QY_TGF_2021.map(count_dif_doy)

    # Sorting a two-dimensional array by dif_doy
    
    dif_doy = datesat_QY_TGF_2021_1.toArray().arraySlice(1,2,3)
    sorted_datesat_2021 =  datesat_QY_TGF_2021_1.toArray().arraySort(dif_doy)

    # Reconstruct the dataset by taking the 6 images with the smallest dif_doy and arrange them in the order of doy.
    
    filter_datesat = sorted_datesat_2021.arraySlice(0,0,6)
    doy = filter_datesat.arraySlice(1,0,1)
    filter_datesat_s = filter_datesat.arraySort(doy)
    
    return filter_datesat_s


def extract_potential_removal_period(filter_datesat_s):
    
    no = [k for k in range(0,6)]
    Name = ['NDVI'+str(k) for k in no]
    
    doy_collection = filter_datesat_s.arraySlice(1,0,1).arrayProject([0]).arrayFlatten([Name])
    NDVI_collection = filter_datesat_s.arraySlice(1,1,2).arrayProject([0]).arrayFlatten([Name])
    zlist = ee.List([])

    for j in range(0,6):
        NDVI = NDVI_collection.select(Name[j]).rename('NDVI')
        doy = doy_collection.select(Name[j]).rename('doy')
      #  t = t_collection.select(Name[j]).rename('t')
       # img = NDVI.addBands(t).addBands(doy)
        img = doy.addBands(NDVI)
        zlist = zlist.add(img)

    PCT_datesat = ee.ImageCollection.fromImages(zlist)
    
    return PCT_datesat



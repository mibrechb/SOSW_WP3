#funcs_ozone = require('users/michaelbrechbuehler/eawag:sos_water/tsm_retrieval/functions_ozone')
import ee
import functions_ozone as funcs_o3

def atmcorr_main_oli(img):
  """ Apply MAIN atmospheric correction (Page et al. 2018) to Landsat-8/9 OLI/OLI-2 sensor. """

  img = ee.Image(img)

  # oli bands
  bands = ['B1','B2','B3','B4','B5', 'B6', 'B7']
  bands_to_remove = ['B8', 'B9', 'B10', 'B11']
  bands_remaining = img.bandNames().removeAll(bands).removeAll(bands_to_remove)

  # tile geometry
  footprint = ee.Image(img).geometry()

  # dem
  #dem_srtm = ee.Image('USGS/SRTMGL1_003').clip(footprint)
  dem_alos = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2').select('DSM').mosaic().clip(footprint)
  DEM = dem_alos

  # ozone
  # load closest available day of gap-filled TOMS/MERGED dataset
  ozone = funcs_o3.ozone_interp
  DU_filt = ozone \
    .filterDate(img.date().advance(-14, 'day'), img.date().advance(14, 'day'))
  DU_nearest = DU_filt \
    .map(lambda du_img: du_img.set('timediff', du_img.date().difference(img.date(), 'day').abs())) \
    .sort('timediff', True) \
    .first()
  DU = DU_nearest
  #DU = ee.Image.constant(300)

  # julian Day
  imgDate = ee.Date(img.get('system:time_start'))
  FOY = ee.Date.fromYMD(imgDate.get('year'),1,1)
  JD = imgDate.difference(FOY,'day').int().add(1)

  # define pi
  pi = ee.Image(3.141592)

  # earth-Sun distance
  d = ee.Image.constant(img.get('EARTH_SUN_DISTANCE'))

  # sun elevation
  SunEl = ee.Image.constant(img.get('SUN_ELEVATION'))

  # sun azimuth
  SunAz = img.select('SAA').multiply(ee.Image(0.01))

  # satellite zenith
  SatZe = img.select('VZA').multiply(ee.Image(0.01))
  cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos()
  sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin()

  # satellite azimuth
  SatAz = img.select('VAA').multiply(ee.Image(0.01))

  # sun zenith
  SunZe = img.select('SZA').multiply(ee.Image(0.01))
  cosdSunZe = SunZe.multiply(pi.divide(ee.Image.constant(180))).cos(); # in degrees
  sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); # in degrees

  # relative azimuth
  RelAz = ee.Image(SunAz)
  cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos()

  # pressure calculation
  P = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)
  Po = ee.Image(1013.25)

  ## radiometric calibration 
  # radiance_mult_bands
  rad_mult = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_MULT_BAND_1')),
    ee.Image(img.get('RADIANCE_MULT_BAND_2')),
    ee.Image(img.get('RADIANCE_MULT_BAND_3')),
    ee.Image(img.get('RADIANCE_MULT_BAND_4')),
    ee.Image(img.get('RADIANCE_MULT_BAND_5')),
    ee.Image(img.get('RADIANCE_MULT_BAND_6')),
    ee.Image(img.get('RADIANCE_MULT_BAND_7'))]
    )).toArray(1)

  # radiance add band
  rad_add = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_ADD_BAND_1')),
    ee.Image(img.get('RADIANCE_ADD_BAND_2')),
    ee.Image(img.get('RADIANCE_ADD_BAND_3')),
    ee.Image(img.get('RADIANCE_ADD_BAND_4')),
    ee.Image(img.get('RADIANCE_ADD_BAND_5')),
    ee.Image(img.get('RADIANCE_ADD_BAND_6')),
    ee.Image(img.get('RADIANCE_ADD_BAND_7'))]
    )).toArray(1)

  # create an empty image to save new radiance bands to
  imgArr = img.select(bands).toArray().toArray(1)
  Ltoa = imgArr.multiply(rad_mult).add(rad_add)

  # esun
  # unclear why values are /10 compared to normal esun values
  ESUN = ee.Image.constant(197.24790954589844)\
  .addBands(ee.Image.constant(201.98426818847656))\
  .addBands(ee.Image.constant(186.12677001953125))\
  .addBands(ee.Image.constant(156.95257568359375))\
  .addBands(ee.Image.constant(96.04714965820312))\
  .addBands(ee.Image.constant(23.8833221450863))\
  .addBands(ee.Image.constant(8.04995873449635)).toArray().toArray(1)
  # # https://bleutner.github.io/RStoolbox/r/2016/01/26/estimating-landsat-8-esun-values
  # ESUN = ee.Image.constant(1895.33) \
  #   .addBands(ee.Image.constant(2004.57)) \
  #   .addBands(ee.Image.constant(1820.75)) \
  #   .addBands(ee.Image.constant(1549.49)) \
  #   .addBands(ee.Image.constant(951.76)) \
  #   .addBands(ee.Image.constant(247.55)) \
  #   .addBands(ee.Image.constant(85.46)).toArray().toArray(1)
  ESUN = ESUN.multiply(ee.Image(1))
  ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands])

  ## ozone correction 
  # ozone coefficients
  # Koontz, A., Flynn, C., Hodges, G., Michalsky, J., Barnard, J., 2013. Aerosol Optical Depth
  # Value-Added Product. U.S. Department of Energy, pp. 32.
  koz = ee.Image.constant(0.0039) \
    .addBands(ee.Image.constant(0.0218)) \
    .addBands(ee.Image.constant(0.1078)) \
    .addBands(ee.Image.constant(0.0608)) \
    .addBands(ee.Image.constant(0.0019)) \
    .addBands(ee.Image.constant(0)) \
    .addBands(ee.Image.constant(0)) \
    .toArray().toArray(1)

  # calculate ozone optical thickness
  Toz = koz.multiply(DU).divide(ee.Image.constant(1000))

  # calculate TOA radiance in the absense of ozone
  Lt = Ltoa.multiply(((Toz)).multiply((ee.Image.constant(1).divide(cosdSunZe)).add(ee.Image.constant(1).divide(cosdSatZe))).exp())

  # rayleigh optical thickness
  bandCenter = ee.Image(443).divide(1000) \
    .addBands(ee.Image(483).divide(1000)) \
    .addBands(ee.Image(561).divide(1000)) \
    .addBands(ee.Image(655).divide(1000)) \
    .addBands(ee.Image(865).divide(1000)) \
    .addBands(ee.Image(1609).divide(1000)) \
    .addBands(ee.Image(2201).divide(1000)) \
    .toArray().toArray(1)

  # create an empty image to save new Tr values to
  Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))))

  ## fresnel Reflection
  # specular reflection (s- and p- polarization states)
  theta_V = ee.Image(0.0000000001)
  sin_theta_j = sindSunZe.divide(ee.Image(1.333))

  theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi))

  theta_SZ = SunZe

  R_theta_SZ_s = (((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

  R_theta_V_s = ee.Image(0.0000000001)

  R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)))

  R_theta_V_p = ee.Image(0.0000000001)

  R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p))

  R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p))

  ## rayleigh scattering phase function
  # sun-sensor geometry
  theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract((sindSunZe).multiply(sindSatZe).multiply(cosdRelAz))

  theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi))

  theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz))

  theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi))

  cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos(); # in degrees

  cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos(); # in degrees

  Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))))

  Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))))

  # rayleigh scattering phase function
  Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos))

  # calulate Lr,
  denom = ee.Image(4).multiply(pi).multiply(cosdSatZe)
  Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom))

  # rayleigh corrected radiance
  Lrc = (Lt.divide(ee.Image(10))).subtract(Lr)
  LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands])

  # rayleigh corrected reflectance
  prc = Lrc.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe))
  prcImg = prc.arrayProject([0]).arrayFlatten([bands])
  #rhorc = prc.arrayProject([0]).arrayFlatten([bands])

  ## aerosol correction
  # bands in nm
  bands_nm = ee.Image(443) \
    .addBands(ee.Image(483)) \
    .addBands(ee.Image(561)) \
    .addBands(ee.Image(655)) \
    .addBands(ee.Image(865)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  # Lam in SWIR bands
  Lam_6 = LrcImg.select('B6')
  Lam_7 = LrcImg.select('B7')

  # calculate aerosol type
  eps = (((((Lam_7).divide(ESUNImg.select('B7'))).log()).subtract(((Lam_6).divide(ESUNImg.select('B6'))).log())).divide(ee.Image(2201).subtract(ee.Image(1609))))#.multiply(water_mask)

  # calculate multiple scattering of aerosols for each band
  Lam = (Lam_7).multiply(((ESUN).divide(ESUNImg.select('B7')))).multiply((eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2201)))).exp())

  # diffuse transmittance
  trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp()

  # compute water-leaving radiance
  Lw = Lrc.subtract(Lam).divide(trans)

  # water-leaving reflectance
  pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))
  pwImg = pw.arrayProject([0]).arrayFlatten([bands])

  # Rrs
  Rrs = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,5))

  # set negatives to None
  Rrs = Rrs.updateMask(Rrs.gt(0))

  # add auxilliary bands back in
  Rrs = Rrs.addBands(img.select(bands_remaining))

  # rename output bands
  # bands_to_rename = ['B1', 'B2', 'B3', 'B4', 'B5']
  # Rrs = ee.Image([
  #   ee.Image(Rrs.select(bands_to_rename).rename(['Rrs_'+band for band in bands_to_rename])), 
  #   ee.Image(Rrs.select(Rrs.bandNames().removeAll(bands_to_rename)))
  #   ])

  # copy image properties to output
  Rrs = Rrs.copyProperties(img) \
    .set('system:time_start', img.get('system:time_start'))

  return(Rrs)

def atmcorr_main_etm(img):
  """ Apply MAIN atmospheric correction (Page et al. 2018) to Landsat-7 ETM+ sensor. """
  
  img = ee.Image(img)

  # oli bands
  bands = ['B1','B2','B3','B4','B5','B7']
  bands_to_remove = ['B6_VCID_1', 'B6_VCID_2', 'B8']
  bands_remaining = img.bandNames().removeAll(bands).removeAll(bands_to_remove)
  
  # tile geometry
  footprint = img.geometry()

  # dem
  #dem_srtm = ee.Image('USGS/SRTMGL1_003').clip(footprint)
  dem_alos = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2').select('DSM').mosaic().clip(footprint)
  DEM = dem_alos

  # ozone
  # load closest available day of gap-filled TOMS/MERGED dataset
  ozone = funcs_o3.ozone_interp
  DU_filt = ozone \
    .filterDate(img.date().advance(-14, 'day'), img.date().advance(14, 'day'))
  DU_nearest = DU_filt \
    .map(lambda du_img: du_img.set('timediff', du_img.date().difference(img.date(), 'day').abs())) \
    .sort('timediff', True) \
    .first()
  DU = DU_nearest
  #DU = ee.Image.constant(300)

  # julian Day
  imgDate = ee.Date(img.get('system:time_start'))
  FOY = ee.Date.fromYMD(imgDate.get('year'),1,1)
  JD = imgDate.difference(FOY,'day').int().add(1)

  # define pi
  pi = ee.Image(3.141592)

  # earth-Sun distance
  d = ee.Image.constant(img.get('EARTH_SUN_DISTANCE'))

  # sun elevation
  SunEl = ee.Image.constant(img.get('SUN_ELEVATION'))

  # sun azimuth
  SunAz = img.select('SAA').multiply(ee.Image(0.01))

  # satellite zenith
  SatZe = img.select('VZA').multiply(ee.Image(0.01))
  cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos()
  sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin()

  # satellite azimuth
  SatAz = img.select('VAA').multiply(ee.Image(0.01))

  # sun zenith
  SunZe = img.select('SZA').multiply(ee.Image(0.01))
  cosdSunZe = SunZe.multiply(pi.divide(ee.Image.constant(180))).cos(); # in degrees
  sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); # in degrees

  # relative azimuth
  RelAz = ee.Image(SunAz)
  cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos()

  # pressure calculation
  P = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)
  Po = ee.Image(1013.25)

  ## radiometric calibration 
  # radiance_mult_bands
  rad_mult = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_MULT_BAND_1')),
    ee.Image(img.get('RADIANCE_MULT_BAND_2')),
    ee.Image(img.get('RADIANCE_MULT_BAND_3')),
    ee.Image(img.get('RADIANCE_MULT_BAND_4')),
    ee.Image(img.get('RADIANCE_MULT_BAND_5')),
    ee.Image(img.get('RADIANCE_MULT_BAND_7'))]
    )).toArray(1)

  # radiance add band
  rad_add = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_ADD_BAND_1')),
    ee.Image(img.get('RADIANCE_ADD_BAND_2')),
    ee.Image(img.get('RADIANCE_ADD_BAND_3')),
    ee.Image(img.get('RADIANCE_ADD_BAND_4')),
    ee.Image(img.get('RADIANCE_ADD_BAND_5')),
    ee.Image(img.get('RADIANCE_ADD_BAND_7'))]
    )).toArray(1)

  # create an empty image to save new radiance bands to
  imgArr = img.select(bands).toArray().toArray(1)
  Ltoa = imgArr.multiply(rad_mult).add(rad_add)

  # esun
  # https://www.usgs.gov/landsat-missions/using-usgs-landsat-level-1-data-product
  ESUN = ee.Image.constant(1970) \
    .addBands(ee.Image.constant(1842)) \
    .addBands(ee.Image.constant(1547)) \
    .addBands(ee.Image.constant(1044)) \
    .addBands(ee.Image.constant(225.7)) \
    .addBands(ee.Image.constant(82.06)).toArray().toArray(1).divide(ee.Image(10))
  ESUN = ESUN.multiply(ee.Image(1))
  ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands])

  ## ozone correction 
  # ozone coefficients
  # Koontz, A., Flynn, C., Hodges, G., Michalsky, J., Barnard, J., 2013. Aerosol Optical Depth
  # Value-Added Product. U.S. Department of Energy, pp. 32.
  koz = ee.Image.constant(0.0218) \
    .addBands(ee.Image.constant(0.1052)) \
    .addBands(ee.Image.constant(0.0536)) \
    .addBands(ee.Image.constant(0.0021)) \
    .addBands(ee.Image.constant(0)) \
    .addBands(ee.Image.constant(0)) \
    .toArray().toArray(1)

  # calculate ozone optical thickness
  Toz = koz.multiply(DU).divide(ee.Image.constant(1000))

  # calculate TOA radiance in the absense of ozone
  Lt = Ltoa.multiply(((Toz)).multiply((ee.Image.constant(1).divide(cosdSunZe)).add(ee.Image.constant(1).divide(cosdSatZe))).exp())

  # rayleigh optical thickness
  bandCenter = ee.Image(483).divide(1000) \
    .addBands(ee.Image(560).divide(1000)) \
    .addBands(ee.Image(662).divide(1000)) \
    .addBands(ee.Image(835).divide(1000)) \
    .addBands(ee.Image(1648).divide(1000)) \
    .addBands(ee.Image(2206).divide(1000)) \
    .toArray().toArray(1)

  # create an empty image to save new Tr values to
  Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))))

  ## fresnel Reflection
  # specular reflection (s- and p- polarization states)
  theta_V = ee.Image(0.0000000001)
  sin_theta_j = sindSunZe.divide(ee.Image(1.333))

  theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi))

  theta_SZ = SunZe

  R_theta_SZ_s = (((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

  R_theta_V_s = ee.Image(0.0000000001)

  R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)))

  R_theta_V_p = ee.Image(0.0000000001)

  R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p))

  R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p))

  ## rayleigh scattering phase function
  # sun-sensor geometry
  theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract((sindSunZe).multiply(sindSatZe).multiply(cosdRelAz))

  theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi))

  theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz))

  theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi))

  cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos(); # in degrees

  cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos(); # in degrees

  Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))))

  Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))))

  # rayleigh scattering phase function
  Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos))

  # calulate Lr,
  denom = ee.Image(4).multiply(pi).multiply(cosdSatZe)
  Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom))

  # rayleigh corrected radiance
  Lrc = (Lt.divide(ee.Image(10))).subtract(Lr)
  LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands])

  # rayleigh corrected reflectance
  prc = Lrc.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe))
  prcImg = prc.arrayProject([0]).arrayFlatten([bands])
  #rhorc = prc.arrayProject([0]).arrayFlatten([bands])

  ## aerosol correction
  # bands in nm
  bands_nm = ee.Image(483) \
    .addBands(ee.Image(560)) \
    .addBands(ee.Image(662)) \
    .addBands(ee.Image(835)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  # Lam in SWIR bands
  Lam_5 = LrcImg.select('B5')
  Lam_7 = LrcImg.select('B7')

  # calculate aerosol type
  eps = (((((Lam_7).divide(ESUNImg.select('B7'))).log()).subtract(((Lam_5).divide(ESUNImg.select('B5'))).log())).divide(ee.Image(2206).subtract(ee.Image(1648))))#.multiply(water_mask)

  # calculate multiple scattering of aerosols for each band
  Lam = (Lam_7).multiply(((ESUN).divide(ESUNImg.select('B7')))).multiply((eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2206)))).exp())

  # diffuse transmittance
  trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp()

  # compute water-leaving radiance
  Lw = Lrc.subtract(Lam).divide(trans)

  # water-leaving reflectance
  pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))
  pwImg = pw.arrayProject([0]).arrayFlatten([bands])

  # Rrs
  Rrs = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,4))

  # set negatives to None
  Rrs = Rrs.updateMask(Rrs.gt(0))

  # add auxilliary bands back in
  Rrs = Rrs.addBands(img.select(bands_remaining))

  # rename output bands
  # bands_to_rename = ['B1', 'B2', 'B3', 'B4']
  # Rrs = ee.Image([
  #   ee.Image(Rrs.select(bands_to_rename).rename(['Rrs_'+band for band in bands_to_rename])), 
  #   ee.Image(Rrs.select(Rrs.bandNames().removeAll(bands_to_rename)))
  #   ])

  # copy image properties to output
  Rrs = Rrs.copyProperties(img) \
    .set('system:time_start', img.get('system:time_start'))

  return(Rrs)

def atmcorr_main_msi(img):
  """ Apply MAIN atmospheric correction (Page et al. 2018) to Sentinel-2A/B MSI sensor. """

  img = ee.Image(img)

  # get platform
  platform = ee.String(img.get('SPACECRAFT_NAME'))

  # msi bands
  bands = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A', 'B11', 'B12']
  bands_to_remove = ['B9', 'B10']
  bands_remaining = img.bandNames().removeAll(bands).removeAll(bands_to_remove)

  # rescale
  scaling_factor = 0.0001
  rescale = img.select(bands).multiply(scaling_factor)

  # tile footprint
  footprint = rescale.geometry()

  # dem
  #dem_srtm = ee.Image('USGS/SRTMGL1_003').clip(footprint)
  dem_alos = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2').select('DSM').mosaic().clip(footprint)
  DEM = dem_alos

  # ozone
  # load closest available day of gap-filled TOMS/MERGED dataset
  ozone = funcs_o3.ozone_interp
  DU_filt = ozone \
    .filterDate(img.date().advance(-14, 'day'), img.date().advance(14, 'day'))
  DU_nearest = DU_filt \
    .map(lambda du_img: du_img.set('timediff', du_img.date().difference(img.date(), 'day').abs())) \
    .sort('timediff', True) \
    .first()
  DU = DU_nearest
  #DU = ee.Image.constant(300)

  # julian Day
  imgDate = ee.Date(img.get('system:time_start'))
  FOY = ee.Date.fromYMD(imgDate.get('year'),1,1)
  JD = imgDate.difference(FOY,'day').int().add(1)

  # define pi
  pi = ee.Image(3.141592)

  # earth-sun distance
  myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
  cosd = myCos.multiply(pi.divide(ee.Image(180))).cos()
  d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

  # sun azimuth
  SunAz = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint)

  # sun zenith
  SunZe = ee.Image.constant(img.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint)
  cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos(); # in degrees
  sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); # in degrees

  # sat zenith
  SatZe = ee.Image.constant(img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint)
  cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos()
  sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin()

  # sat azimuth
  SatAz = ee.Image.constant(img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint)

  # relative azimuth
  RelAz = SatAz.subtract(SunAz)
  cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos()

  # pressure
  P = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)
  Po = ee.Image(1013.25)

  # esun
  ESUN = ee.Image(ee.Array([ee.Image(img.get('SOLAR_IRRADIANCE_B1')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B2')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B3')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B4')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B5')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B6')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B7')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B8')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B8A')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B11')),
                    ee.Image(img.get('SOLAR_IRRADIANCE_B12'))]
                    )).toArray().toArray(1)

  ESUN = ESUN.multiply(ee.Image(1))
  ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands])

  # create empty array for the images
  imgArr = rescale.select(bands).toArray().toArray(1)

  # pTOA to Ltoa
  Ltoa = imgArr.multiply(ESUN).multiply(cosdSunZe).divide(pi.multiply(d.pow(2)))

  # band centers
  bandCenter_s2a = ee.Image(444).divide(1000) \
    .addBands(ee.Image(496).divide(1000)) \
    .addBands(ee.Image(560).divide(1000)) \
    .addBands(ee.Image(664).divide(1000)) \
    .addBands(ee.Image(704).divide(1000)) \
    .addBands(ee.Image(740).divide(1000)) \
    .addBands(ee.Image(782).divide(1000)) \
    .addBands(ee.Image(835).divide(1000)) \
    .addBands(ee.Image(865).divide(1000)) \
    .addBands(ee.Image(1613).divide(1000)) \
    .addBands(ee.Image(2202).divide(1000)) \
    .toArray().toArray(1)

  bandCenter_s2b = ee.Image(442).divide(1000) \
    .addBands(ee.Image(492).divide(1000)) \
    .addBands(ee.Image(559).divide(1000)) \
    .addBands(ee.Image(665).divide(1000)) \
    .addBands(ee.Image(703).divide(1000)) \
    .addBands(ee.Image(739).divide(1000)) \
    .addBands(ee.Image(779).divide(1000)) \
    .addBands(ee.Image(833).divide(1000)) \
    .addBands(ee.Image(864).divide(1000)) \
    .addBands(ee.Image(1610).divide(1000)) \
    .addBands(ee.Image(2185).divide(1000)) \
    .toArray().toArray(1)

  bandCenter = ee.Image(
        ee.Algorithms.If(platform.equals('Sentinel-2A'), bandCenter_s2a, bandCenter_s2b)
    )

  # ozone coefficients
  # Koontz, A., Flynn, C., Hodges, G., Michalsky, J., Barnard, J., 2013. Aerosol Optical Depth
  # Value-Added Product. U.S. Department of Energy, pp. 32.
  koz_s2a = ee.Image(0.0040) \
    .addBands(ee.Image(0.0244)) \
    .addBands(ee.Image(0.1052)) \
    .addBands(ee.Image(0.0516)) \
    .addBands(ee.Image(0.0208)) \
    .addBands(ee.Image(0.0112)) \
    .addBands(ee.Image(0.0079)) \
    .addBands(ee.Image(0.0021)) \
    .addBands(ee.Image(0.0019)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  koz_s2b = ee.Image(0.0037) \
    .addBands(ee.Image(0.0223)) \
    .addBands(ee.Image(0.1027)) \
    .addBands(ee.Image(0.0505)) \
    .addBands(ee.Image(0.0212)) \
    .addBands(ee.Image(0.0112)) \
    .addBands(ee.Image(0.0085)) \
    .addBands(ee.Image(0.0022)) \
    .addBands(ee.Image(0.0021)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  koz = ee.Image(
        ee.Algorithms.If(platform.equals('Sentinel-2A'), koz_s2a, koz_s2b)
    )

  # calculate ozone optical thickness
  Toz = koz.multiply(DU).divide(ee.Image(1000))

  # calculate TOA radiance in the absense of ozone
  Lt = Ltoa.multiply(((Toz)).multiply((ee.Image(1).divide(cosdSunZe)).add(ee.Image(1).divide(cosdSatZe))).exp())

  # rayleigh optical thickness
  Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))))

  # specular reflection (s- and p- polarization states)
  theta_V = ee.Image(0.0000000001)
  sin_theta_j = sindSunZe.divide(ee.Image(1.333))

  theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi))

  theta_SZ = SunZe

  R_theta_SZ_s = (((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

  R_theta_V_s = ee.Image(0.0000000001)

  R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)))

  R_theta_V_p = ee.Image(0.0000000001)

  R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p))

  R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p))

  # sun-sensor geometry
  theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract((sindSunZe).multiply(sindSatZe).multiply(cosdRelAz))

  theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi))

  theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz))

  theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi))

  cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos(); # in degrees

  cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos(); # in degrees

  Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))))

  Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))))

  # rayleigh scattering phase function
  Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos))

  # rayleigh radiance contribution
  denom = ee.Image(4).multiply(pi).multiply(cosdSatZe)
  Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom))

  # rayleigh corrected radiance
  Lrc = Lt.subtract(Lr)
  LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands])
  prcImg = (Lrc.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))
  prcImg = prcImg.arrayProject([0]).arrayFlatten([bands])

  ## aerosol Correction
  # bands in nm
  bands_nm_s2a = ee.Image(444) \
    .addBands(ee.Image(496)) \
    .addBands(ee.Image(560)) \
    .addBands(ee.Image(664)) \
    .addBands(ee.Image(703)) \
    .addBands(ee.Image(740)) \
    .addBands(ee.Image(782)) \
    .addBands(ee.Image(835)) \
    .addBands(ee.Image(865)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  bands_nm_s2b = ee.Image(442) \
    .addBands(ee.Image(492)) \
    .addBands(ee.Image(559)) \
    .addBands(ee.Image(665)) \
    .addBands(ee.Image(703)) \
    .addBands(ee.Image(739)) \
    .addBands(ee.Image(779)) \
    .addBands(ee.Image(833)) \
    .addBands(ee.Image(864)) \
    .addBands(ee.Image(0)) \
    .addBands(ee.Image(0)) \
    .toArray().toArray(1)

  bands_nm = ee.Image(
        ee.Algorithms.If(platform.equals('Sentinel-2A'), bands_nm_s2a, bands_nm_s2b)
    )

  # Lam in SWIR bands
  Lam_10 = LrcImg.select('B11'); # = 0
  Lam_11 = LrcImg.select('B12'); # = 0

  # calculate aerosol type
  eps = ((((Lam_11).divide(ESUNImg.select('B12'))).log()).subtract(((Lam_10).divide(ESUNImg.select('B11'))).log())).divide(ee.Image(2190).subtract(ee.Image(1610)))

  # calculate multiple scattering of aerosols for each band
  Lam = (Lam_11).multiply(((ESUN).divide(ESUNImg.select('B12')))).multiply((eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2190)))).exp())

  # diffuse transmittance
  trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp()

  # cmpute water-leaving radiance
  Lw = Lrc.subtract(Lam).divide(trans)

  # water-leaving reflectance
  pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))

  # remote sensing reflectance
  Rrs = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,9))

  # set negatives to None
  Rrs = Rrs.updateMask(Rrs.select('B1').gt(0))

  # add auxilliary bands back in
  Rrs = Rrs.addBands(img.select(bands_remaining))

  # rename output bands
  # bands_to_rename = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A']
  # Rrs = ee.Image([
  #   ee.Image(Rrs.select(bands_to_rename).rename(['Rrs_'+band for band in bands_to_rename])), 
  #   ee.Image(Rrs.select(Rrs.bandNames().removeAll(bands_to_rename)))
  #   ])

  # copy image properties to output
  Rrs = Rrs.copyProperties(img) \
    .set('system:time_start', img.get('system:time_start'))

  return(Rrs)
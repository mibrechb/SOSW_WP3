import ee
import functions_masking as funcs_masking
import functions_atmcorr_main as funcs_ac

def load_sr_imcoll(sensor, **kwargs):
  if sensor == 'oli': sr_imcoll = prep_sr_oli(**kwargs)
  elif sensor == 'msi': sr_imcoll = prep_sr_msi(**kwargs)
  elif sensor == 'etm': sr_imcoll = prep_sr_etm(**kwargs)
  return(ee.ImageCollection(sr_imcoll))

def load_rrs_imcoll(sensor, **kwargs):
  if sensor == 'oli': rrs_imcoll = prep_rrs_oli(**kwargs)
  elif sensor == 'msi': rrs_imcoll = prep_rrs_msi(**kwargs)
  elif sensor == 'etm': rrs_imcoll = prep_rrs_etm(**kwargs)
  return(ee.ImageCollection(rrs_imcoll))

def prep_rrs_msi(
      start_date, end_date,
      bounds=None, 
      cld_filt_thresh=75,
      qa_band='cs_cdf', clear_thresh=0.75,
      watermask=False
      ):
  """ Loads and processes Sentinel-2A/B MSI L1C collection from Rtoa to Rrs"""
  
  def _apply_cloudmask_cloudscore_msi(img):
    """ Apply Cloud Score+ cloudmask to Sentinel-2 MSI. """
    img = ee.Image(img)
    cld_probability = ee.Image(img.get('s2_csp')).select(qa_band).rename(qa_band)
    is_cloud = cld_probability.gte(clear_thresh).Not().rename('is_cloud')
    cld_bands = ee.Image([cld_probability, is_cloud])
    img_masked = img.select('.*B.*').updateMask(is_cloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  filter_msi = ee.Filter.And(
    ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds),
    ee.Filter.neq('system:index', '20160313T033552_20160313T034514_T47QRC'),
    ee.Filter.notNull(['MEAN_INCIDENCE_ZENITH_ANGLE_B5'])
  )

  s2cloudless = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
  s2_csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')

  ic_msi_rtoa = ee.ImageCollection('COPERNICUS/S2_HARMONIZED') \
    .filter(filter_msi) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring')

  ic_msi_rtoa_masked = ee.ImageCollection(ee.Join.saveFirst('s2_csp').apply(**{
        'primary': ic_msi_rtoa,
        'secondary': s2_csp,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    })) \
    .map(_apply_cloudmask_cloudscore_msi)
    
  if watermask: 
    ic_msi_rtoa_masked = ic_msi_rtoa_masked.map(funcs_masking.apply_watermask)
  
  ic_msi_rrs = ic_msi_rtoa_masked.map(funcs_ac.atmcorr_main_msi)
  return(ic_msi_rrs)

def prep_rrs_oli(start_date, end_date,
                 bounds=None,
                 cld_filt_thresh=75,
                 watermask=False
                 ):
  """ Loads and processes Landsat-8/9 OLI/OLI-2 collection from Rad to Rrs. """

  def _apply_cloudmask_oli(img):
    """ Apply QA_PIXEL cloudmask to Landsat-8/9 OLI/OLI-2. """
    
    def _bitwiseExtract(value, fromBit, toBit=None):
      if toBit == None: toBit = fromBit
      maskSize = ee.Number(1).add(toBit).subtract(fromBit)
      mask = ee.Number(1).leftShift(maskSize).subtract(1)
      return value.rightShift(fromBit).bitwiseAnd(mask)

    img = ee.Image(img)
    qa = img.select('QA_PIXEL')
    dilatedCloud =  _bitwiseExtract(qa, 1).rename('dilated_cloud')
    cirrus =        _bitwiseExtract(qa, 2).rename('cirrus')
    cloud =         _bitwiseExtract(qa, 3).rename('cloud')
    cloudShadow =   _bitwiseExtract(qa, 4).rename('shadow')
    clear =         _bitwiseExtract(qa, 6).rename('clear')
    isCloud = ((cloud.Not()) \
                .And(cirrus.Not()) \
                #.And(dilatedCloud.Not()) \
                .And(cloudShadow.Not()) \
                    .Not()).rename('is_cloud')
    cld_bands = ee.Image([dilatedCloud, cirrus, cloud, cloudShadow, isCloud])
    img_masked = img.select('.*B.*').updateMask(isCloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  filter_oli = ee.Filter.And(
    ee.Filter.lt('CLOUD_COVER', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds)
  )

# Top-of_atmosphere radiance (Rad)
  ic_oli1_rad = ee.ImageCollection("LANDSAT/LC08/C02/T1") \
    .filter(filter_oli) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') 
  ic_oli2_rad = ee.ImageCollection("LANDSAT/LC09/C02/T1") \
    .filter(filter_oli) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') 

  ic_oli_rad = ee.ImageCollection(ic_oli1_rad.merge(ic_oli2_rad)).sort('system:time_start')
  ic_oli_rad_masked = ic_oli_rad \
    .map(_apply_cloudmask_oli)
  
  if watermask: 
    ic_oli_rad_masked = ic_oli_rad_masked.map(funcs_masking.apply_watermask)

  ic_oli_rrs = ic_oli_rad_masked.map(funcs_ac.atmcorr_main_oli)
  return(ic_oli_rrs)

def prep_rrs_etm(start_date, end_date,
                 bounds=None,
                 cld_filt_thresh=75,
                 watermask=False
                 ):
  """ Loads and processes Landsat-7 ETM+ collection from Rad to Rrs. """

  def _apply_cloudmask_etm(img):
    """ Apply QA_PIXEL cloudmask to Landsat-7 ETM+. """
    
    def _bitwiseExtract(value, fromBit, toBit=None):
      if toBit == None: toBit = fromBit
      maskSize = ee.Number(1).add(toBit).subtract(fromBit)
      mask = ee.Number(1).leftShift(maskSize).subtract(1)
      return value.rightShift(fromBit).bitwiseAnd(mask)

    img = ee.Image(img)
    qa = img.select('QA_PIXEL')
    dilatedCloud =  _bitwiseExtract(qa, 1).rename('dilated_cloud')
    cirrus =        _bitwiseExtract(qa, 2).rename('cirrus')
    cloud =         _bitwiseExtract(qa, 3).rename('cloud')
    cloudShadow =   _bitwiseExtract(qa, 4).rename('shadow')
    clear =         _bitwiseExtract(qa, 6).rename('clear')
    isCloud = ((cloud.Not()) \
                .And(cirrus.Not()) \
                #.And(dilatedCloud.Not()) \
                .And(cloudShadow.Not()) \
                    .Not()).rename('is_cloud')
    cld_bands = ee.Image([dilatedCloud, cirrus, cloud, cloudShadow, isCloud])
    img_masked = img.select('.*B.*').updateMask(isCloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  filter_etm = ee.Filter.And(
    ee.Filter.lt('CLOUD_COVER', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds)
  )

# Top-of_atmosphere radiance (Rad)
  ic_etm_rad = ee.ImageCollection("LANDSAT/LE07/C02/T1") \
    .filter(filter_etm) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') \
    .sort('system:time_start')
  
  ic_etm_rad_masked = ic_etm_rad \
    .map(_apply_cloudmask_etm)
  
  if watermask: 
    ic_etm_rad_masked = ic_etm_rad_masked.map(funcs_masking.apply_watermask)

  ic_etm_rrs = ic_etm_rad_masked.map(funcs_ac.atmcorr_main_etm)
  return(ic_etm_rrs)

def prep_sr_msi(
      start_date, end_date,
      bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90), 
      cld_filt_thresh=75,
      qa_band='cs_cdf', clear_thresh=0.75,
      watermask=None
      ):
  """ Loads and processes Sentinel-2A/B MSI L2A SR collection. """
  
  def _apply_msi_offset(img):
    """ Apply scaling offset to Sentinel-2 MSI L2A SR collection. """
    img = ee.Image(img)
    img_scaled = img.select(['B.*', 'WVP']).multiply(0.0001)
    return(img.addBands(**{'srcImg': img_scaled, 'overwrite': True}))

  def _apply_cloudmask_cloudscore_msi(img):
    """ Apply Cloud Score+ cloudmask to Sentinel-2 MSI. """
    img = ee.Image(img)
    cld_probability = ee.Image(img.get('s2_csp')).select(qa_band).rename(qa_band)
    is_cloud = cld_probability.gte(clear_thresh).Not().rename('is_cloud')
    cld_bands = ee.Image([cld_probability, is_cloud])
    img_masked = img.select('.*B.*').updateMask(is_cloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  filter_msi = ee.Filter.And(
    ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds)
  )

  s2cloudless = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
  s2_csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')

  ic_msi_sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
    .filter(filter_msi) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') \
    .map(_apply_msi_offset)

  ic_msi_sr_masked = ee.ImageCollection(ee.Join.saveFirst('s2_csp').apply(**{
        'primary': ic_msi_sr,
        'secondary': s2_csp,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    })) \
    
  if watermask=='swm': 
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_swm_watermask)
  elif watermask=='index':
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_index_watermask)
  elif watermask=='scl':
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_scl_watermask)
  
  return(ic_msi_sr_masked)

def prep_sr_oli(start_date, end_date,
                bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90),
                cld_filt_thresh=75,
                watermask=None
                ):
  """ Loads and processes Landsat-8/9 OLI/OLI-2 SR collection. """

  def _apply_cloudmask_oli(img):
    """ Apply QA_PIXEL cloudmask to Landsat-8/9 OLI/OLI-2. """
    
    def _bitwiseExtract(value, fromBit, toBit=None):
      if toBit == None: toBit = fromBit
      maskSize = ee.Number(1).add(toBit).subtract(fromBit)
      mask = ee.Number(1).leftShift(maskSize).subtract(1)
      return value.rightShift(fromBit).bitwiseAnd(mask)

    img = ee.Image(img)
    qa = img.select('QA_PIXEL')
    dilatedCloud =  _bitwiseExtract(qa, 1).rename('dilated_cloud')
    cirrus =        _bitwiseExtract(qa, 2).rename('cirrus')
    cloud =         _bitwiseExtract(qa, 3).rename('cloud')
    cloudShadow =   _bitwiseExtract(qa, 4).rename('shadow')
    clear =         _bitwiseExtract(qa, 6).rename('clear')
    isCloud = ((cloud.Not()) \
                .And(cirrus.Not()) \
                #.And(dilatedCloud.Not()) \
                .And(cloudShadow.Not()) \
                    .Not()).rename('is_cloud')
    cld_bands = ee.Image([dilatedCloud, cirrus, cloud, cloudShadow, isCloud])
    img_masked = img.select('.*B.*').updateMask(isCloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  def _apply_scale_factors(image):
    image = ee.Image(image)
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(optical_bands, None, True).addBands(
        thermal_bands, None, True
    )

  filter_oli = ee.Filter.And(
    ee.Filter.lt('CLOUD_COVER', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds)
  )

# Top-of_atmosphere radiance (Rad)
  ic_oli1_sr = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
    .filter(filter_oli) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') 
  ic_oli2_sr = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2") \
    .filter(filter_oli) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') 

  ic_oli_sr= ee.ImageCollection(ic_oli1_sr.merge(ic_oli2_sr)).sort('system:time_start')
  ic_oli_sr_masked = ic_oli_sr \
    .map(_apply_cloudmask_oli) \
    .map(_apply_scale_factors)
  
  if watermask=='qa': 
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_masking.apply_qa_watermask)
  elif watermask=='swm':
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_masking.apply_swm_watermask)
  
  return(ic_oli_sr_masked)

def prep_sr_etm(start_date, end_date,
                bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90),
                cld_filt_thresh=75,
                watermask=False
                ):
  """ Loads and processes Landsat-7 ETM+ SR collection. """

  def _apply_cloudmask_etm(img):
    """ Apply QA_PIXEL cloudmask to Landsat-7 ETM+. """
    
    def _bitwiseExtract(value, fromBit, toBit=None):
      if toBit == None: toBit = fromBit
      maskSize = ee.Number(1).add(toBit).subtract(fromBit)
      mask = ee.Number(1).leftShift(maskSize).subtract(1)
      return value.rightShift(fromBit).bitwiseAnd(mask)

    img = ee.Image(img)
    qa = img.select('QA_PIXEL')
    dilatedCloud =  _bitwiseExtract(qa, 1).rename('dilated_cloud')
    cirrus =        _bitwiseExtract(qa, 2).rename('cirrus')
    cloud =         _bitwiseExtract(qa, 3).rename('cloud')
    cloudShadow =   _bitwiseExtract(qa, 4).rename('shadow')
    clear =         _bitwiseExtract(qa, 6).rename('clear')
    isCloud = ((cloud.Not()) \
                .And(cirrus.Not()) \
                #.And(dilatedCloud.Not()) \
                .And(cloudShadow.Not()) \
                    .Not()).rename('is_cloud')
    cld_bands = ee.Image([dilatedCloud, cirrus, cloud, cloudShadow, isCloud])
    img_masked = img.select('.*B.*').updateMask(isCloud.Not())
    return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
  
  def _apply_scale_factors(image):
    image = ee.Image(image)
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(optical_bands, None, True).addBands(
        thermal_bands, None, True
    )

  filter_etm = ee.Filter.And(
    ee.Filter.lt('CLOUD_COVER', cld_filt_thresh),
    ee.Filter.date(start_date, end_date),
    ee.Filter.bounds(bounds)
  )

# Top-of_atmosphere radiance (Rad)
  ic_etm_sr = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2") \
    .filter(filter_etm) \
    .map(lambda im: im.set('datestring', im.date().format('YYYY-MM-dd'))) \
    .distinct('datestring') \
    .sort('system:time_start')
  
  ic_etm_sr_masked = ic_etm_sr \
    .map(_apply_cloudmask_etm) \
    .map(_apply_scale_factors)
  
  if watermask=='qa': 
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_masking.apply_qa_watermask)
  elif watermask=='swm':
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_masking.apply_swm_watermask)

  return(ic_etm_sr_masked)
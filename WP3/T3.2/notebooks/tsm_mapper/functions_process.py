import ee
import functions_masking as funcs_masking
import functions_atmcorr_main as funcs_ac

def load_rrs_imcoll(sensor, **kwargs):
  if sensor == 'oli': rrs_imcoll = prep_rrs_oli(**kwargs)
  elif sensor == 'msi': rrs_imcoll = prep_rrs_msi(**kwargs)
  elif sensor == 'etm': rrs_imcoll = prep_rrs_etm(**kwargs)
  return(ee.ImageCollection(rrs_imcoll))

def prep_rrs_msi(
      start_date, end_date, bounds, 
      cld_filt_thresh=75,
      qa_band='cs_cdf', clear_thresh=0.75,
      mask_water=False
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
    
  if mask_water: ic_msi_rtoa_masked = ic_msi_rtoa_masked.map(funcs_masking.apply_watermask)
  
  ic_msi_rrs = ic_msi_rtoa_masked.map(funcs_ac.atmcorr_main_msi)
  return(ic_msi_rrs)

def prep_rrs_oli(bounds,
                 start_date, end_date,
                 cld_filt_thresh=75,
                 qa_band='cs_cdf', clear_thresh=0.75,
                 mask_water=False
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
  
  if mask_water: ic_oli_rad_masked = ic_oli_rad_masked.map(funcs_masking.apply_watermask)

  ic_oli_rrs = ic_oli_rad_masked.map(funcs_ac.atmcorr_main_oli)
  return(ic_oli_rrs)

def prep_rrs_etm(bounds,
                 start_date, end_date,
                 cld_filt_thresh=75,
                 qa_band='cs_cdf', clear_thresh=0.75,
                 mask_water=False
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
  
  if mask_water: ic_etm_rad_masked = ic_etm_rad_masked.map(funcs_masking.apply_watermask)

  ic_etm_rrs = ic_etm_rad_masked.map(funcs_ac.atmcorr_main_etm)
  return(ic_etm_rrs)
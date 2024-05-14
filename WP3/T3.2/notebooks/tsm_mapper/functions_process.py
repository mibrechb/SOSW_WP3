import ee
import functions_masking as funcs_masking
import functions_atmcorr_main as funcs_ac
import functions_metadata as funcs_metadata
import functions_turbidity as funcs_turb

def load_sr_imcoll(sensor, add_indices=False, add_ratios=False, **kwargs):
  if sensor == 'oli': sr_imcoll = prep_sr_oli(**kwargs)
  elif sensor == 'msi': sr_imcoll = prep_sr_msi(**kwargs)
  elif sensor == 'etm': sr_imcoll = prep_sr_etm(**kwargs)
  sr_imcoll = sr_imcoll.map(lambda img: ee.Image(img).set('sensor', sensor))
  if add_indices: sr_imcoll = sr_imcoll.map(funcs_turb.calc_indices)
  if add_ratios: sr_imcoll = sr_imcoll.map(funcs_turb.calc_ratios)
  return(ee.ImageCollection(sr_imcoll))

def load_rrs_imcoll(sensor, add_indices=False, add_ratios=False, **kwargs):
  if sensor == 'oli': rrs_imcoll = prep_rrs_oli(**kwargs)
  elif sensor == 'msi': rrs_imcoll = prep_rrs_msi(**kwargs)
  elif sensor == 'etm': rrs_imcoll = prep_rrs_etm(**kwargs)
  rrs_imcoll = rrs_imcoll.map(lambda img: ee.Image(img).set('sensor', sensor))
  if add_indices: rrs_imcoll = rrs_imcoll.map(funcs_turb.calc_indices)
  if add_ratios: rrs_imcoll = rrs_imcoll.map(funcs_turb.calc_ratios)
  return(ee.ImageCollection(rrs_imcoll))

def prep_rrs_msi(
      start_date, end_date,
      bounds=None, 
      cld_filt_thresh=75, cld_buffer=0,
      qa_band='cs_cdf', clear_thresh=0.75,
      watermask=False, water_buffer=0,
      harmonize_bnames=True
      ):
  """ Loads and processes Sentinel-2A/B MSI L1C collection from Rtoa to Rrs"""
   
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
    .map(funcs_masking.apply_csp_cloudmask(qa_band, clear_thresh, buffer_dist=cld_buffer))
    
  if watermask=='swm':
    ic_msi_rtoa_masked = ic_msi_rtoa_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_msi_rtoa_masked = ic_msi_rtoa_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))

  ic_msi_rrs = ic_msi_rtoa_masked.map(funcs_ac.atmcorr_main_msi)

  if harmonize_bnames:
    ic_msi_rrs = ic_msi_rrs.map(funcs_metadata.apply_bname_harmonization('msi', is_rrs=True))
  
  ic_msi_rrs = ic_msi_rrs.map(lambda img: ee.Image(img).set('platform', 'Sentinel-2'))

  return(ic_msi_rrs)

def prep_rrs_oli(start_date, end_date,
                 bounds=None,
                 cld_filt_thresh=75, cld_buffer=0,
                 watermask=False, water_buffer=0,
                 harmonize_bnames=True
                 ):
  """ Loads and processes Landsat-8/9 OLI/OLI-2 collection from Rad to Rrs. """
  
  def _set_platform(img):
    platform = ee.String(ee.Algorithms.If(
      img.getString('SPACECRAFT_ID').equals('LANDSAT_8'),
      'Landsat-8', 
      'Landsat-9'
      ))
    img = ee.Image(img).set('platform', platform)
    return(img)
    
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
    .map(funcs_masking.apply_qa_cloudmask(buffer_dist=cld_buffer))
  
  if watermask=='qa': 
    ic_oli_rad_masked = ic_oli_rad_masked.map(funcs_masking.apply_qa_watermask(buffer_dist=water_buffer))
  elif watermask=='swm':
    ic_oli_rad_masked = ic_oli_rad_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_oli_rad_masked = ic_oli_rad_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))

  ic_oli_rrs = ic_oli_rad_masked.map(funcs_ac.atmcorr_main_oli)

  if harmonize_bnames:
    ic_oli_rrs = ic_oli_rrs.map(funcs_metadata.apply_bname_harmonization('oli', is_rrs=True))

  ic_oli_rrs = ic_oli_rrs.map(_set_platform)

  return(ic_oli_rrs)

def prep_rrs_etm(start_date, end_date,
                 bounds=None,
                 cld_filt_thresh=75, cld_buffer=0,
                 watermask=False, water_buffer=0,
                 harmonize_bnames=True
                 ):
  """ Loads and processes Landsat-7 ETM+ collection from Rad to Rrs. """
  
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
    .map(funcs_masking.apply_qa_cloudmask(buffer_dist=cld_buffer))

  if watermask=='qa': 
    ic_etm_rad_masked = ic_etm_rad_masked.map(funcs_masking.apply_qa_watermask(buffer_dist=water_buffer))
  elif watermask=='swm':
    ic_etm_rad_masked = ic_etm_rad_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_etm_rad_masked = ic_etm_rad_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))

  ic_etm_rrs = ic_etm_rad_masked.map(funcs_ac.atmcorr_main_etm)

  if harmonize_bnames:
    ic_etm_rrs = ic_etm_rrs.map(funcs_metadata.apply_bname_harmonization('etm', is_rrs=True))

  ic_etm_rrs = ic_etm_rrs.map(lambda img: ee.Image(img).set('platform', 'Landsat-7'))

  return(ic_etm_rrs)

def prep_sr_msi(
      start_date, end_date,
      bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90), 
      cld_filt_thresh=75, cld_buffer=0,
      qa_band='cs_cdf', clear_thresh=0.75,
      watermask=None, water_buffer=0,
      harmonize_bnames=True
      ):
  """ Loads and processes Sentinel-2A/B MSI L2A SR collection. """
    
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
    .map(funcs_metadata.apply_msi_scaling)

  ic_msi_sr_masked = ee.ImageCollection(ee.Join.saveFirst('s2_csp').apply(**{
        'primary': ic_msi_sr,
        'secondary': s2_csp,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    })) \
    .map(funcs_masking.apply_csp_cloudmask(qa_band, clear_thresh, buffer_dist=cld_buffer))
    
  if watermask=='swm': 
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))
  elif watermask=='scl':
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_masking.apply_scl_watermask(buffer_dist=water_buffer))
  
  if harmonize_bnames:
    ic_msi_sr_masked = ic_msi_sr_masked.map(funcs_metadata.apply_bname_harmonization('msi'))

  ic_msi_sr_masked = ic_msi_sr_masked.map(lambda img: ee.Image(img).set('platform', 'Sentinel-2'))

  return(ic_msi_sr_masked)

def prep_sr_oli(start_date, end_date,
                bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90),
                cld_filt_thresh=75, cld_buffer=0,
                watermask=None, water_buffer=0,
                harmonize_bnames=True,
                ):
  """ Loads and processes Landsat-8/9 OLI/OLI-2 SR collection. """
 
  def _set_platform(img):
    platform = ee.String(ee.Algorithms.If(
      img.getString('SPACECRAFT_ID').equals('LANDSAT_8'),
      'Landsat-8', 
      'Landsat-9'
      ))
    img = ee.Image(img).set('platform', platform)
    return(img)

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
    .map(funcs_masking.apply_qa_cloudmask(buffer_dist=cld_buffer)) \
    .map(funcs_metadata.apply_l2c2_scaling)
  
  if watermask=='qa': 
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_masking.apply_qa_watermask(buffer_dist=water_buffer))
  elif watermask=='swm':
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))

  if harmonize_bnames:
    ic_oli_sr_masked = ic_oli_sr_masked.map(funcs_metadata.apply_bname_harmonization('oli'))

  ic_oli_sr_masked = ic_oli_sr_masked.map(_set_platform)

  return(ic_oli_sr_masked)

def prep_sr_etm(start_date, end_date,
                bounds=ee.Algorithms.GeometryConstructors.BBox(-180, -90, 180, 90),
                cld_filt_thresh=75, cld_buffer=0,
                watermask=False, water_buffer=0,
                harmonize_bnames=True
                ):
  """ Loads and processes Landsat-7 ETM+ SR collection. """

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
    .map(funcs_masking.apply_qa_cloudmask(buffer_dist=cld_buffer)) \
    .map(funcs_metadata.apply_l2c2_scaling)
  
  if watermask=='qa': 
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_masking.apply_qa_watermask(buffer_dist=water_buffer))
  elif watermask=='swm':
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_masking.apply_swm_watermask(buffer_dist=water_buffer))
  elif watermask=='index':
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_masking.apply_index_watermask(buffer_dist=water_buffer))

  if harmonize_bnames:
    ic_etm_sr_masked = ic_etm_sr_masked.map(funcs_metadata.apply_bname_harmonization('etm'))

  ic_etm_sr_masked = ic_etm_sr_masked.map(lambda img: ee.Image(img).set('platform', 'Landsat-7'))

  return(ic_etm_sr_masked)
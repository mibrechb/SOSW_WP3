import ee
import functions_masking as funcs_masking
import constants as c

def load_st_imcoll(start, end, roi, **kwargs):
  """ Loads Landsat L2C2 Surface Temperature ImageCollection. """

  def _apply_scale_factors(im):
    """ Apply scaling-factors/offsets. """
    optical_bands = im.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_bands = im.select('ST_B.*').multiply(0.00341802).add(149.0)
    uncertainty = im.select('ST_QA').multiply(0.01)
    cloud_distance = im.select('ST_CDIST').multiply(0.01)
    return im \
      .addBands(optical_bands, None, True) \
      .addBands(thermal_bands, None, True) \
      .addBands(uncertainty, None, True) \
      .addBands(cloud_distance, None, True)

  def _harmonize_bandnames(im):
    """ Harmonize bandnames to colornames. """
    sensor = im.getString('SPACECRAFT_ID')
    band_lut = ee.Dictionary(ee.Algorithms.If(
      ee.List(['LANDSAT_8', 'LANDSAT_9']).contains(sensor),
      c.bands_oli, # if OLI
      c.bands_tm_etm # if TM/ETM+
      )
    )
    opt_bands_og = [band_lut.get('blue'), band_lut.get('green'), \
                    band_lut.get('red'), band_lut.get('nir'), band_lut.get('swir1')]
    opt_bands = ['SR_blue', 'SR_green', 'SR_red', 'SR_nir', 'SR_swir1']
    qa_bands = ['QA_PIXEL', 'ST_QA', 'ST_EMIS', 'ST_CDIST']
    time = im.get('system:time_start')
    optical = im.select(opt_bands_og, opt_bands)
    thermal = im.select([band_lut.get('st')], ['ST'])
    return im \
      .addBands(optical, None, True) \
      .addBands(thermal, None, True) \
      .select(opt_bands+['ST']+qa_bands) \
      .set('system:time_start', time)

  filters = ee.Filter.And(
    ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'),
    ee.Filter.date(start, end),
    ee.Filter.bounds(roi)
  )

  ic_l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2").filter(filters)
  ic_l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filter(filters)
  ic_l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filter(filters) 
  ic_l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filter(filters)
  ic_l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filter(filters)
  
  ic_lst = ic_l4.merge(ic_l5).merge(ic_l7).merge(ic_l8).merge(ic_l9) \
    .map(_apply_scale_factors) \
    .map(_harmonize_bandnames) \
    .map(funcs_masking.apply_qa_cloudmask) \
    .map(funcs_masking.apply_qa_watermask)
    # .cast(**{
    #   'bandTypes': {'ST': 'uint16', 'QA_PIXEL': 'uint16', 'ST_QA': 'uint16', 'ST_EMIS': 'uint16', 'ST_CDIST': 'uint16'},
    #   'bandOrder': ['ST', 'QA_PIXEL', 'ST_QA', 'ST_EMIS', 'ST_CDIST']}) \
  
  return(ee.ImageCollection(ic_lst))
import ee
import constants as c


def apply_l2c2_scaling(img):
  """ Apply Landsat L2C2 scaling factors and offsets. """
  img = ee.Image(img)
  opticalBands = img.select('SR_B.').multiply(0.0000275).add(-0.2)
  thermalBands = img.select('ST_B.*').multiply(0.00341802).add(149.0)
  return img.addBands(opticalBands, None, True).addBands(thermalBands, None, True)

def apply_msi_scaling(img):
  """ Apply scaling offset to Sentinel-2 MSI L2A SR collection. """
  img = ee.Image(img)
  img_scaled = img.select(['B.*', 'WVP']).multiply(0.0001)
  return(img.addBands(**{'srcImg': img_scaled, 'overwrite': True}))

def apply_bname_harmonization(sensor, is_rrs=False):
  """ Rename bands to harmonized names. """
  if sensor == 'msi':
    bn = c.bands_msi
  elif sensor == 'etm':
    bn = c.bands_etm
  elif sensor =='oli':
    bn = c.bands_oli

  def _wrap_all(img):
    img = ee.Image(img)
    img = img.select(
      [bn['blue'], bn['green'], bn['red'], bn['nir'], \
       bn['swir1'], bn['swir2'], 'is_water', 'is_cloud'],
      ['blue', 'green', 'red', 'nir', \
       'swir1', 'swir2', 'is_water', 'is_cloud']
    )
    return(img)
  
  def _wrap_rrs(img):
    img = ee.Image(img)
    img = img.select(
      [bn['blue'], bn['green'], bn['red'], bn['nir'], \
       'is_water', 'is_cloud'],
      ['blue', 'green', 'red', 'nir', \
       'is_water', 'is_cloud']
    )
    return(img)

  if is_rrs:
    return(_wrap_rrs)
  else:
    return(_wrap_all)
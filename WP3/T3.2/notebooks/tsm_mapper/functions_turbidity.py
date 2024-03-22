import ee

def rename_landsat_l2c2_bands(img):
  platform = img.getString('SPACECRAFT_ID')
  band_dict = ee.Dictionary({
    'LANDSAT_4': [
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
      ['B_blue', 'B_green', 'B_red', 'B_nir', 'B_swir1', 'B_swir2', 'B_st', 'QA_PIXEL']],
    'LANDSAT_5': [
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
      ['B_blue', 'B_green', 'B_red', 'B_nir', 'B_swir1', 'B_swir2', 'B_st', 'QA_PIXEL']],
    'LANDSAT_7': [
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
      ['B_blue', 'B_green', 'B_red', 'B_nir', 'B_swir1', 'B_swir2', 'B_st', 'QA_PIXEL']],
    'LANDSAT_8': [
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
      ['B_blue', 'B_green', 'B_red', 'B_nir', 'B_swir1', 'B_swir2', 'B_st', 'QA_PIXEL']],
    'LANDSAT_9': [
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
      ['B_blue', 'B_green', 'B_red', 'B_nir', 'B_swir1', 'B_swir2', 'B_st', 'QA_PIXEL']]
  })
  bands_l2c2 = ee.List(band_dict.get(platform)).get(0)
  bands_generic = ee.List(band_dict.get(platform)).get(1)
  return img.select(bands_l2c2, bands_generic)

def oli_calc_tsm_markert_generic(img):
  """ Calculate TSM/SPM from Landsat SR (L2C2) imagery renamed with generic bandnames
      using Markert algorithm (Markert et al. 2018). """
  timestamp = img.get('system:time_start')
  b_g = img.select('B_green')
  b_r = img.select('B_red')
  ratio = (b_r.divide(b_g)).log()
  log_tsm = ratio.expression('a*e**(b*X+c)',
              {'X': ratio,
               'a': 1.90353307,
               'b': 1.44788939,
               'c': 0.62996462,
               'e': 2.718281828459045,
              })
  tsm = log_tsm.exp() \
    .rename(['add_tsm_markert']) \
    .set('system:time_start', timestamp)
  return img.addBands(tsm)

def oli_calc_tsm_markert_l2c2(img):
  """ Calculate TSM/SPM from Landsat SR (L2C2) imagery
      using Markert algorithm (Markert et al. 2018). """
  timestamp = img.get('system:time_start')
  platform = img.getString('SPACECRAFT_ID')
  band_dict = ee.Dictionary({
    'LANDSAT_4': ['SR_B2', 'SR_B3'],
    'LANDSAT_5': ['SR_B2', 'SR_B3'],
    'LANDSAT_7': ['SR_B2', 'SR_B3'],
    'LANDSAT_8': ['SR_B3', 'SR_B4'],
    'LANDSAT_9': ['SR_B3', 'SR_B4']
  })
  b_g = img.select(ee.List(band_dict.get(platform)).getString(0))
  b_r = img.select(ee.List(band_dict.get(platform)).getString(1))
  ratio = b_r.divide(b_g).log()
  log_tsm = ratio.expression('a*e**(b*X+c)',
              {'X': ratio,
               'a': 1.90353307,
               'b': 1.44788939,
               'c': 0.62996462,
               'e': 2.718281828459045,
              })
  tsm = log_tsm.exp() \
    .rename(['add_tsm_markert']) \
    .set('system:time_start', timestamp)
  return img.addBands(tsm)

def calc_spm_nechad(img):
  """ Calculate SPM using generic Nechad et al. (2010) algorithm.
      (https://doi.org/10.1016/j.rse.2009.11.022) """

  def _merge_bands(image, previous):
    return ee.Image(previous).addBands(image)

  # define coefficients based on platform
  # A-coeff.: Nechad et al. (2010) - Table 2
  # C-coeff.: Nechad et al. (2010) - Table 1
  coeff_dict = ee.Dictionary(ee.Algorithms.If(
      img.getString('SPACECRAFT_NAME'), 
      {'A_dict': {'B4': 355.85, 'B5': 468.13, 'B6': 1664.72, 'B7': 1802.62, 'B8': 2379.77, 'B8A': 2971.93}, 
       'C_dict': {'B4': 0.1728, 'B5': 0.1872, 'B6': 0.1973, 'B7': 0.2050, 'B8': 0.2103, 'B8A': 0.2115}
       }, #MSI
      ee.Algorithms.If(
          img.getString('SPACECRAFT_ID').equals('LANDSAT_7'), 
          {'A_dict': {'B3':342.56, 'B4':2379.77}, 
           'C_dict': {'B3':0.1719, 'B4':0.2103}
           }, #ETM+
          {'A_dict': {'B4': 289.29, 'B5': 2971.93}, 
           'C_dict': {'B4': 0.1686, 'B5': 0.2115}
           }, #OLI
        )
      )
    )

  def _calc_spm(band):
    A_Coeff = ee.Dictionary(coeff_dict.get('A_dict')).getNumber(band)
    C_Coeff = ee.Dictionary(coeff_dict.get('C_dict')).getNumber(band)
    Rrs = img.select([band])
    tsm = Rrs.expression('(A*pi*Rrs)/(1-pi*Rrs/C)+B',
            {'Rrs': Rrs,
             'A': A_Coeff,
             'B': 0,
             'C': C_Coeff,
             'pi': 3.14159265359,
            })
    return(tsm.rename(ee.String('add_spm_').cat(ee.String(band)).cat('_nechad')))

  # add spm bands
  bands  = ee.List(ee.Dictionary(coeff_dict.get('A_dict')).keys())
  img_tsm = ee.Image(ee.List(bands.map(_calc_spm)).iterate(_merge_bands, ee.Image())).slice(1)
  return(img.addBands(img_tsm))

def calc_tur_nechad(img):
  """ Calculate Turbidity using generic Nechad et al. (2009) algorithm.
      (https://doi.org/10.1117/12.830700) """

  def _merge_bands(image, previous):
    return ee.Image(previous).addBands(image)

  # define coefficients based on platform
  # A-coeff.: Nechad et al. (2009) - Table 2
  # C-coeff.: Nechad et al. (2010) - Table 1 (https://doi.org/10.1016/j.rse.2009.11.022)
  coeff_dict = ee.Dictionary(ee.Algorithms.If(
      img.getString('SPACECRAFT_NAME'),
      {'A_dict': {'B4': 282.95, 'B5': 333.50, 'B6': 1149.39, 'B7': 1304.10, 'B8': 1693.52, 'B8A': 2109.35}, 
       'C_dict': {'B4': 0.1728, 'B5': 0.1872, 'B6': 0.1973, 'B7': 0.2050, 'B8': 0.2103, 'B8A': 0.2115}
       }, #MSI
      ee.Algorithms.If(
          img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
          {'A_dict': {'B3':273.32, 'B4':1693.52}, 
           'C_dict': {'B3':0.1719, 'B4':0.2103}
           }, #ETM+
          {'A_dict': {'B4': 235.32, 'B5': 2109.35}, 
           'C_dict': {'B4': 0.1686, 'B5': 0.2115}
           } #OLI
        )
      )
    )

  def _calc_tur(band):
    A_Coeff = ee.Dictionary(coeff_dict.get('A_dict')).getNumber(band)
    C_Coeff = ee.Dictionary(coeff_dict.get('C_dict')).getNumber(band)
    Rrs = img.select([band])
    tur = Rrs.expression('(A*pi*Rrs)/(1-pi*Rrs/C)+B',
            {'Rrs': Rrs,
             'A': A_Coeff,
             'B': 0,
             'C': C_Coeff,
             'pi': 3.14159265359,
            })
    return(tur.rename(ee.String('add_tur_').cat(ee.String(band)).cat('_nechad')))

  # add tur bands
  bands  = ee.List(ee.Dictionary(coeff_dict.get('A_dict')).keys())
  img_tsm = ee.Image(ee.List(bands.map(_calc_tur)).iterate(_merge_bands, ee.Image())).slice(1)
  return(img.addBands(img_tsm))

def calc_tur_dogliotti(img):
  """ Calculate Turbidity from Rrs using switching algorithm
      by Dogliotti et al. (2015) (https'://doi.Org/10.1016/j.rse.2014.09.020) 
      and generic TUR coefficients by Nechad et al. (2009) (https://doi.org/10.1117/12.830700)."""

  # define coefficients based on platform
  # A-coeff.: Nechad et al. (2009) - Table 2
  # C-coeff.: Nechad et al. (2010) - Table 1 (https://doi.org/10.1016/j.rse.2009.11.022)
  coeff_dict = ee.Dictionary(ee.Algorithms.If(
    img.getString('SPACECRAFT_NAME'), ee.Dictionary({
      'A_dict': {'B4': 282.95, 'B5': 333.50,}, 
      'C_dict': {'B4': 0.1728, 'B5': 0.1872}
      }), #MSI
    ee.Algorithms.If(
        img.getString('SPACECRAFT_ID').equals('LANDSAT_7'), ee.Dictionary({
          'A_dict': {'B3':273.32, 'B4':1693.52}, 
          'C_dict': {'B3':0.1719, 'B4':0.2103}
          }), #ETM+
        ee.Dictionary({
          'A_dict': {'B4': 235.32, 'B5': 2109.35}, 
          'C_dict': {'B4': 0.1686, 'B5': 0.2115}
          }) #OLI
      )
    )
  )

  band_red = ee.Dictionary(coeff_dict.get('A_dict')).keys().get(0)
  A_red = ee.Dictionary(coeff_dict.get('A_dict')).getNumber(band_red)
  C_red = ee.Dictionary(coeff_dict.get('C_dict')).getNumber(band_red)
  band_nir = ee.Dictionary(coeff_dict.get('A_dict')).keys().get(1)
  A_nir = ee.Dictionary(coeff_dict.get('A_dict')).getNumber(band_nir)
  C_nir = ee.Dictionary(coeff_dict.get('C_dict')).getNumber(band_nir)

  def _calc_tur(img):
    # calculating turbidity estimates based on red and nir band
    T_red = img.expression('(A_red*pi*Rrs_red)/(1-pi*Rrs_red/C_red)+B_red',
            {'Rrs_red': img.select([band_red]),
             'A_red': A_red,
             'B_red': 0,
             'C_red': C_red,
             'pi': 3.14159265359,
            })
    T_nir = img.expression('(A_nir*pi*Rrs_nir)/(1-pi*Rrs_nir/C_nir)+B_nir',
            {'Rrs_nir': img.select([band_nir]),
             'A_nir': A_nir,
             'B_nir': 0,
             'C_nir': C_nir,
             'pi': 3.14159265359,
            })
    # blending turbidity estimates
    l_low = 0.05
    l_high = 0.07
    tur_blended = img.expression('(1-w)*T_red+w*T_nir',
            {'w': img.select([band_red]).unitScale(l_low, l_high).where(img.select([band_red]).gt(l_high), 1).where(img.select([band_red]).lt(l_low), 0),
             'T_red': T_red,
             'T_nir': T_nir,
            })
    return(tur_blended)

  img_tur = _calc_tur(img)
  return(img.addBands(img_tur.rename(ee.String('add_tur_dogliotti'))))

def calc_indices(img):
  """ Calculate spectral indices used for sediment retrieval. """
  img = ee.Image(img)
  bands = ee.List(ee.Algorithms.If(
      img.getString('SPACECRAFT_NAME'), # if MSI
      ['B2', 'B3', 'B4', 'B5'], #MSI
      ee.Algorithms.If(
          img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
          ['.*B1', '.*B2', '.*B3', '.*B4'], #ETM+
          ['.*B2', '.*B3', '.*B4', '.*B5']) #OLI
  ))
  b_red = bands.getString(2)
  b_green = bands.getString(1)
  b_blue = bands.getString(0)
  b_nir = bands.getString(3)
  # Normalized Suspended Material Index for OLI and MSI.
  nsmi = img.expression('(red+green-blue)/(red+green+blue)', {
    'red': img.select(b_red),
    'green': img.select(b_green),
    'blue': img.select(b_blue)}).rename('add_nsmi')
  # Calculate Normalized Difference Suspended Sediment Index for OLI and MSI (blue/nir-index).
  ndssi = img.normalizedDifference([b_blue, b_nir]).rename('add_ndssi')
  ndti = img.normalizedDifference([b_red, b_green]).rename('add_ndti')
  return(img.addBands([nsmi, ndssi, ndti]))

def calc_band_ratios(img):
  # Calculate band ratios used for sediment retrieval. #
  img = ee.Image(img)
  bands = ee.List(ee.Algorithms.If(
      img.getString('SPACECRAFT_NAME'), # if MSI
      ['B2', 'B3', 'B4', 'B5'], #MSI
      ee.Algorithms.If(
          img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
          ['.*B1', '.*B2', '.*B3', '.*B4'], #ETM+
          ['.*B2', '.*B3', '.*B4', '.*B5']) #OLI
    )
  )
  b_red = bands.getString(2),
  b_green = bands.getString(1)
  b_blue = bands.getString(0)
  b_nir = bands.getString(3)
  bratio_br = img.select(b_blue).divide(img.select(b_red)).rename('add_bratio_br')
  bratio_bn = img.select(b_blue).divide(img.select(b_nir)).rename('add_bratio_bn')
  bratio_gr = img.select(b_green).divide(img.select(b_red)).rename('add_bratio_gr')
  bratio_rgb = ((ee.Image.constant(1).divide(img.select(b_red))) \
                    .subtract((ee.Image.constant(1).divide(img.select([b_green]))))).multiply(img.select(b_blue)).rename('add_bratio_rgb')
  return(img.addBands([bratio_br, bratio_bn, bratio_gr, bratio_rgb]))
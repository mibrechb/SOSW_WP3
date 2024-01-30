import ee

def apply_swm_watermask(img):
  """ Apply global threshold water mask based on Sentinel Water Mask (SWM) index 
      (Milczarek et al. 2017). """
  bands = ee.List(ee.Algorithms.If(
      img.getString('SPACECRAFT_NAME'), # if MSI
      ['B2', 'B3', 'B8', 'B11'], #MSI
      ee.Algorithms.If(
         img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
         ['B1', 'B2', 'B4', 'B5'], #ETM+
         ['B2', 'B3', 'B5', 'B6']) #OLI
  ))
  b_blue = bands.get(0)
  b_green = bands.get(1)
  b_nir = bands.get(2)
  b_swir1 = bands.get(3)
  img = ee.Image(img)
  swm = img.expression('(blue+green)/(nir+swir1)', {
    'blue': img.select([b_blue]),
    'green': img.select([b_green]),
    'nir': img.select([b_nir]),
    'swir1': img.select([b_swir1])}).rename('swm')
  watermask = swm.gt(1.5).rename('watermask')
  return(img.updateMask(watermask).addBands([swm, watermask]))

def apply_watermask(img):
    """ Apply water mask based on Fan et al. 2022 (10.1109/IGARSS46834.2022.9883460) method applied on the
        Sentinel Water Mask (SWM) index (Milczarek et al. 2017)."""

    def _otsu(img):
        """ Returns optimal threshold value based on histogram and Otsu (1979) method 
        for image segmentation. Expects cloud-free input with clear land/water difference (e.g. MNDWI). """
        histogram = ee.Dictionary(img.reduceRegion(**{
            'reducer': ee.Reducer.histogram(100),
            'geometry': img.geometry(), 
            'scale': 20,
            'bestEffort': True
            }).values().get(0))

        counts = ee.Array(histogram.get('histogram'))
        means = ee.Array(histogram.get('bucketMeans'))
        size = means.length().get([0])
        total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
        sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
        mean = sum.divide(total)
        indices = ee.List.sequence(1, size)
    
        def _bss(i):
            aCounts = counts.slice(0, 0, i)
            aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
            aMeans = means.slice(0, 0, i)
            aMean = aMeans.multiply(aCounts) \
                .reduce(ee.Reducer.sum(), [0]).get([0]) \
                .divide(aCount)
            bCount = total.subtract(aCount)
            bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount)
            return aCount.multiply(aMean.subtract(mean).pow(2)).add(
            bCount.multiply(bMean.subtract(mean).pow(2)))

        # return threshold value corresponding to the maximum BSS
        bss = indices.map(_bss) 
        return means.sort(bss).get([-1])

    bands = ee.List(ee.Algorithms.If(
        img.getString('SPACECRAFT_NAME'), # if MSI
        ['B2', 'B3', 'B8', 'B11'], #MSI
        ee.Algorithms.If(
            img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
            ['B1', 'B2', 'B4', 'B5'], #ETM+
            ['B2', 'B3', 'B5', 'B6']) #OLI
    ))
    b_blue = bands.get(0)
    b_green = bands.get(1)
    b_nir = bands.get(2)
    b_swir1 = bands.get(3)
    img = ee.Image(img)
    swm = img.expression('(blue+green)/(nir+swir1)', {
        'blue': img.select([b_blue]),
        'green': img.select([b_green]),
        'nir': img.select([b_nir]),
        'swir1': img.select([b_swir1])}).rename('swm')
        
    wi = swm
    swir = ee.List(ee.Algorithms.If(
        img.getString('SPACECRAFT_NAME'), # if MSI
        'B11', #MSI
        ee.Algorithms.If(
            img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
            'B5', #ETM+
            'B6') #OLI+
    ))
    r_swir = img.select(swir)
    t_swir = 0.1
    t_wir = 0.4
    t_s =   _otsu(wi)
    wi_water = wi.updateMask(wi.gt(t_s))
    wi_land =  wi.updateMask(wi.lt(t_s))
    # m1 =    wi_water.reduceRegion(**{'reducer': ee.Reducer.median(), 'bestEffort':True}).values().getNumber(0)
    # std1 =  wi_water.reduceRegion(**{'reducer': ee.Reducer.stdDev(), 'bestEffort':True}).values().getNumber(0)
    # m2 =    wi_land.reduceRegion(**{'reducer': ee.Reducer.median(), 'bestEffort':True}).values().getNumber(0)
    # std2 =  wi_land.reduceRegion(**{'reducer': ee.Reducer.stdDev(), 'bestEffort':True}).values().getNumber(0)
    # t_pure =    ee.List([(m1.add(std1)).multiply(0.5), t_s]).reduce(ee.Reducer.max())
    # t_mixed =   ee.List([(m2.add(std2)).multiply(0.5), t_s]).reduce(ee.Reducer.min())
    # c1 =    wi.gt(ee.Image.constant(t_pure)).rename('c1')
    # p1 =    (wi.gt(ee.Image.constant(t_mixed))).And(wi.lte(ee.Image.constant(t_pure))).rename('p1')
    # p2 =    p1.And(r_swir.lt(t_swir)).rename('p2')
    # wir =   wi.subtract(wi.reduceNeighborhood(**{'reducer': ee.Reducer.min(), 'kernel': ee.Kernel.square(**{'radius': 2})}))
    # c2 =    c1.Or(p2.And(wir.gt(ee.Image.constant(t_wir)))).rename('c2')
    mask = wi.gt(t_s).rename('watermask')
    return(img.updateMask(mask).addBands(mask))
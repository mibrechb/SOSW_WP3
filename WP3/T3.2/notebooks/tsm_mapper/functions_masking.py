import ee
import constants as c

def apply_growth(img, distance):
    """ Applies an efficient buffering method to a 0/1 image to grow 1. """
    img = ee.Image(img)
    nbh_sz = (ee.Number(distance).divide(30).add(3)).ceil()
    d = (img.fastDistanceTransform(nbh_sz).sqrt() \
         .multiply(ee.Image.pixelArea().sqrt()))
    return d.lte(distance).rename('buffer')

def apply_erosion(img, distance):
    """ Applies an efficient neg. buffering method to a binary image to erode 1. """
    img = ee.Image(img)
    return apply_growth(img.Not(), distance).Not()

def apply_csp_cloudmask(qa_band, clear_thresh, buffer_dist=0):
    """ Apply Cloud Score+ cloudmask to Sentinel-2 MSI. """
    def _wrap(img):
        img = ee.Image(img)
        cld_probability = ee.Image(img.get('s2_csp')).select(qa_band).rename(qa_band)
        is_cloud = cld_probability.gte(clear_thresh).Not().rename('is_cloud')
        is_cloud_buffered = apply_growth(is_cloud, buffer_dist).rename('is_cloud')
        is_cloud = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_cloud_buffered,
            is_cloud
        ))
        cld_bands = ee.Image([cld_probability, is_cloud])
        img_masked = img.select('.*B.*').updateMask(is_cloud.Not())
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
    return(_wrap)

def apply_qa_cloudmask(buffer_dist=0):
    """ Apply QA_PIXEL cloudmask to Landsat Collection 2 Level 2. """
    def _wrap(img):
        def _bitwiseExtract(value, fromBit, toBit=None):
            if toBit == None: toBit = fromBit
            maskSize = ee.Number(1).add(toBit).subtract(fromBit)
            mask = ee.Number(1).leftShift(maskSize).subtract(1)
            return value.rightShift(fromBit).bitwiseAnd(mask)

        img = ee.Image(img)
        qa = img.select('QA_PIXEL')
        dilated_cloud =  _bitwiseExtract(qa, 1).rename('dilated_cloud')
        cirrus =         _bitwiseExtract(qa, 2).rename('cirrus')
        cloud =          _bitwiseExtract(qa, 3).rename('cloud')
        cloud_shdw =     _bitwiseExtract(qa, 4).rename('cloud_shdw')
        clear =          _bitwiseExtract(qa, 6).rename('clear')
        is_cloud = cloud \
                    .Or(cirrus) \
                    .Or(dilated_cloud) \
                    .Or(cloud_shdw) \
                        .rename('is_cloud')
        
        is_cloud_buffered = apply_growth(is_cloud, buffer_dist).rename('is_cloud')
        is_cloud = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_cloud_buffered,
            is_cloud
        ))
        cld_bands = ee.Image([dilated_cloud, cirrus, cloud, cloud_shdw, is_cloud])
        img_masked = img.select('.*B.*').updateMask(is_cloud.Not())
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(cld_bands))
    return(_wrap)

def apply_qa_watermask(buffer_dist=0):
    """ Apply Landsat water mask based on QA flag. """
    def _wrap(img):
        def _bitwiseExtract(value, fromBit, toBit=None):
            if toBit == None: toBit = fromBit
            maskSize = ee.Number(1).add(toBit).subtract(fromBit)
            mask = ee.Number(1).leftShift(maskSize).subtract(1)
            return value.rightShift(fromBit).bitwiseAnd(mask)
        
        img = ee.Image(img)
        qa = img.select('QA_PIXEL')
        is_water =  _bitwiseExtract(qa, 7).eq(1).rename('is_water')
        is_water_buffered = apply_erosion(is_water, buffer_dist).rename('is_water')
        is_water = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_water_buffered,
            is_water
        ))
        img_masked = img.select('.*B.*').updateMask(is_water)
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(is_water))
    return(_wrap)

def apply_scl_watermask(buffer_dist=0):
    """ Apply Sentinel-2 water mask based on L2A Scene Classification (SCL) flag. """
    def _wrap(img):
        img = ee.Image(img)
        is_water = img.select('SCL').eq(6).rename('is_water')
        is_water_buffered = apply_erosion(is_water, buffer_dist).rename('is_water')
        is_water = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_water_buffered,
            is_water
        ))
        img_masked = img.select('.*B.*').updateMask(is_water)
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(is_water))
    return(_wrap)

def apply_index_watermask_dethier(buffer_dist=0):
    """ Apply index based water mask based on Dethier et al. (10.1126/science.abn7980) """
    def _wrap(img):
        bands = ee.Dictionary(ee.Algorithms.If(
            img.getString('SPACECRAFT_NAME'), # if MSI
            c.bands_msi, #MSI
            ee.Algorithms.If(
                img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
                c.bands_etm, #ETM+
                c.bands_oli #OLI
            )
        ))
        img = ee.Image(img)
        mndwi = img.expression('(green-swir1)/(green+swir1)', {
            'green': img.select([bands.get('green')]),
            'swir1': img.select([bands.get('swir1')])
        }).rename('mndwi')
        is_water = mndwi.gt(0) \
            .And(img.select([bands.get('swir1')]).lt(0.05)) \
            .And((img.select([bands.get('blue')]).add(img.select([bands.get('green')]))).lt(0.5)) \
            .rename('is_water')
        is_water_buffered = apply_erosion(is_water, buffer_dist).rename('is_water')
        is_water = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_water_buffered.rename('is_water'),
            is_water
        ))
        img_masked = img.select('.*B.*').updateMask(is_water)
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(is_water))
    return(_wrap)

def apply_index_watermask(buffer_dist=0):
    """ Apply index based water mask based on MNDWI threshold. """
    def _wrap(img):
        bands = ee.Dictionary(ee.Algorithms.If(
            img.getString('SPACECRAFT_NAME'), # if MSI
            c.bands_msi, #MSI
            ee.Algorithms.If(
                img.getString('SPACECRAFT_ID').equals('LANDSAT_7'),
                c.bands_etm, #ETM+
                c.bands_oli #OLI
            )
        ))
        img = ee.Image(img)
        mndwi = img.expression('(green-swir1)/(green+swir1)', {
            'green': img.select([bands.get('green')]),
            'swir1': img.select([bands.get('swir1')])
        }).rename('mndwi')
        is_water = mndwi.gt(0).rename('is_water')
        is_water_buffered = apply_erosion(is_water, buffer_dist).rename('is_water')
        is_water = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_water_buffered,
            is_water
        ))
        img_masked = img.select('.*B.*').updateMask(is_water)
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(is_water))
    return(_wrap)

def apply_swm_watermask(buffer_dist=0):
    """ Apply water mask based on Fan et al. 2022 (10.1109/IGARSS46834.2022.9883460) method applied on the
        Sentinel Water Mask (SWM) index (Milczarek et al. 2017)."""
    def _wrap(img):
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
        is_water = wi.gt(t_s).rename('is_water')
        is_water_buffered = apply_erosion(is_water, buffer_dist).rename('is_water')
        is_water = ee.Image(ee.Algorithms.If(
            ee.Number(buffer_dist).neq(0),
            is_water_buffered,
            is_water
        ))
        img_masked = img.select('.*B.*').updateMask(is_water)
        return(img.addBands(**{'srcImg': img_masked, 'overwrite': True}).addBands(is_water))
    return(_wrap)
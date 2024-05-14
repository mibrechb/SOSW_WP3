import ee
import constants as c

def apply_growth(img, distance):
    """ Applies an efficient buffering method to a 0/1 image to grow 1. """
    img = img.unmask()
    nbh_sz = (ee.Number(distance).divide(30).add(3)).ceil()
    d = (img.fastDistanceTransform(nbh_sz).sqrt() \
         .multiply(ee.Image.pixelArea().sqrt()))
    return d.lte(distance).rename('buffer')

def apply_erosion(img, distance):
    """ Applies an efficient neg. buffering method to a binary image to erode 1. """
    img = img.unmask()
    return apply_growth(img.Not(), distance).Not()

def apply_qa_cloudmask(buffer_dist=0):
    """ Apply QA_PIXEL cloudmask to Landsat Collection 2 Level 2 with optional buffer. """
    def _wrap(img):
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
                    .And(dilatedCloud.Not()) \
                    .And(cloudShadow.Not()) \
                        .Not()).rename('is_cloud')
        
        buffer = apply_growth(isCloud, buffer_dist).rename('buffer')
        isCloud = ee.Image(ee.Algorithms.If(
            buffer_dist,
            isCloud.Or(buffer).rename('is_cloud'),
            isCloud
        ))
        cld_bands = ee.Image([dilatedCloud, cirrus, cloud, cloudShadow, isCloud])
        opt_masked = img.select('SR_.*').updateMask(isCloud.Not().And(buffer.Not()))
        st_masked = img.select('SR_.*').updateMask(isCloud.Not().And(buffer.Not()))
        return img \
            .addBands(opt_masked, None, True) \
            .addBands(st_masked, None, True) \
            .addBands(cld_bands, None, True)
    return _wrap

def apply_qa_watermask(buffer_dist=0):
    """ Apply Landsat water mask based on QA flag with optional negative buffer. """
    def _wrap(img):
        def _bitwiseExtract(value, fromBit, toBit=None):
            if toBit == None: toBit = fromBit
            maskSize = ee.Number(1).add(toBit).subtract(fromBit)
            mask = ee.Number(1).leftShift(maskSize).subtract(1)
            return value.rightShift(fromBit).bitwiseAnd(mask)
        img = ee.Image(img)
        qa = img.select('QA_PIXEL')
        isWater =  _bitwiseExtract(qa, 7).eq(1).rename('is_water')
        buffer = apply_erosion(isWater, buffer_dist).rename('is_water')
        isWater = ee.Image(ee.Algorithms.If(
            buffer_dist,
            isWater.Or(buffer).rename('is_water'),
            isWater
        ))
        return img.updateMask(isWater) \
                    .addBands([isWater])
    return _wrap

def apply_index_watermask(buffer_dist=0):
    """ Apply index based water mask based on Dethier et al. (10.1126/science.abn7980)
        with optional negative buffer."""
    def _wrap(img):
        img = ee.Image(img)
        mndwi = img.expression('(green-swir1)/(green+swir1)', {
            'green': img.select('SR_green'),
            'swir1': img.select('SR_swir1')}).rename('mndwi')
        isWater = mndwi.gt(0) \
            .And(img.select('SR_swir1').lt(0.05)) \
            .And((img.select('SR_blue').add(img.select('SR_green'))).lt(0.5)) \
            .rename('is_water')
        buffer = apply_erosion(isWater, buffer_dist).rename('is_water')
        isWater = ee.Image(ee.Algorithms.If(
            buffer_dist,
            isWater.Or(buffer).rename('is_water'),
            isWater
        ))
        return img.updateMask(isWater) \
                    .addBands(isWater)
    return _wrap

def apply_swm_watermask(buffer_dist=0):
    """ Apply water mask based on Fan et al. 2022 (10.1109/IGARSS46834.2022.9883460) method applied on the
        Sentinel Water Mask (SWM) index (Milczarek et al. 2017)."""
    def _wrap(im):
        def _otsu(im):
            """ Returns optimal threshold value based on histogram and Otsu (1979) method 
            for image segmentation. Expects cloud-free input with clear land/water difference (e.g. MNDWI). """
            histogram = ee.Dictionary(im.reduceRegion(**{
                'reducer': ee.Reducer.histogram(100),
                'geometry': im.geometry(), 
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

        # SWM-Index based on Marta Milczarek et. al (2017)
        swm = im.expression('(blue+green)/(nir+swir1)', {
            'blue': im.select('SR_blue'),
            'green': im.select('SR_green'),
            'nir': im.select('SR_nir'),
            'swir1': im.select('SR_swir1')
            }).rename('swm')
        
        mndwi = im.expression('(green-swir1)/(green+swir1)', {
            'green': im.select('SR_green'),
            'swir1': im.select('SR_swir1')
            }).rename('mndwi')

        # 2.1. Identifying the categories ‘potentially water’ and ‘certainly water’
        wi = mndwi
        r_swir = im.select('SR_swir1')
        t_swir = 0.1
        t_wir = 0.4
        t_s =   _otsu(wi)
        wi_water = wi.updateMask(wi.gt(t_s))
        wi_land =  wi.updateMask(wi.lt(t_s))

        # 2.2 Refining the categories 'potentially water' and 'certainly water'
        m1 =    wi_water.reduceRegion(**{'reducer': ee.Reducer.median(), 'bestEffort':True}).values().getNumber(0)
        std1 =  wi_water.reduceRegion(**{'reducer': ee.Reducer.stdDev(), 'bestEffort':True}).values().getNumber(0)
        m2 =    wi_land.reduceRegion(**{'reducer': ee.Reducer.median(), 'bestEffort':True}).values().getNumber(0)
        std2 =  wi_land.reduceRegion(**{'reducer': ee.Reducer.stdDev(), 'bestEffort':True}).values().getNumber(0)
        t_pure =    ee.List([(m1.add(std1)).multiply(0.5), t_s]).reduce(ee.Reducer.max())
        t_mixed =   ee.List([(m2.add(std2)).multiply(0.5), t_s]).reduce(ee.Reducer.min())
        c1 =    wi.gt(ee.Image.constant(t_pure)).rename('c1')
        p1 =    (wi.gt(ee.Image.constant(t_mixed))).And(wi.lte(ee.Image.constant(t_pure))).rename('p1')
        p2 =    p1.And(r_swir.lt(t_swir)).rename('p2')
        wir =   wi.subtract(wi.reduceNeighborhood(**{'reducer': ee.Reducer.min(), 'kernel': ee.Kernel.square(**{'radius': 2})}))
        c2 =    c1.Or(p2.And(wir.gt(ee.Image.constant(t_wir)))).rename('c2')

        # 2.3 Extracting multi-scale water based on spatial analysis
        combined = c2.Or(p2)
        connected_components = combined.connectedComponents(**{
            'connectedness': ee.Kernel.square(1), #eight-connected
            'maxSize': 128
            }
        )
        index_criteria = connected_components.mask(wi).reduceConnectedComponents(
            reducer=ee.Reducer.mean(),
            labelBand='labels'
        ).gt(t_s)
        isWater = index_criteria.rename('is_water')
        buffer = apply_erosion(isWater, buffer_dist).rename('is_water')
        isWater = ee.Image(ee.Algorithms.If(
            buffer_dist,
            isWater.Or(buffer).rename('is_water'),
            isWater
        ))
        return(im.updateMask(wi_water).addBands(wi_water))
    return(_wrap)
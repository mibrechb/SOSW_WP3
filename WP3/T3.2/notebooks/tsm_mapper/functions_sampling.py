import ee

def sample_image(stations, sample_area=30*30):
  """Samples image bands at provided virtual stations and with the provided sample area."""
  def wrap(im):
    def _get_roi_coverage(geometry):
      """ Computes a property with relative area of clear pixels covering ROI """
      total_area = geometry.area(0.5) # calculate roi area in m^2
      unmasked_area = im.select('B4.*').mask() \
        .selfMask() \
        .multiply(ee.Image.pixelArea()) \
        .reduceRegion(**{
          'reducer': ee.Reducer.sum(), 
          'geometry': geometry, 
          'scale': scale,
          'maxPixels': 1e15,
          'bestEffort': True
        }).values().get(0)
      px_ratio = ee.Number(unmasked_area).divide(total_area)
      roi_coverage = px_ratio.multiply(100)
      return(roi_coverage)
      
    def _get_stats(station):
      """ Calculate areal stats over virtual station. """
      im_subset = im.select(['B.*', 'add.*'])
      geometry = station.geometry().buffer(ee.Number(sample_area).sqrt().divide(2), 0.01).bounds()
      station_id = station.getString('station_id')
      
      reducers = ee.Reducer.mean().combine(**{
        "reducer2": ee.Reducer.stdDev(), "sharedInputs": True}).combine(**{
        "reducer2": ee.Reducer.minMax(), "sharedInputs": True}).combine(**{
        "reducer2": ee.Reducer.median(**{"maxBuckets":500, "minBucketWidth": 0.125}), "sharedInputs": True})#.combine({
        #"reducer2": ee.Reducer.skew(), "sharedInputs":  true}).combine({
        #"reducer2": ee.Reducer.kurtosis(), "sharedInputs": True}).combine({
        #"reducer2": ee.Reducer.percentile({"percentiles": [05,25,75,95], "maxBuckets": 500, "minBucketWidth": 0.125}), "sharedInputs": true})
      
      cld_stats_500m = im.select('is_cloud').rename('cloudiness_500m').reduceRegion(**{
        'reducer': ee.Reducer.mean(), 'geometry': geometry.buffer(500), 'scale': scale, 'tileScale': 2
      })
    
      im_stats = im_subset.reduceRegion(**{
        'reducer': reducers, 'geometry': geometry, 'scale': scale, 'tileScale': 2
      })
      
      feature = ee.Feature(geometry, {}) \
        .set('system:time_start', im.get('system:time_start')) \
        .set('timestamp', ee.Date(im.get('system:time_start'))) \
        .set('station_id', station_id) \
        .set('roi_coverage', _get_roi_coverage(geometry)) \
        .set(ee.Dictionary(im_stats)) \
        .set(ee.Dictionary(cld_stats_500m))
      return(feature)
    
    CLOUD_COVER = ee.Number(ee.Algorithms.If(
      im.getNumber('CLOUDY_PIXEL_PERCENTAGE'),
      im.getNumber('CLOUDY_PIXEL_PERCENTAGE'),
      im.getNumber('CLOUD_COVER')
      ))
    platform = ee.String(ee.Algorithms.If(
      im.getString('SPACECRAFT_NAME'),
      im.getString('SPACECRAFT_NAME').slice(0, -1).toUpperCase(),
      im.getString('SPACECRAFT_ID').replace('_', '-')
      ))
    scale = ee.Number(ee.Algorithms.If(
      ee.String(platform).equals('SENTINEL-2'), 10, 30))

    fc_samples = stations \
      .map(_get_stats) \
      .map(lambda feature: feature.set('CLOUD_COVER', CLOUD_COVER).set('platform', platform))
    return(fc_samples)
  return(wrap)
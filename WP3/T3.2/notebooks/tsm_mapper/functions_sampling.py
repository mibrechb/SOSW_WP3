import ee

def get_matchups(fc_station, ic_rs, max_diff=1):
    """ Matches FeatureCollection with closest match from ImageCollection. """
    geometry = ee.FeatureCollection(fc_station).geometry()
    ic_rs = ic_rs.filter(ee.Filter.bounds(geometry))
    max_diff_filter = ee.Filter.maxDifference(**{
      'difference': max_diff * 24 * 60 * 60 * 1000,
      'leftField': 'system:time_start',
      'rightField': 'system:time_start'
    })
    save_best_join = ee.Join.saveBest(**{
      'matchKey': 'bestImage',
      'measureKey': 'timeDiff'
    })
    fc_matchups = save_best_join.apply(fc_station, ic_rs, max_diff_filter)
    return fc_matchups

def get_matchup_sample(buffer_dist=15):
    """ Sample matched image at feature geometry and add aggregated value as property. """
    def _wrap(feature):
      feature = ee.Feature(feature)
      match_img = ee.Image(feature.get('bestImage'))
      cloud_cover = ee.Algorithms.If(
        match_img.getString('sensor').equals('msi'),
        match_img.getNumber('CLOUDY_PIXEL_PERCENTAGE'),
        match_img.getNumber('CLOUD_COVER')
      )
      dt_img = match_img.date()
      dt_insitu = ee.Date.parse('YYYY-MM-dd HH:mm:ss', feature.get('dt_utc'))
      td_days = dt_img.difference(dt_insitu, 'day')
      geometry = feature.geometry().buffer(buffer_dist)

      reducers = ee.Reducer.mean() \
        .combine(**{"reducer2": ee.Reducer.stdDev(), "sharedInputs": True}) \
        .combine(**{"reducer2": ee.Reducer.minMax(), "sharedInputs": True}) \
        .combine(**{"reducer2": ee.Reducer.median(**{
           "maxBuckets":500, "minBucketWidth": 0.125}), "sharedInputs": True}) #\
        #.combine(**{"reducer2": ee.Reducer.skew(), "sharedInputs":  True}) \
        #.combine(**{"reducer2": ee.Reducer.kurtosis(), "sharedInputs": True}) \
        #.combine(**{"reducer2": ee.Reducer.percentile(
          #**{"percentiles": [5,25,75,95], "maxBuckets": 500, "minBucketWidth": 0.125}), "sharedInputs": True})
      samples_agg = match_img.reduceRegion(**{
         'reducer':   reducers, 
         'geometry':  geometry,
         'scale':     30, 
         'tileScale': 2
         })
      
      feature = feature \
          .set('dt_utc', dt_insitu) \
          .set('match_values', samples_agg) \
          .set('bestImage', match_img.get('system:index')) \
          .set('cloud_cover', cloud_cover) \
          .set('match_dt_utc', ee.Date(match_img.get('system:time_start')).format()) \
          .set('match_td_days', td_days) \
          .set('platform', match_img.getString('platform'))
      return(feature)
    return(_wrap)

def get_sample(geom, buffer_dist=15):
    """ Sample image as FeatureCollection at provided geometry and add aggregated value as property. """
    def _wrap(img):
      img = ee.Image(img)
      cloud_cover = ee.Algorithms.If(
        img.getString('sensor').equals('msi'),
        img.getNumber('CLOUDY_PIXEL_PERCENTAGE'),
        img.getNumber('CLOUD_COVER')
      )
      dt_utc_img = img.date()
      geometry = ee.Geometry(geom).buffer(buffer_dist)

      reducers = ee.Reducer.mean() \
        .combine(**{"reducer2": ee.Reducer.stdDev(), "sharedInputs": True}) \
        .combine(**{"reducer2": ee.Reducer.minMax(), "sharedInputs": True}) \
        .combine(**{"reducer2": ee.Reducer.median(**{
           "maxBuckets":500, "minBucketWidth": 0.125}), "sharedInputs": True}) #\
        #.combine(**{"reducer2": ee.Reducer.skew(), "sharedInputs":  True}) \
        #.combine(**{"reducer2": ee.Reducer.kurtosis(), "sharedInputs": True}) \
        #.combine(**{"reducer2": ee.Reducer.percentile(
          #**{"percentiles": [5,25,75,95], "maxBuckets": 500, "minBucketWidth": 0.125}), "sharedInputs": True})
      samples_agg = img.reduceRegion(**{
         'reducer':   reducers, 
         'geometry':  geometry,
         'scale':     30, 
         'tileScale': 2
         })
      
      feature = ee.Feature(None, {
        'dt_utc': dt_utc_img,
        'match_values': samples_agg,
        'cloud_cover': cloud_cover,
        'platform': img.getString('platform'),
        'sensor': img.getString('sensor')
        })
      return(feature)
    return(_wrap)
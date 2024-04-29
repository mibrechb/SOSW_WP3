import ee

def process_station(ic_rs, max_diff=1):
    def wrap(fc_station):
        fc_matchups = get_matchups(fc_station, ic_rs, max_diff).map(get_sample)
        return fc_matchups
    return wrap

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

def get_sample(feature):
    """ Sample matched image at feature geometry and add aggregated value as property. """
    feature = ee.Feature(feature)
    match_img = ee.Image(feature.get('bestImage'))
    time_diff = feature.get('timeDiff')
    geometry = feature.geometry().buffer(15)
    value = feature.get('value')
    samples_agg = match_img.reduceRegion(reducer=ee.Reducer.median(), geometry=geometry)
    feature = feature \
        .set('values_eo', samples_agg) \
        .set('bestmatch', match_img.get('system:index')) \
        .set('timediff', time_diff) \
        .set('sensor', match_img.getString('SPACECRAFT_ID')) \
        .set('timestamp', ee.Date(match_img.get('system:time_start')).format())
    return(feature)
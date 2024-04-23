# Band-LUT for Landsat TM (L4, L5) and ETM+ (L7)
bands_tm_etm = { 
    'blue':     '.*B1', 
    'green':    '.*B2', 
    'red':      '.*B3', 
    'nir':      '.*B4',  
    'swir1':    '.*B5',
    'swir2':    '.*B7', 
    'st':       '.*B6'
}

# Band-LUT for Landsat OLI (L8, L9)
bands_oli = { 
    'ultra_blue': '.*B1', 
    'blue':     '.*B2', 
    'green':    '.*B3', 
    'red':      '.*B4', 
    'nir':      '.*B5',  
    'swir1':    '.*B6',
    'swir2':    '.*B7',
    'st':       '.*B10' 
}
bands_msi = {
    'aerosols': 'B1', 
    'blue':     'B2', 
    'green':    'B3', 
    'red': 'B4', 
    're1': 'B5', 
    're2': 'B6', 
    're3': 'B7', 
    'nir': 'B8', 
    're4': 'B8A', 
    'water_vapour': 'B9', 
    'swir1': 'B11', 
    'swir2': 'B12'
}

bands_etm = { 
    'blue':     '.*B1', 
    'green':    '.*B2', 
    'red':      '.*B3', 
    'nir':      '.*B4',  
    'swir1':    '.*B5',
    'swir2':    '.*B7',  
}

bands_oli = { 
    'ultra_blue': '.*B1', 
    'blue':     '.*B2', 
    'green':    '.*B3', 
    'red':      '.*B4', 
    'nir':      '.*B5',  
    'swir1':    '.*B6',
    'swir2':    '.*B7', 
}
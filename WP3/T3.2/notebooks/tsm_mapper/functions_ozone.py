import ee

# Dataset selection #
ozone = ee.ImageCollection('TOMS/MERGED')
data = ozone.filter(ee.Filter.date('1999-01-01', '2030-01-01'))

# Linear interpolation over time #
# Add a band containing timestamp to each image
# This will be used to do pixel-wise interpolation later

def func_btl(image):
  timeImage = image.metadata('system:time_start').rename('timestamp')
  # The time image doesn't have a mask.
  # We set the mask of the time band to be the same as the first band of the image
  timeImageMasked = timeImage.updateMask(image.mask().select(0))
  return image.addBands(timeImageMasked)

data = data.map(func_btl)

# Specify the time-window
# This will determine how much backward and forward are we willing to
# look for an unmasked pixel in the time-series
days = 14

# For each image in the collection, we need to find all images
# before and after the specified time-window

# This is accomplished using Joins
# We need to do 2 joins
# Join 1: Join the collection with itself to find all images before each image
# Join 2: Join the collection with itself to find all images after each image

# We first define the filters needed for the join

# Define a maxDifference filter to find all images within the specified days
# The filter needs the time difference in milliseconds
# Convert days to milliseconds
millis = ee.Number(days).multiply(1000*60*60*24)
maxDiffFilter = ee.Filter.maxDifference(**{
  'difference': millis,
  'leftField': 'system:time_start',
  'rightField': 'system:time_start'
})

# We need a lessThanOrEquals filter to find all images after a given image
# This will compare the given image's timestamp against other images' timestamps
lessEqFilter = ee.Filter.lessThanOrEquals(**{
  'leftField': 'system:time_start',
  'rightField': 'system:time_start'
})

# We need a greaterThanOrEquals filter to find all images before a given image
# This will compare the given image's timestamp against other images' timestamps
greaterEqFilter = ee.Filter.greaterThanOrEquals(**{
  'leftField': 'system:time_start',
  'rightField': 'system:time_start'
})

# Apply the joins

# For the first join, we need to match all images that are after the given image.
# To do this we need to match 2 conditions
# 1. The resulting images must be within the specified time-window of target image
# 2. The target image's timestamp must be lesser than the timestamp of resulting images
# Combine two filters to match both these conditions
filter1 = ee.Filter.And(maxDiffFilter, lessEqFilter)
# This join will find all images after, sorted in descending order
# This will gives us images so that closest is last
join1 = ee.Join.saveAll(**{
  'matchesKey': 'after',
  'ordering': 'system:time_start',
  'ascending': False})

join1Result = join1.apply(**{
  'primary': data,
  'secondary': data,
  'condition': filter1
})
# Each image now as a property called 'after' containing
# all images that come after it within the time-window
#print(join1Result.first())

# Do the second join now to match all images within the time-window
# that come before each image
filter2 = ee.Filter.And(maxDiffFilter, greaterEqFilter)
# This join will find all images before, sorted in ascending order
# This will gives us images so that closest is last
join2 = ee.Join.saveAll(**{
  'matchesKey': 'before',
  'ordering': 'system:time_start',
  'ascending': True})

join2Result = join2.apply(**{
  'primary': join1Result,
  'secondary': join1Result,
  'condition': filter2
})

# Do the interpolation

# We now write a function that will be used to interpolate all images
# This function takes an image and replaces the masked pixels
# with the interpolated value from before and after images.

def interpolateImages(image):
  image = ee.Image(image)
  # We get the list of before and after images from the image property
  # Mosaic the images so we a before and after image with the closest unmasked pixel
  beforeImages = ee.List(image.get('before'))
  beforeMosaic = ee.ImageCollection.fromImages(beforeImages).mosaic()
  afterImages = ee.List(image.get('after'))
  afterMosaic = ee.ImageCollection.fromImages(afterImages).mosaic()

  # Interpolation formula
  # y = y1 + (y2-y1)*((t – t1) / (t2 – t1))
  # y = interpolated image
  # y1 = before image
  # y2 = after image
  # t = interpolation timestamp
  # t1 = before image timestamp
  # t2 = after image timestamp
  
  # We first compute the ratio (t – t1) / (t2 – t1)

  # Get image with before and after times
  t1 = beforeMosaic.select('timestamp').rename('t1')
  t2 = afterMosaic.select('timestamp').rename('t2')

  t = image.metadata('system:time_start').rename('t')

  timeImage = ee.Image.cat([t1, t2, t])

  timeRatio = timeImage.expression('(t - t1) / (t2 - t1)', {
    't': timeImage.select('t'),
    't1': timeImage.select('t1'),
    't2': timeImage.select('t2'),
  })
  # You can replace timeRatio with a constant value 0.5
  # if you wanted a simple average
  
  # Compute an image with the interpolated image y
  interpolated = beforeMosaic \
    .add((afterMosaic.subtract(beforeMosaic).multiply(timeRatio)))
  # Replace the masked pixels in the current image with the average value
  result = image.unmask(interpolated)
  return(result.copyProperties(image, ['system:time_start']))

# map() the function to interpolate all images in the collection
interpolatedCol = ee.ImageCollection(join2Result.map(interpolateImages))

# Show results
ozone_interp = interpolatedCol.select('ozone')
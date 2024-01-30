def apply_msi_scaling(image):
  """ Apply MSI L1C/L2A scaling factor. """
  opticalBands = image.select('B.*').divide(10000)
  return image.addBands(opticalBands, None, True)


def apply_oli_scaling(image):
  """ Apply OLI L2C2 scaling factors and offsets. """
  opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, None, True).addBands(thermalBands, None, True)
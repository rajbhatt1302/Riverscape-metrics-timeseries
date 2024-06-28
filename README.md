# Riverscape-metrics-timeseries

**Floodplain segmentation**
Floodplain segmentation into DGOs was implemented in a Python environment in Google Colab. The code file is Floodplain_segmentation_into_DGOs.ipynb

**Riverscape metrics time series generation**
Google Earth Engine has implemented riverscape metrics time series generation of pixel count and water mask perimeter for each date. The code has been provided in the timeseries.txt file. 
Later, the time series data was exported to .csv files where the pixel count was multiplied with LANDSAT's pixel area to compute the river metric(s) area for each date. Then, the river metric(s) area was divided by the length of each DGO to derive the river metric(s) width for that particular DGO. Further, the Braiding Index for each date was computed by dividing the perimeter of the water mask with 2x segment length of each DGO. 



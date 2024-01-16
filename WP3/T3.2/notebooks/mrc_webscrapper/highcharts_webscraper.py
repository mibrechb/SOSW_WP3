# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:05:16 2023

@author: Michael
"""

import pandas as pd

from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium import webdriver
from selenium.webdriver import Chrome, ChromeOptions
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities

def get_local_safe_setup():
    options = ChromeOptions() 
    options.add_argument("--disable-blink-features")
    options.add_argument("--disable-blink-features=AutomationControlled")
    options.add_argument("--disable-infobars")
    options.add_argument("--disable-popup-blocking")
    options.add_argument("--disable-notifications")

    driver = Chrome(desired_capabilities = options.to_capabilities())

    return driver

if __name__ == "__main__":
    # Create local setup like here: https://gist.github.com/theDestI/aa21a0e721b06a74bd58a0a391d96e8f
    driver = get_local_safe_setup() 
    
    # Define url
    url = 'https://portal.mrcmekong.org/time-series/chart?ts=0a6a205f58eb436787a0ebd26bab8fcb&sd=1987-09-14&ed=2003-12-13'
    
	# Call url through selenium driver
    driver.get(url)

# 	# Wait until element with highcharts graph appears
#     WebDriverWait(driver, 10).until(
#         EC.presence_of_element_located((By.CLASS_NAME, "app-time-series-chart"))
#     )

	# Parse dates and values in those dates
    dates = driver.execute_script('return $("highcharts-chart").highcharts().series[0].data.map(x => x.series).map(x => x.xData)[0].map(x => new Date(x).toISOString())')
    values = driver.execute_script('return $("highcharts-chart").highcharts().series[0].data.map(x => x.series).map(x => x.yData)[0]')

    df = pd.DataFrame({'Dates': dates, 'Temperatures': values })

	# Store to dataframe and save it in current folder
    #df.to_csv('./output.csv', index=False)
    print(df)
    
	# Don't forget ot quit!
    driver.quit()
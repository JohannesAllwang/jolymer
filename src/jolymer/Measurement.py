# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:14:39 2020

@author: xcill
"""

from pathlib import Path
import os
import getpass
import datetime as dt

path = Path(__file__).parent

class Measurement:
    
    rawdatapath = os.path.join(path, '../../rawdata')
    rawdatapath = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\rawdata\\"
    figures_path = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\figures\\"
    
    def get_mdt(self):
        "returns date in datetime format and boolean for if time is includet or not"
        year = int(self.datestring[0:4])
        month = int(self.datestring[4:6])
        day = int(self.datestring[6:8])
        
        if self.timestring!= None:
            hour, minute = [int(x) for x in self.timestring.split(':')]
            date = dt.datetime(year, month, day, hour=hour, minute=minute)
            return date, True
        else:
            date = df.date(year, month, day)
            return date, False

    def get_sample_age(self):
        pdt = self.sample.get_pdt()[0]
        mdt = self.get_mdt()[0]
        age = mdt - pdt

        return age

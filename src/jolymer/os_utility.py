#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:15:52 2020

@author: johannes
"""


import os
import sys
import zipfile
import shutil
import subprocess


def unzip(source_filename, dest_dir):
    with zipfile.ZipFile(source_filename) as zf:
        for member in zf.infolist():
            words = member.filename.split('/')
            path = dest_dir
            zf.extract(member, path)


def create_path(path):
    try:
        os.mkdir(path)
    except OSError:
        # print ("Creation of the directory %s failed" % path)
        pass
    else:
        # print ("Successfully created the directory %s " % path)
        pass


def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    if sys.platform == 'linux':
        cmd = [f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=atsas_bin)

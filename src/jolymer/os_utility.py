# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:15:52 2020

@author: johannes
"""

from configparser import ConfigParser

import os
import sys
import zipfile
import shutil
import subprocess
import logging

config = ConfigParser()

# config["DEFAULT"] = {
#     'datapath': '/home/johannes/data'
#     }

config.read(os.path.expanduser('~/.config/jolymer/jolymer.ini'))
# config.get('measurement_path')
config_data = config[ 'johannes' ]
print(config_data.get('measurement_path'))

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(
        os.path.expanduser('~/.local/share/jolymer/jolymer.log'))
formatter = logging.Formatter("%(levelname)s %(message)s")
file_handler.setFormatter(formatter)
# file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

# logging.basicConfig(
#         level=logging.DEBUG,
#         format="%(levelname)s %(message)s",
#         filename=os.path.expanduser(
#             '~/.local/share/jolymer/jolymer.log')
#         )

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
        logging.info("Creation of the directory %s failed" % path)
        pass
    else:
        logging.info("Successfully created the directory %s " % path)
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

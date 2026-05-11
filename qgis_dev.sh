#!/bin/bash

source /home/lemap/umep/bin/activate
/home/lemap/umep/bin/python3 -m pip install pip -U --break-system-packages

pb_tool deploy -y

# On force PYTHONPATH à pointer d'abord sur le venv
export PYTHONPATH=/home/lemap/umep/lib/python3.14/site-packages:/usr/share/qgis/python:/usr/share/qgis/python/plugins:$PYTHONPATH

export QGIS_PREFIX_PATH=/usr
export GDAL_FILENAME_IS_UTF8=YES

qgis
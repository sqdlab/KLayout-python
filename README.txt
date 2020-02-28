This project contains code developed to design superconducting devices. There are two versions of the code, one made to work with the KLayout Python interpreter, and one that runs in a Python 3 installation.

## Installation

### KLayout interpreter

To successfully use this library you should:
1. download this library from github/shamil777/KLayout-python
2. open KLayout, press 'Ctrl + F5'(Windows) 'F5'(MacOS) to open macros editor.
3. In left panel choose the "Python" tab.
4. Right-click to the left panel and choose 'Add Location'. Enter path to the 'KLayout-python' directory.
5. Create KLAYOUT_PYTHONPATH environment variable (this variable is used by KLayout python interpreter). 
	This step depends on your OS so google 'how to' do this step.
7. Locate and launch KLAYOUT_PYTHONPATH.py. Assign output to the KLAYOUT_PYTHONPATH environment variable.
8. Restart KLayout. 
9. Open macro editor and try to launch some example from "Klayout-python/Examples" folder.
10. Congrats! Have fun. If not -> leave an issue.

### Python package

1. Install the KLayout Python module from https://github.com/KLayout/klayout/wiki/KLayout-Python-Module

2. Now you can use the ClassLib module.

## Tutorials

A very basic tutorial is found in the Tutorials folder as Python Notebook files.
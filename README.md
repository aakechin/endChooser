# endChooser
Tool endChooser chooses the best variant of the last two nucleotides of allele-specific primer

## Requirements
### Linux
endChooser works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython`
* **argparse** - you can install it with: `sudo apt-get install python3-argparse`
* **easygui** - for use of GUI-version: `sudo apt-get install python3-easygui`

### Windows
For use on windows download and install python3.4.3 from www.python.org/downloads/release/python-343/ (remember to check "Add python to PATH" during installation). 
If you do not have Visual Studio C++ already installed, download and install it from http://landinghub.visualstudio.com/visual-cpp-build-tools/ or from https://www.microsoft.com/ru-ru/download/details.aspx?id=48145. 
Visual Studio C++ may ask you to update .Net Framework. Download and install it from https://www.microsoft.com/ru-ru/download/details.aspx?id=39257 (generally for Windows XP users). 
After that open command line and install the followong packages with respective commands:
* **Biopython** - with: `pip install biopython`
* **argparse** - you can install it with: `pip install argparse`
* **easygui** - for use of GUI-version: `pip install easygui`

### Mac OS
For use on Mac OS download and install python3.6 from www.python.org/downloads/
After that install the followong packages with respective commands in command line:
* **Biopython** - with: `pip install biopython`
* **argparse** - you can install it with: `pip install argparse`
* **easygui** - for use of GUI-version: `pip install easygui`

## Installation
endChooser does not require any installations

## Use
You can see all parameters for command-line version with:
`python endChooser.py -h`
Or you can start GUI-version with:
`python endChooser_gui.py`
## Parameters
```
-h, --help            show this help message and exit
  --fasta-file FASTAFILE, -ff FASTAFILE
                        multi-fasta-file with highlighted variable position in
                        the format of [A/T], where A is a reference allele and
                        T is an alternative
```
## Citation
Manuscript is prepared. Now you can cite it by link.

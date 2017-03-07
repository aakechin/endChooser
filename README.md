# endChooser
Tool endChooser chooses the best variant of the last two nucleotides of allele-specific primer

## Requirements
### Linux
endChooser works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython`
* **argparse** - you can install it with: `sudo apt-get install python3-argparse`

### Windows
For use on windows download and install python3.5 from www.python.org/downloads/release/python350/ (remember to check "Add python to PATH")
If you do not have Visual Studio C++ already installed, download and install it from landinghub.visualstudio.com/visual-cpp-build-tools/
After that install the followong packages with respective commands in command line:
* **Biopython** - with: `pip install biopython`
* **argparse** - you can install it with: `pip install argparse`

### Mac OS
For use on Mac OS download and install python3.6 from www.python.org/downloads/
After that install the followong packages with respective commands in command line:
* **Biopython** - with: `pip install biopython`
* **argparse** - you can install it with: `pip install argparse`

## Installation
endChooser does not require any installations

## Use
You can see all parameters with 
```
python3 endChooser.py -h
```
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

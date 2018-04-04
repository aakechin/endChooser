# endChooser
Tool endChooser chooses the best variant of the mismatch amplification mutation assay (MAMA) allele-specific primers.

## Requirements
### Linux
endChooser works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython`
* **primer3** - you can install it with: `sudo apt-get install python3-primer3`
* **easygui** - for use of GUI-version: `sudo apt-get install python3-easygui`

### Windows
Unfortunately, primer3 does not work on Windows, so endChooser can not be used under Windows, too.

### Mac OS
For use on Mac OS download and install python3.6 from www.python.org/downloads/
After that install the followong packages with respective commands in command line:
* **Biopython** - with: `pip install biopython`
* **primer3** - you can install it with: `pip install primer3-py`
* **easygui** - for use of GUI-version: `pip install easygui`

## Installation
endChooser does not require any installations

## Use
Go to the directory where program was installed with:
`cd endChooser-master`
You can start GUI-version with:
`python endChooser_gui.py`
After that, choose fasta-file with sequence for which you want to design allele-specific-primers. It should have the format like the example file "seq.fa". Variable position should be designated like [A/T] or (A/C). When you select a file,click Open. You will get a new windows for choosing parameters of primer design. When you fill in values, click OK, and process of primer design will start.
## Citation
Manuscript is prepared. Now you can cite it by link.

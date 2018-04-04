# endChooser
Tool endChooser chooses the best variant of the mismatch amplification mutation assay (MAMA) allele-specific primers.

## Requirements
### Linux
endChooser works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython`
* **primer3** - you can install it with: `sudo apt-get install python3-primer3`
* **argparse** - for use command-line parameters: `sudo apt-get install python3-argparse`
* **myvariant** - for checking if primers cover SNPs from dbSNP: `sudo apt-get install python3-myvariant`
* **xlsxwriter** - for writing output to EXCEL-file: `sudo apt-get install python3-xlsxwriter`

### Windows
For use endChooser in Windows, you need to install Cygwin (https://cygwin.com/install.html). When you download the install package, run it. 1) Press "Next"; 2) Choose "Install from Internet" (default), press "Next"; 3) Choose "Root Directory" (you can leave it default path), press "Next"; 4) Choose "Local Package Directory" (you can leave it default path), press "Next"; 5) Choose your proxy settings (in the most of cases, it is default "Use System Proxy Settings"); 6) Wait, while installation gets full list of Download Sites; 7) Choose preferable Download Site (e.g. you can choose the first one), press "Next"; 8) Wait, while installation gets full list of additional packages; 9) In the selection list "View" choose "Not installed" and in the "Search" field write "python3". Find in the list package "python3: Py3K language interpreter" and click on the value of column "New". This will select "Python3" for installation. After that, select in the similar way package with name "make: The GNU version of the 'make' utility". Do not forget to write "make" in the Search field. Press "Next" and "Next". Wait, while installation is completed.
Then, run istalled Cygwin and install additional python packages with the following commands in the command line:
* **Biopython**: `pip3 install biopython`
* **primer3**: `pip3 install primer3-py`
* **argparse**: `pip3 install argparse`
* **myvariant**: `pip3 install myvariant`
* **xlsxwriter**: `pip3 install xlsxwriter`

### Mac OS
For use on Mac OS download and install python3.6 from www.python.org/downloads/
After that install the followong packages with respective commands in command line:
* **Biopython**: `pip3 install biopython`
* **primer3**: `pip3 install primer3-py`
* **argparse**: `pip3 install argparse`
* **myvariant**: `pip3 install myvariant`
* **xlsxwriter**: `pip3 install xlsxwriter`

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

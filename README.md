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
For use endChooser in Windows, you need to install Cygwin (https://cygwin.com/install.html). When you download the install package, run it. 1) Press "Next"; 2) Choose "Install from Internet" (default), press "Next"; 3) Choose "Root Directory" (you can leave it default path), press "Next"; 4) Choose "Local Package Directory" (you can leave it default path), press "Next"; 5) Choose your proxy settings (in the most of cases, it is default "Use System Proxy Settings"); 6) Wait, while installation gets full list of Download Sites; 7) Choose preferable Download Site (e.g. you can choose the first one), press "Next"; 8) Wait, while installation gets full list of additional packages; 9) In the selection list "View" choose "Category" and in the "Search" field write "python". Find in the list category "Python" and click on the "Default", so it will be changed onto "Install". This will select all Python packages for installation. After that, change "View" field from "Category" onto "Not installed" and write "make" in the Seach field. Select packages with names "make: The GNU version of the 'make' utility". Press "Next" and "Next". Wait, while installation is completed. In the end you can add icons for Cygwin on the desktop.  
Then, run istalled Cygwin and install additional python packages with the following commands in the command line:
* **Biopython**: `pip3 install biopython`
* **primer3**: `pip3 install primer3-py`
* **argparse**: `pip3 install argparse`
* **myvariant**: `pip3 install myvariant`
* **xlsxwriter**: `pip3 install xlsxwriter`  
#### Use of endChooser in Cygwin command line
Download endChooser package to some directory (button "Clone or download", "Download ZIP"). It is simpler to download it to the home directory of Cygwin folder (e.g., for me, it is C:\cygwin64\home\Andrey\). Unzip it to the current directory. Now, you should start Cygwin (if you have stopped it earlier) and go to the directory with unzipped endChooser. For this use the following commands:  
`cd endChooser-master` (**Tip**: while you are typing some command or file/directory name in the command line, you can press "Tab" button on the keyboard and it will be completed automatically. If there are several commands or file/directory name that fit written letters, press "Tab" two times, and it will write all possible variants).  
If you have downloaded endChooser to some windows disk (e.g., C:\Users\Andrey\Downloads\), you will need to go to this directory (after unzipping endChooser archive):  
`cd /cygdrive/c/Users/Andrey/Downloads/`  
Now you can run endChooser (this will start endChooser with the example file):  
`python3 endChooser.py -fa seq.fa -alen 150 -alendev 10 -plen 20 -plendev 4 -mtemp 60 -mtempdev 5 -dg 3000`  
#### endChooser arguments and parameters
Some parameters have default values and if they fit your needs, you can skip them. All list of parameters is below:  
|Parameter | Description |
| --- | --- |
| -fa | Fasta-file that contains sequence with variable position designated like [A/G] or (A/G). You can see example file seq.fa |
| -alen | Amplicon length (Default: 150 bp). This value includes primers' lengths |
| -alendev | Amplicon length deviation (Default: 10 bp). So length of primer will be +- this value |
| -plen | Optimal (in your opinion) primer length (Default: 20 bp) |
| -plendev | Optimal primer length deviation (Default: 4 bp) |
| -mtemp | Optimal melting temperature (Default: 60 degrees Celsius) |
| -mtempdev | Optimal melting temperature deviation (Default: 5 degrees Celsius) |
| -dg | Minimal dG of dimer and hairpin formation (Default: 3000 kcal/mol) that you can accept |
| -min | Use this parameter, if you want amplicons to have possibly minimal length |
| -bd | Use this parameter, if blast search of your amplicons has been already done by endChooser (e.g. you have already started endChooser without this argument, and now you want to change some parameters and run endChooser with the same fasta-file again. So endChooser will take result of previous blast of input sequences) |
| -ref | Fasta-file with human reference genome sequence. Use this parameter only if you want to check primers for covering SNPs from dbSNP |

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

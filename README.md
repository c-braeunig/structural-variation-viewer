# Structural Variation Viewer

A suite of Python-based tools to visualize genomic variation, with focus on structural variants. 
The suite consists of a number of command-line scripts as well as a graphical user interface for visual data exploration. 
__SVread.py__ and __SVsee.py__ extract variant data from VCF files and visualize them respectively. __SVsee.py__ can be display additional annotation data corresponding to the used VCF file and a heat map with selectable window size.
Both functionalities are boundled in __SVpipe.py__ to allow for a one-step process between data and figure.
__SVpipe.py__'s functionality is extended with a GUI in SVgui.py. 

## Prerequisites
The SVx.py programs require*:
- __Python 3.7.6__
- __Matplotlib 3.2.1__
- __Tkinter 8.6__
- __Numpy 1.18.3__

The variant mapping script requires*:
- __samtools 1.9__
- __bedtools 2.29.2__
- __ngmlr 0.2.7__
- __sniffles 1.0.11__

\* *The scripts may work also with older version of the respective programs*

The [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installers provide quick and easy installation instructions and bring the other prerequisites.
Alternatively, you can install Python 3 with your preferred method and then additionally download the other prerequisites.

Once all the required software is in place, simply download the program(s) of interest.

# Usage
When relying on the command-line programs, you can choose to move step-wise from input VCF to output figure (using __SVread.py__ and __SVsee.py__) or to take the direct path (using __SVpipe.py__). The former allows for additional manipulation before visualization, the latter is more intuitive and quicker.

### Variant mapping
The following programs were developed in a usage context of existing assemblies that were to be investigated for structural variants. The complete mapping process can be handled at once using the enclosed variant mapping script __full_SV_pipe.sh__. The bash script first uses samtools's _faidx_ to make an index of the 'query' assembly (in FASTA format) that is later to be aligned against the reference assembly. Based on the query index, __makebed.py__ generates a BED file of simulated sequencing reads by designating distinct sequence section coordinates based on the window size and step size parameters given. BEDtools's _getfasta_ then extracts these fake reads from the query assembly. _ngmlr_ then maps these reads against the reference assembly. Based on the alignments, _sniffles_ calls the variants, producing a VCF-formatted file.

### Using __SVread.py__ and __SVsee.py__
In the location where the programs are saved execute 
```
$ python SVread.py -h
```
This will print all available options and settable parameters
```
usage: SVread.py [-h] [--vcf VCF]
                 [--sort {asc_size,des_size,type,type_asc_size,type_des_size,CR}]
                 [--ov] [-C CONTIGNAME] [-T TYPE [TYPE ...]]
                 [-S CONTTYPE [CONTTYPE ...]] [--ulim ULIM] [--llim LLIM]
                 [--seq] [--nosize] [--circos] [--out OUT]

Parse and derive SV information from .vcf output of variant calling software;
--ov generates a superficial overview of the listed SV types and their number
of occurrences. Using --sort, --C, --T, --S, --llim and --ulim, the SVs of
intereste can be isolated and printed to stdout or to a file ( --out) in a
bed-like file format. For further information on the vcf format consult
https://samtools.github.io/hts-specs/VCFv4.2.pdf

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             Provide variant call file in VCF format
  --sort {asc_size,des_size,type,type_asc_size,type_des_size,CR}
                        Provide sorting criterium to organize SVs; SVs sorted
                        by position on chromosome/contig by default; no
                        sorting recommended for downstream SVsee use; Use "CR"
                        when filtering by "BND" to sort target locations
  --ov                  Choose to display an overview of SVs in the provided
                        vcf file; --ov only will print assembly-level
                        information; using --C option will print single-
                        contig-level information
  -C CONTIGNAME, --contigname CONTIGNAME
                        Provide contig/chromosome ID by which to filter SV
                        records; multiple IDs given space-separated; use only
                        one for visualization
  -T TYPE [TYPE ...], --type TYPE [TYPE ...]
                        Provide type(s) by which to filter SV records i.e.
                        'DEL'; multiple types given space-separated; 'BND' for
                        complex rearrangments with breakends
  -S CONTTYPE [CONTTYPE ...], --conttype CONTTYPE [CONTTYPE ...]
                        Complementary type filter, same usage as --T
  --ulim ULIM           Provide upper length limit for SVs [bp] (inclusive)
  --llim LLIM           Provide lower length limit for SVs [bp] (inclusive)
  --seq                 Select to write sequence data into single fasta file;
                        requires --out
  --nosize              Select to not display the end position and size [bp]
                        of each SV in line
  --circos              Circos mode; generates karyotype and link files in
                        current directory for circos for BND-type SVs to show
                        inter-chromosomal events; compatible w/ -C option
  --out OUT             Provide desired name of bed-formatted output file
```
__SVread.py__ relies on a VCF file and can filter and sort the variant calls. If you want to gain a rough overview of the variant call data execute ```python SVread.py --ov```. This will print the IDs of the contigs in the assembly and the numbers of variants called on each contig. The breadth of the overview can be narrowed down to a given contig with ```python SVread.py --ov --C [contig ID]```. __SVread.py__ can print overviews for several contigs at once, but visualization only works with individual contigs. Several filters and sorting options allow to only visualize variant subsets of interest. Option order is not relevant. Save the filtered variant calls to an output file with ```--out```.

The filtered variant calls can then be fed to __SVsee.py__ for visualization.
To get an idea of the available options and parameters execute
```
$ python SVsee.py -h
```
This will print
```
usage: SVsee.py [-h] [--IN IN] [--vcf VCF] [--gff MOD] [--bed MOD]
                [--blast MOD] [-R ANNO_TYPE [ANNO_TYPE ...]]
                [-RC ANNO_CONT [ANNO_CONT ...]] [-RL ANNO_LOWER]
                [-RU ANNO_UPPER] [--heat HEAT] [--reg REG REG] [--out OUT]
                [--outfmt {pdf,svg,png}]

Single-chromsome visualization of SV distribution based on SVread.py output

optional arguments:
  -h, --help            show this help message and exit
  --IN IN               Provide input file (SVread.py output file)
  --vcf VCF             Provide VCF file corresponding to input file
  --gff MOD             Provide contig annotation as gff file; must have .gff
                        or .GFF extension
  --bed MOD             Provide contig annotation as bed file; must have .bed
                        or .BED extension
  --blast MOD           Provide contig annotation as blast output file (outfmt 6)
  -R ANNOTYPE [ANNOTYPE ...], --annotype ANNOTYPE [ANNOTYPE ...]
                        Provide type(s) of annotations to be included 
                        (only applicable on GFF annotation files)
  -RC ANNO_CONT [ANNO_CONT ...], --anno_cont ANNO_CONT [ANNO_CONT ...]
                        Complementary type filter, selected type(s) will be
                        excluded, multiple arguments given space-separated
                        (only applicable on GFF annotation files)
  -RL ANNO_LOWER, --anno_lower ANNO_LOWER
                        Provide lower length limit for annotation features (in bp)
  -RU ANNO_UPPER, --anno_upper ANNO_UPPER
                        Provide upper length limit for annotation features (in bp)
  --heat HEAT           Generate heatmap with given window size (in bp)
  --reg REG REG         Specify region bounds (in bp) for more zoomed-in look;
                        bounds given whitespace-separated i.e. 1000 1500
  --out OUT             Provide desired name of output file
  --outfmt {pdf,svg,png}
                        Select output format, default=png
```
After providing __SVread.py__'s output file and the VCF file again, you can additionally provide annotations in the form of GFF, BED or (tabular) output (w/ ```--gff```, ```--bed```, ```--blast``` respectively) as well as a heat map to visualize variant density more intiutively (```--heat```). The corresponding graphs are added below the variant graph within the same figure. The view of the selected contig can be restricted with ```--reg```. Finally save the figure with ```--out``` and ```-outfmt```.

### Using SVpipe.py
__SVpipe.py__ combines the functionalities of __SVread.py__ and __SVsee.py__ into a single script, fed with a VCF file and immediately producing a figure. 
Execute to see all options and parameters
```
$ python SVpipe.py -h
```
```
usage: SVpipe.py [-h] [--vcf VCF] [-C CONTIGNAME] [-T TYPE [TYPE ...]]
                 [-S CONTTYPE [CONTTYPE ...]] [--ulim ULIM] [--llim LLIM]
                 [--gff MOD] [--bed MOD] [--blast MOD]
                 [-R ANNO_TYPE [ANNO_TYPE ...]]
                 [-RC ANNO_CONT [ANNO_CONT ...]] [-RL ANNO_LOWER]
                 [-RU ANNO_UPPER] [--heat HEAT] [--reg REG REG] [--out OUT]
                 [--outfmt {pdf,svg,png}] [--ov]

Full structural variant visualization from raw VCF file to figure

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             Provide VCF file
  -C CONTIGNAME, --contigname CONTIGNAME
                        Provide contig/chromosome criterion by which to filter
                        SV records; multiple criteria given space-separated
  -T TYPE [TYPE ...], --type TYPE [TYPE ...]
                        Provide type criterion by which to filter SV records
                        i.e. 'DEL'; multiple criteria given space-separated;
                        'BND' for complex rearrangments with breakends
  -S CONTTYPE [CONTTYPE ...], --conttype CONTTYPE [CONTTYPE ...]
                        Verbose type filter: types entered will be omitted.
                        Same usage as --tfil
  --ulim ULIM           Provide upper length limit for SVs [bp] (inclusive)
  --llim LLIM           Provide lower length limit for SVs [bp] (inclusive)
  --gff MOD             Provide contig annotation as gff file; must have .gff
                        or .GFF extension
  --bed MOD             Provide contig annotation as bed file; must have .bed
                        or .BED extension
  --blast MOD           Provide contig annotation as blast output file (outfmt
                        6)
  -R ANNOTYPE [ANNOTYPE ...], --annotype ANNOTYPE [ANNOTYPE ...]
                        Provide type(s) of annotations to be included 
                        (only applicable on GFF annotation files)
  -RC ANNO_CONT [ANNO_CONT ...], --anno_cont ANNO_CONT [ANNO_CONT ...]
                        Complementary type filter...analogous to -S
                        (only applicable on GFF annotation files)
  -RL ANNO_LOWER, --anno_lower ANNO_LOWER
                        Provide lower length limit for annotation features (in bp)
  -RU ANNO_UPPER, --anno_upper ANNO_UPPER
                        Provide upper length limit for annotation features (in bp)
  --heat HEAT           Generate heatmap with given window size (in bp)
  --reg REG REG         Specify region bounds for more zoomed-in look; bounds
                        given whitespace-separated i.e. 1000 1500
  --out OUT             Provide desired name of graphic output
  --outfmt {pdf,svg,png}
                        Select output format, default=png
  --ov                  Choose to display an overview of SVs in the provided
                        VCF file; --ov only will print assembly-level
                        information; using --C option will print single-
                        contig-level information
``` 
__SVpipe.py__ features all relevant options and flags from __SVread.py__ and __SVsee.py__.

For example
```
$ python SVpipe.py --vcf OEVEvOLAT_genomic.vcf \
                   --gff OLAT_genomic.gff \
                   -C NC_019867.2
                   -T DEL INS \
                   --llim 100 \
                   --ulim 10000 \
                   -R exon CDS mRNA \
                   -RL 100 \
                   -RL 5000 \
                   --heat 1000000 \
                   --out OEVEvOLAT.png
```
will render deletions and insertions between 100 and 10000 bp in length on contig NC_019867.2 from the file _OEVEvOLAT_genomic.vcf_ and all exon, CDS and mRNA features from _OLAT\_genomic.gff_ between 100 and 5000 bp in length, along with a SV density heatmap with window size of 1 mb. 

![Example SVpipe.py render](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/OEVEvOLAT.png?raw=true)

### GUI-based exploration with SVgui.py
While they get the job done, the scripts above don't let look you at the graphs before saving them and as a result fine tuning parameters can be a bit laborious. __SVgui.py__ allows you to look at the data, pan, zoom and view a number of contigs at once and only save the figures you want. For this, it makes use of a tkinter-based graphical user interface and as a result requires a terminal with graphical capabilities. For testing, I used Anaconda3's Anaconda prompt locally. 

```$ python SVgui.py ``` will start SVgui.py and open the main window.

![Main Interface Window](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_main_2.JPG?raw=true)

The main interface window features a textbox for showing and recording information about the used files and shown data (yellow), a parameter pane (blue) and a tip pane (green).
Load a VCF file by clicking ```Load```, followed by ```Variant file```. Then select a contig from the dropdown menu at the top of the parameter pane. With  ```Load``` and ```Annotation file``` a GFF, BED or blast (tabular) output file can be loaded for annotation. This can also be done later on though.

![Contig Dropdown Menu](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_contig.JPG)

Once selected an overview of the SV types and numbers on the selected contig will be printed to the textbox and the various parameters option interfaces updated. Hover over individual option interfaces to receive information in the tip pane. 

![After contig selection](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_double_typeselection.JPG?raw=true)

Now you can set the same parameters as seen in __SVpipe.py__.

![Parameter selection](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_both_selected.JPG?raw=true)

Finally, clicking ```Draw``` will trigger the spawning of a child window and generation of the figure therein.

![Figure window](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_figure_window.JPG?raw=true)

When a GFF-formatted annotation file was loaded, left-click on individual features will spawn a child window that records the attributes of the clicked features and allows for saving the attributes to file.

![Feature attribute window](https://github.com/c-braeunig/structural-variation-viewer/blob/main/images/GUI_anno_info_box.JPG?raw=true)

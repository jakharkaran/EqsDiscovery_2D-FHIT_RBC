# Equation discovery for 2D-FHIT and RBC

Discover subgrid-scale (SGS) momentum flux closures for 2D forced homogeneous isotropic turbulence (FHIT) and 2D turbulent Rayleigh-Bénard convection (RBC), and SGS heat flux closures for RBC.

## Installation

If you are interested in using this, please clone this repository and run,
```
pip install -e ./
```

#### [[project website]](http://pedram.rice.edu/team/)
<img src="docs/repo_template.png" width="250">

## Table of contents
* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Experiments](#Experiments)
    * [Case 1](#Case-1)
* [Citation](#Citation)
* [References](#References)

## Introduction
Here, we apply a common equation-discovery technique with expansive libraries to learn closures for SGS fluxes from filtered direct numerical simulations (FDNS) of 2D-FHIT and RBC across several common filters types (Gaussian, box, Gaussian + box, sharp-spectral) and filter sizes. 


<ul>
<li>Do the first item,</li>
<li>Test the second item,</li>
<li>Investigate the third item,</li>
<li>Code the fourth item.</li>
</ul>

## Requirements
<!-- These are examples,
	add or remove as appropriate -->

- Matlab R2016+
- python 3.6
	- [scipy](https://pypi.org/project/scipy/)
	- [numpy](https://pypi.org/project/numpy/)
- [TensorFlow 2](https://www.tensorflow.org/install)
- [Keras 2.3.1](https://pypi.org/project/Keras/)

## Experiments
### Case 1
Case 1 is disscused here [Case 1 Location](./experiments/case1) 

open matlab
```
matlab -nodisplay -nosplash
```

Run the main file
```
python main_example.py
```

Post process
```
python post_example.py
```

Python code

```python
def myfun():
   print('Hello!')
```

Latex 

```bibtex
@article { Lubis_AMS_2021,
      author = {Sandro W. Lubis and Pedram Hassanzadeh},
      title = {An Eddy–Zonal Flow Feedback Model for Propagating Annular Modes},
      journal = {Journal of the Atmospheric Sciences},
      year = {2021},
      publisher = {American Meteorological Society},
      address = {Boston MA, USA},
      volume = {78},
      number = {1},
      doi = {10.1175/JAS-D-20-0214.1},
      pages= {249 - 267},
      url = "https://journals.ametsoc.org/view/journals/atsc/78/1/jas-d-20-0214.1.xml"
}
```

## Citation
- Lubis, Sandro W., and Pedram Hassanzadeh. " An Eddy–Zonal Flow Feedback Model for Propagating Annular Modes", Journal of the Atmospheric Sciences 78, 1 (2021).([url](https://doi.org/10.1175/JAS-D-20-0214.1))<details><summary>BibTeX</summary><pre>
@article { Lubis_AMS_2021,
      author = {Sandro W. Lubis and Pedram Hassanzadeh},
      title = {An Eddy–Zonal Flow Feedback Model for Propagating Annular Modes},
      journal = {Journal of the Atmospheric Sciences},
      year = {2021},
      publisher = {American Meteorological Society},
      address = {Boston MA, USA},
      volume = {78},
      number = {1},
      doi = {10.1175/JAS-D-20-0214.1},
      pages= {249 - 267},
      url = "https://journals.ametsoc.org/view/journals/atsc/78/1/jas-d-20-0214.1.xml"
}</pre></details>

## References
For a guide to markdown syntax see  

```
https://www.markdownguide.org/basic-syntax/
```



MGCAMB v4.0
===========
## Modified Growth with CAMB
This is the official repository for the MGCAMB v4.0 patch.  Below is an introduction and the instructions to install and run the code.  

## Table of contents
- [MGCAMB v4.0](#mgcamb-v40)
  - [Modified Growth with CAMB](#modified-growth-with-camb)
  - [Table of contents](#table-of-contents)
  - [1. Introduction](#1-introduction)
    - [Structure of the code](#structure-of-the-code)
    - [Citing MGCAMB](#citing-mgcamb)
  - [2. How to install](#2-how-to-install)
  - [3. How to run](#3-how-to-run)
  - [4. What's new](#4-whats-new)
  - [5. Known issues](#5-known-issues)
  - [6. Authors List](#6-authors-list)



## 1. Introduction
Modified Growth with CAMB (MGCAMB) is a patch for the Einstein Boltzmann solver [CAMB](https://github.com/cmbant/CAMB) that intrdouces phenomenological Modifications of Growth (MG) along with dynamical Dark Energy (DE). It includes several phenomenological parametrizations. For instance:

- the mu, gamma parametrization, defined as
<p align="center">
<img src="img/mu_gamma.png" width="350" title="mu gamma parametrization" />
</p>

- the mu, Sigma parametrization, defined as
<p align="center">
<img src="img/mu_sigma.png" width="350" title="mu sigma parametrization" />
</p>

- the Q,R parametrization, defined as
<p align="center">
<img src="img/q_r.png" width="350" title="Q R parametrization" />
</p>

MGCAMB is implemented in the latest version of [CosmoMC](https://github.com/cmbant/CosmoMC) and can be used with [Cobaya](https://cobaya.readthedocs.io/en/latest/). The [MGCosmoMC](https://github.com/sfu-cosmo/MGCosmoMC) and [MGCobaya](https://github.com/sfu-cosmo/MGCobaya) codes are also available for download. 

### Structure of the code
The new MGCAMB patch is structured as in the figure.

<p align="center">
<img src="img/MGCAMB_flowchart.png" width="1000" title="MGCAMB code structure" />
</p>

The parameters in  [``` params_MG.ini ``` ](inifiles/params_MG.ini) are used to run the code and follow the structure above. Please note that dynamical DE is supported in the ``` pure_MG_models ```, where DE perturbations could also be included, and ``` cubic-spline model ```. 



### Citing MGCAMB
If you use MGCAMB for your scientific work, please cite the following papers:

* *New MGCAMB tests of gravity with CosmoMC and Cobaya*\
    Zhuangfei Wang, Seyed Hamidreza Mirpoorian, Levon Pogosian, Alessandra Silvestri, Gong-Bo Zhao\
    [arXiv:2305.05667 [astro-ph.CO]](https://arxiv.org/abs/2305.05667)


* *MGCAMB with massive neutrinos and dynamical dark energy*   
    Alex Zucca, Levon Pogosian, Alessandra Silvestri, and Gong-Bo Zhao  
    [arXiv:1901.05956 [astro-ph.CO]](https://arxiv.org/abs/1901.05956)
    
    
* *Testing Gravity with CAMB and CosmoMC*  
    Alireza Hojjati, Levon Pogosian, Gong-Bo Zhao,  
    [arXiv:1106.4543 [astro-ph.CO]](https://arxiv.org/abs/1106.4543), [JCAP 1108:005,2011](http://iopscience.iop.org/article/10.1088/1475-7516/2011/08/005)
    
    
* *Searching for modified growth patterns with tomographic surveys*  
    Gong-Bo Zhao, Levon Pogosian, Alessandra Silvestri, Joel Zylberberg,  
    [arXiv:0809.3791 [astro-ph]](http://arxiv.org/abs/0809.3791), [Phys. Rev. D 79, 083513](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.79.083513)


as well as the original CAMB [paper](http://arxiv.org/abs/astro-ph/9911177) and [code](https://github.com/cmbant/CAMB).

## 2. How to install
To install MGCAMB on your machine, simply run
```bash
git clone https://github.com/sfu-cosmo/MGCAMB_v4.git
cd MGCAMB_v4/fortran/
make 
```

## 3. How to run
To run MGCAMB, first modify the  [``` params_MG.ini ``` ](inifiles/params_MG.ini) file according to which models you want to use. Then run

```bash
./camb ../inifiles/params.ini
```



## 4. What's new
New features with this new version of the code:
- implemented the QSA models where only dark matter is coupled to scalar field
- added background model based on parameterization of the DE density Omega_x, as opposed to w
- implemented a non-parametric reconstruction of mu, Sigma and Omega_x based on cubic-spline interpolation
- implemented the direct mu, Sigma parametrization for different models
- the option of including DE perturbations to achieve the consistency with CAMB when working with $w \neq −1$ background models. (controlled by parameter `MGDE_pert`)
- added a python wrapper to run MGCAMB from the python interface\
  (make sure to use function ``` set_mgparams ``` alongside other functions e.g. ``` set_cosmology ``` or ``` set_params ``` in your script.)
- the code has been restructured and upgraded to CAMB version 1.3.5 
- added compatibility with Cobaya [(MGCobaya)](https://github.com/sfu-cosmo/MGCobaya)
  
The MG and DE parametrizations along with the computation of the quantities related to the perturbations are introduced in the file [``` mgcamb.f90 ```](fortran/mgcamb.f90).

## 5. Known issues
- The debug mode (`make Debug`) is not well-tested for this version. Please run the code using the normal option described in sections [2](#2-how-to-install) and [3](#3-how-to-install) .
- (**Fixed**) The setting of Dark energy interface and the reading issue of X_arr for `DE_model == 3` in Fortran core. Thanks to Kushal Lodha. 
- The model of effective Newton's constant (``mugamma_par = 3``) is not fully developed.


## 6. Authors List
Main Developers:
- Zhuangfei Wang (Email: zhuangfei_wang@sfu.ca)
- Alex Zucca (Email: alexzucca90@gmail.com)

Original Code Developers:
* [Gong-Bo Zhao](http://english.nao.cas.cn/aboutus/directors/202103/t20210321_265637.html)
* Alireza Hojjati
* [Levon Pogosian](http://www.sfu.ca/%7Elevon/)
* [Alessandra Silvestri](http://wwwhome.lorentz.leidenuniv.nl/%7Esilvestri/Home.html)



Repo created and maintained by Zhuangfei Wang. If you find any bugs in the code, please contact Zhuangfei Wang at zhuangfei_wang@sfu.ca .

<p align="center">
    <a href="http://www.sfu.ca/physics.html"><img src="https://www.sfu.ca/content/sfu/sfuwest/cc-2015/registration/sponsored-students/jcr:content/main_content/image_0.img.640.medium.jpg/1430033284084.jpg" width="380px" height="200px" ></a>
    <a href="http://www.sfu.ca/physics/cosmology/"><img src="https://avatars0.githubusercontent.com/u/7880410?s=280&v=4" width="200px"></a>
</p>

# Programs for quantifying tooth movements for the aligner performance study reported in [1].

## Data definition
**dataDefinitionsCaseStudy.m**

Program includes definition of data and registration parameters.

## Registration and evaluation program

**surfaceMeshRegistrationAllStates.m**

**Input data** are 

(1) weekly oral scans for first nine weeks and for one scan at the end of treatment
```
n=1:11, where n=1  is Week 0 before treatment
              n=10 is Week 9 and 
              n=11 is at the end of treatment
```
(2) semiautomatic segmentation of crowns in the end state at Week 9 (endN=10=refN).

**Processing** includes

(a) manual definition of occlusion plane in end state (endN=10=refN).

(b) rigid and then non-rigid iterative closest point (ICP) registration
rigid:      global alignment of all teeth
non-rigid:  local alignment of individual teeth

(c) evaluation of displacements from state n->endN->simN for the crowns based on non-rigid motion, i.e. non-rigid versus rigid ICP results defined only for moving mesh hence need to deform reference 
and displacements components are then in moving space. The processing is done in two runs:
```
run0: do n=[1:refN-1 refN+1:N]   endN=refN
      propagate semiautomatic teeth segmentation and occlusion plane
      from state endN=refN to all states n=[1:refN-1 refN+1:N] 
      via surface registration
       
run1: do n=1 endN=[2:N]
      determine motion from state n=1 to endN=2:N
      and difference between endN and simN
      via surface registration
```

## Supporting software

**Freecad** was used to visualize and cut the surface meshes.

**OnyxCeph3**<sup>TM</sup> (Image Instruments GmbH, Chemnitz, Germany) was employed to semiautomatically segment the individual crowns at Week 9.

**nricp.m** from github.com/charlienash/nricp was used for non-rigid ICP registration.

**Kabsch.m** from www.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm was employed for fitting a rigid transformation to the non-rigid motion field of a crown.


May 2024, Christine Tanner, Biomaterials Science Center, University of Basel

## Reference
Please cite [1] if you use this software.

[1] I. Filippon, C. Tanner, J. A. von Jackowski, G. Schulz, T. Toepper, B. Mueller.
Aligner-induced tooth movements in three dimensions using clinical data of two patients, 
Oral, 4(4), 487-504, 2024
https://doi.org/10.3390/oral4040039

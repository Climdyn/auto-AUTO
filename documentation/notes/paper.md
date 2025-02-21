---
title: "auto-AUTO: A Python Layer for Automatically Running the AUTO-07p Continuation Software"
tags:
- Bifurcation Analysis
- Continuation
- Dynamical Systems
- Python
authors:
- name: Jonathan Demaeyer
  affiliation: "1"
  orcid: 0000-0002-5098-404X
- name: Oisín Hamilton
  affiliation: "1, 2"
  orcid: 0000-0002-0447-1657
affiliations:
  - name: Meteorological and Climatological Information Service, Royal Meteorological Institute of Belgium, Brussels, Belgium
    index: 1
  - name: Université catholique de Louvain
    index: 2
data: 20 February 2025
bibliography: auto2.bib
---



# Summary

auto-AUTO (or AUTO$^2$) is a python package that acts as an intermediate layer between the user and AUTO-07p continuation 
software [@doedel2007] (AUTO from here on). auto-AUTO automates the running of AUTO by monitoring the continuation, 
while also keeping track of bifurcation points. 
The package can automatically continue along these branching points, and continue branching along further bifurcation points 
to attempt to construct full bifurcation trees. To achieve this in a reliable way, AUTO computations are managed by auto-AUTO 
by setting appropriate stopping conditions, such as meeting other branches, looping branches, and specified bifurcations. 
These stopping conditions supplement the usual AUTO bounds and ensure that branches of detected solutions are unique and well 
defined.

In addition, auto-AUTO provides a comprehensive and documented Python API to investigate properties of the computed 
continuations and solutions. auto-AUTO can be run in Jupyter 
notebooks [@Kluyver2016jupyter] and the results can be plotted with Matplotlib [@Hunter2007] using predefined plotting methods.
This facilitates the integration of the results from AUTO with other Python workflows. 



# Statement of Need

AUTO is a highly optimised and tested continuation code base, and for this reason is one possible choice for use in 
bifurcation analysis studies [when the number of variable is higher than 10.](https://youtu.be/4n8iGysPgus?si=LxoJamGKU7JQ8Be5&t=3601)
However, the use of AUTO in continuation analysis requires a steep learning curve as the existing Python integration 
is limited. 
In addition, there is limited support for automated tracking and continuation of branching points in AUTO.
Therefore, when faced with problems of great complexity (typically when the dimensionality of the problem is greater than 10), 
the user must engage in a long and tedious analysis, restarting the computations many times, with the risk of loosing tracks of 
computation and missing key elements in the process.
This package aims to solve this problem, by automating the running of AUTO, and by systematizing the investigation of branching points.

Tipping points have recently received increased attention, in part due to the concern about possible tipping points being crossed in the climate system, due to climate change [@lenton2008, @wunderling2021, @armstrongmckay2022]. Bifurcation tipping (B-tipping) [@ashwin2012] is one possible mechanism for tipping. Conceptual models, with a relatively low number of degrees of freedom, are still used extensively in increasing our knowledge about tipping points [@wunderling2021]. For this reason, robust but easy to use continuation software is of great importance in this analysis.






# Existing Alternatives

* BifrucationKit.jl [@veltz2020]
* PyDSTool [@clewley2007]
* Other useful tools can be found at [https://dsweb.siam.org/Software](https://dsweb.siam.org/Software)



# Acknowledgments

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under 
the Marie Sklodowska–Curie grant agreement no. 956170. 
In addition, funding has been provided through the "Fédération Wallonie-Bruxelles" with the 
instrument "Fonds Spéciaux de Recherche".



#  References






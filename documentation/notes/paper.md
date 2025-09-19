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
  - name: UCLouvain --- Earth and Life Institute, Belgium
    index: 2
data: 20 February 2025
bibliography: auto2.bib
---



# Summary

auto-AUTO (or AUTO$^2$) is a Python package that acts as an intermediate layer between the user and AUTO-07p continuation 
software [@doedel2007] (AUTO from here on). auto-AUTO automates the running of AUTO by monitoring the continuation, 
while also keeping track of bifurcation points. 
The package can automatically continue along these branching points, and continue branching along further bifurcation points 
to attempt to construct full bifurcation trees. To achieve this in a reliable way, AUTO computations are managed by auto-AUTO 
by setting appropriate stopping conditions, such as meeting other branches, looping branches, and specified bifurcations. 
These stopping conditions supplement the usual AUTO bounds and ensure that branches of detected solutions are unique and well 
defined.

In addition, auto-AUTO provides a comprehensive and documented Python API to investigate properties of the computed 
continuations and solutions. auto-AUTO can be run in Jupyter 
notebooks [@Kluyver2016jupyter] and the results can be plotted with Matplotlib [@Hunter2007] using predefined 
plotting methods (see Figure \ref{fig1}).
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

auto-AUTO is currently used by the authors to investigate the bifurcations and stability of a coupled land-atmosphere model [@hamilton2025] 
using the `qgs` model [@demaeyer2020]. This is part of a wider research project on the concept of 
climate tipping points [@lenton2008; @ashwin2012; @wunderling2021; @armstrongmckay2022], and which typically 
involves bifurcation analysis in high-dimensional dynamical systems[^1].
For this reason, robust but easy to use continuation software is of great importance in this analysis.

[^1]: Although conceptual models, with a relatively low number of degrees of freedom, are still used extensively in increasing our 
knowledge about tipping points [@wunderling2021]. 

![Example of plots using auto-AUTO functionalities to study a high-dimensional model (from the `RP1982` example notebook studying the Reinhold and Pierrehumbert 1982 model [@reinhold1982]). \label{fig1}](bif1.png)

# Existing Alternatives

Over time, and alongside the development of AUTO, many tools were introduced to help the users deal with the complexity 
of the output data stream of AUTO. We can mention for 
example [XPPAUT](https://sites.pitt.edu/~phase/bard/bardware/xpp/xpp.html) [@ermentrout2012], which provides a 
frontend to AUTO, making it easy to continue the solutions to boundary value problems as well as 
equilibria, fixed points, and limit cycles [@ermentrout2007].

Another more recent development is [PyCoBi](https://github.com/pyrates-neuroscience/PyCoBi), which 
"provides a Python interface to Auto-07p, allowing for a more intuitive usage of Auto-07p commands within Python scripts."
[@gast2025]

Other bifurcation analysis software must also be mentioned here:

* BifurcationKit.jl [@veltz2020]
* PyDSTool [@clewley2007]
* Other useful tools can be found at [https://dsweb.siam.org/Software](https://dsweb.siam.org/Software)

However, to the best of the authors' knowledge, these tools do not provide automatic continuation features relevant to the 
problems raised in the [Statement of Need](#statement-of-need). PyCoBi proposes limited automatic generation for 
codimension-2 bifurcation diagrams (regime diagrams), a feature that is also planned in auto-AUTO future developments.

# Future Developments

Future development plans for auto-AUTO include:

* Automatic continuation of regime diagrams, i.e., codimension-1 bifurcation branch tracking and detection of 
  codimension-2 bifurcations by AUTO (regime diagrams)
* Investigation of the formalization and abstraction of automatic continuation using graph theory, leveraging the functionalities already 
  available in auto-AUTO, and allowing for a more insightful construction of the diagrams

However, while the latter proposed development is an obvious way to pursue the development of auto-AUTO, 
it is presently out of the scope of the authors' research.

# Acknowledgments

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under 
the Marie Sklodowska–Curie grant agreement no. 956170. 
In addition, funding has been provided through the "Fédération Wallonie-Bruxelles" with the 
instrument "Fonds Spéciaux de Recherche".
The authors thank the three referees for their constructive review and their suggestions that helped improve the code 
and the manuscript.



#  References






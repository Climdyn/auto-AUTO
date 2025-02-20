
Rationale behind auto-AUTO
==========================

auto-AUTO redefines the notion of `branch` found in |AUTO|. A branch in auto-AUTO regroups both a `forward` and a `backward`
continuations together. A base class :class:`~auto2.continuations.base.Continuation` was designed to include shared
methods and attributes for the different types (fixed points, periodic orbits, ?).
For the moment, there are two main classes of continuations:

* :class:`~auto2.continuations.fixed_points.FixedPointContinuation`: to perform continuation of fixed points based on initial data or on a AUTOSolution
  object.
* :class:`~auto2.continuations.periodic_orbits.PeriodicOrbitContinuation`: to perform continuation of periodic orbits based on initial data or on a AUTOSolution
  object.

These two classes form the building blocks used by the classes devoted to the construction of the tree of bifurcation
and regime diagrams. For the moment, only bifurcation diagrams are possible, with the class:

* :class:`~auto2.diagrams.bifurcations.BifurcationDiagram`

Regime diagrams class development is planned in the future.

These diagrams classes include the logic to expand and compute the bifurcation tree, and will perform these computations
up to pre-specified level.

auto-AUTO has been designed with bifurcations of continuous-time dynamical system in mind. algebraïc and discrete-time
systems could work but have not been tested. Homoclinic continuations are out of the scope of auto-AUTO.

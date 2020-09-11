---
title: 'pfsspy: A Python package for potential field source surface modelling'
tags:
  - Python
  - astronomy
authors:
  - name: David Stansby
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Anthony Yeates
    orcid: 0000-0002-2728-4053
    affiliation: 2
  - name: Samuel Badman
    affiliation: "3, 4"
affiliations:
 - name: Mullard Space Science Laboratory, University College London, Holmbury St. Mary, Surrey RH5 6NT, UK
   index: 1
 - name: Department of Mathematical Sciences, Durham University, Durham, DH1 3LE, UK
   index: 2
 - name: Physics Department, University of California, Berkeley, CA 94720-7300, USA
   index: 3
 - name: Space Sciences Laboratory, University of California, Berkeley, CA 94720-7450, USA
   index: 4
date: 9 September 2020
bibliography: paper.bib
---

# Summary
Magnetic fields play a crucial role in the dynamics and evolution of our Sun
and other stars. A common method used to model the magnetic field of the Sun
and other stars is the potential field source surface (PFSS) model [@Altschuler1969; @Schatten1969].
The PFSS equations assume that there is zero electrical current in the domain of
interest, leading to the equations
\begin{equation}
	\nabla \cdot \mathbf{B} = 0
\end{equation}
\begin{equation}
	\nabla \times \mathbf{B} = 0
\end{equation}
These are solved in a spherical shell between 1 solar radius (1$R_{\odot}$) and
a configurable outer radius called the 'source surface' ($R_{ss}$). Boundary
conditions are given by the radial component of $\mathbf{B}$ on inner boundary
and the imposed condition of a purely radial field on the outer boundary at $R_{ss}$.

Historically, the most widely used software package for performing PFSS extrapolations
within the solar physics community is `pfsspack`. This forms part of the larger `SolarSoft`
library for solar physics, written in the proprietary Interactive Data Language (IDL) programming
language. As more solar physicists use python instead of IDL [@Bobra2020], there is a need to provide similar functionality
within solar physics python ecosystem, a niche that `pfsspy` fills.


# pfsspy
`pfsspy` is a Python package for solving the PFSS equations, and carrying out
other common related tasks such as tracing magnetic field lines through the
solution, importing various magnetic field data sources, and visualising all of this data.

`pfsspy` is designed to closely integrate with other packages in the astronomical and solar physics Python
ecosystem. Coordinate aware input and ouptut maps are created with the SunPy package [@sunpy], and the package is fully integrated with the coordinate and unit frameworks present in astropy [@astropy]. This makes it easy to combine
magnetic fields and field lines calculated in `pfsspy` with other data sources.
As an example, figure \autoref{fig:example} shows magnetic field lines overplotted
on an image of a large active region on the Sun.

![An image of the Sun taken by SDO/AIA at 193 angstroms, with magnetic field lines calculated by pfsspy overplotted in white.\label{fig}](pfsspy.pdf)

The solar physics community has already made use of `pfsspy` in a number of
works, from interpreting observations from Parker Solar Probe [@Bale2019, @Badman2020],
interpreting observations of coronal mass ejections [@Maguire2020], and drawing
links between the Sun and the solar wind [@Stansby2020]. We hope that it provides
a useful resource for the community in the future.


# Acknowledgements

David Stansby acknowledges STFC grants ... and ... for financial support.

# References

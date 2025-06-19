# BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian Inference and Descriptors-Based Machine Learning

![BayeSAF workflow](img/workflow.png)
**_BayeSAF_ framework**: workflow description. We propose an innovative methodology for formulating physicochemical surrogate mixtures that emulate the properties of real jet fuel or drive alternative jet fuel design processes. First, we develop an extensive hydrocarbon property database (992 chemical compounds from 7 hydrocarbon molecular groups) encompassing individual components of candidate sustainable aviation fuels (SAFs) based on available experimental measurements and descriptors-based machine learning (DB-ML) techniques, addressing: (i) 8 lumped properties (critical temperature T<sub>c</sub>, critical pressure P<sub>c</sub>, critical molar volume V<sub>c</sub>, critical density ρ<sub>c</sub>, critical compressibility factor Ζ<sub>c</sub>, acentric factor ω, and derived cetane number _DCN_); 7 temperature-dependent thermophysical properties (liquid-phase density ρ<sub>l</sub>, liquid-phase dynamic viscosity μ<sub>l</sub>, liquid-phase isobaric specific heat C<sub>p,l</sub>, liquid-phase thermal conductivity k<sub>l</sub>, vapor pressure p<sub>v</sub>, surface tension σ, and latent heat of vaporization H<sub>v</sub>). Thereafter, the MATLAB<sup>©</sup>-based _BayeSAF_ algorithm resorts to the hydrocarbon property database and Bayesian inference techniques to formulate surrogate mixtures for alternative jet fuel emulation and design. 

**USER-DEFINED INPUTS:** 

(i) target physicochemical properties (real fuel emulation or jet fuel design)

(ii) number of surrogate mixture components, N<sub>s</sub>
(iii) corresponding hydrocarbon molecular groups, ℋ

**What is _BayeSAF_ inferring then?**

_BayeSAF_ infers **λ** = [**X**; **n<sub>C</sub>**; **η<sub>Β</sub><sup>*</sup>**], i.e., a set of (3N<sub>s</sub>-1) compositional parameters - molar fractions, carbon atom counts, and normalized topochemical atom indices for branching - that univocally identifies mixtures comprising chemical compounds from the hydrocarbon property database. 

[The following is an excerpt from `report.pdf`. View the document to see the full analysis and conclusions.]

# Introduction

The most common way to treat cancer in a patient is through radiation therapy. A radioactive isotope is
injected into the patient and causes the DNA of the cancerous cells to break down, resulting in death of
the cancerous cells. In a paper by Neira, Gago-Arias, Guiu-Souto, and Pardo-Montero (2020), a model is
introduced for the “kinetics of populations of tumor cells” undergoing continuous radiation therapy. The
model consists of three categories of cells: healthy (N), sub-lethally damaged (N<sub>s</sub>), and doomed cells (N<sub>d</sub>).
All cells reproduce with different rates logistically toward a fixed saturation population for the entire tumor.
Healthy cells damaged by radiation can become sub-lethally damaged or lethally damaged (doomed). Sub-
lethally damaged cells can repair themselves and become healthy, whereas doomed cells cannot and eventually
die. Since the radiation applied is continuous, sub-lethally damaged cells can also be further damaged and
become doomed.

This simple model has a potential for exciting applications. However, the paper does not thoroughly
investigate the dynamics of the model. Additionally, two simplifications are offered to the model that appear
to change the dynamics of the model potentially significantly. In this report, we investigate the dynamics of
this continuous dose rate model and its two simplifications. Using parameter values obtained from figures
in the model’s paper we perform simulations of all three models and compare their dynamics. Further,
we attempt to reformulate the ODEs of the model as stochastic differential equations (SDEs) and perform
stochastic simulations of the results.

<p align="center">
  <img src="https://github.com/rbottoms18/tumor-model/blob/main/img/model_flow_chart_opaque.png" width="400"/>
</p>

# Code

`main.m` contains functions for modeling the tumor models and radiaiton doses, as well as plotting the results found in `report.pdf`.

`sde_1.m`, `sde_2.m`, and `sde_3.m` contain three attempts at modeling stochastic differential equations for the model.

`model_fp.nb` contains Mathematica code for computing fixed point and stabilities of fixed points for both the general and simplified models.

`stochastic_matrix.nb` contains Mathematica code for computing the matrix used in the formulation of the SDEs.


# Results

The simplified model that Neira, Gago-Arias, Guiu-Souto, and Pardo-Montero (2020) described mimicked the behavior of the general model almost
exactly. As seen from the fixed point analysis, despite the simplified model compressing fixed curves down to
single points, these fixed points followed the fixed curves. This kept the dynamics of the two models almost
exact, so that no qualitative changes occurred in one but not the other. Additionally, for both models the
computed bifurcation value r of the radiation for the fixed point at the origin remained the same. In light
of this, our conclusion is that the simplified model that the authors used is a good approximation of the
general model.

On the other hand, the LQ limit reduction the authors performed does not accurately reflect the behavior
of either the general or simplified models. The LQ limit version has a single fixed point at the origin that
remains an attractor for all feasible values for a, b, and r. Under this, with only the slightest presence of
radiation that would not deter the general model, the LQ limit will predict that the tumor will be eventually
eradicated. This makes the LQ model an inadequate reduction, and one may wonder what insight is gained
through their comparison to the LQ model (Lea and Catcheside 1942, Kellerer and Rossi 1974, Fowler 1989).

Lastly, we completed our goal of computing an equivalent set of stochastic differential equations for the
general model, but we failed to simulate these equations numerically despite three different attempts. The
variance matrix revealed itself to be much larger than anticipated, rendering it unwieldy. Additionally, the
matrix V in two of our simulation attempts either contained infinite entries or became singular, preventing
successful simulation.

# Acknowledgements

This assignment was completed as a quarter long project for AMATH 536 Mathematical Modeling of Cancer with Prof. Mamis
at the University of Washington, Spring 2024. The writeup of the assignment can be read in `report.pdf`.

Link to the paper discussed: https://iopscience.iop.org/article/10.1088/1361-6560/aba21d/meta

---
math: true
title: "Estimating the Age of the Earth with the Heat Equation"
date: 2026-06-08 # TODO: change date when the time comes!!!!!
categories: [Applied Math]
tags: [diffusion equation, geophysics]     # TAG names should always be lowercase
---

## Introduction 
In this tutorial, we'll see boundary effects and source terms in action as we use the heat equation to estimate the age of the Earth. The presentation here draws largely from {% cite Korner1988 %} and {% cite TS1990 %}. For an excellent physical discussion of heat conduction in the Earth, see chapter 4 of {% cite Lowrie2007 %}, and for historical discussions see {% cite Burchfield1990 %}.

## Kelvin's Estimate
In the nineteenth century, William Thomson (known better as Lord Kelvin)
attempted to approximately compute the age of the Earth using the tools
of thermal physics. As described in {% cite Korner1988 %}, his approach was
roughly based on the assumptions below:

1.  at a reference time $t=0$ early in the life of the solar system, the
    entire Earth is at a uniform reference temperature
    $T_{0}=1200^{\circ}\text{C}$ (roughly, a representative value for the melting temperature of
    rock);

2.  the temperature at the Earth's surface is at a reference value
    $0^{\circ}$C;

3.  the Earth's curvature has a negligible effect on the distribution of
    temperature throughout the crust;

4.  the inside of the Earth is (on average) not moving, so that heat can
    only be transported via conduction;

5.  the Earth's crust is materially homogeneous, so the thermal
    conductivity is constant throughout the surface;

6.  there are no chemical reactions driving the production of heat
    within the Earth.

Note that chapter 4 of {% cite Lowrie2007 %} says that Kelvin's value of $T_{0}$ is larger than that in {% cite Korner1988 %}, but as it turns out this doesn't affect the order of magnitude of our age estimate (!). Under these simplifying assumptions, the Earth can roughly be modelled by the vertical half-space

$$
\left\{ z\ \vert \  0\leq z\leq \infty\right\}.
$$

We treat $z=0$ as the Earth's surface, and as $z\rightarrow +\infty$ we get closer and closer to the planet core (technically, it may be more natural to use the negative half-space $(-\infty,0]$ instead, but this introduces an extra negative sign that's awkward to deal with). Since conduction is assumed to be the only heat transport mechanism, we can determine the temperature of the Earth $T(z,t)$ by solving an initial-boundary value problem for the heat equation:

$$
    \begin{equation}
        \label{eqn:kelvin_problem}
        \left\{
            \begin{aligned}
                T_{t} - \kappa T_{zz} &= 0
                \\
                T\big\vert_{z=0} &= 0
                \\
                T\big\vert_{t=0} &= T_{0}
            \end{aligned}
        \right.
    \end{equation}
$$

where $\kappa$ is a material constant determined from geological observations. According to {% cite TS1990 %}, an accepted value for $\kappa$ is 

$$
\kappa = 6\times 10^{-3} \ \text{cm}^2/\text{s}.
$$

Now, how can we use a solution of \eqref{eqn:kelvin_problem} to estimate the age of the Earth? Again we appeal to geological field observations: an empirical value for the temperature gradient at Earth's surface at the present day is given by 

$$
\gamma = T_{z}\big\vert_{z=0}(t) = \left(3\times 10^{-4} \right)^{\circ}\text{C}/\text{cm}. 
$$

In geophysics, $\gamma$ is called the **geothermal gradient**. See the table below for a summary of all parameter values in the model.  

<div align="center">
  <img width="700" height="300" alt="kelvin table" src="assets/img/kelvin_param_homog.png"/>
</div>

&nbsp; &nbsp; &nbsp; &nbsp; Since $\gamma$ is known, we can solve \eqref{eqn:kelvin_problem} to obtain $T(t)$, compute $T_z\big\vert_{z=0}$, set $\gamma=T_z\big\vert_{z=0}$, and solve for $t$ to determine the age of the Earth. Let's put this plan into action! First, we solve the problem \eqref{eqn:kelvin_problem} using the **method of reflection** (also known as the **method of images**, especially in electromagnetic theory). This involves introducing an odd extension of $T(z,t)$ toall of $z\in \mathbb{R}$ according to

$$
    \widetilde{T}(z,t) \doteq \begin{cases} T(z,t) \quad & z\geq 0
    \\
    T(-z,t) \quad & z<0. 
    \end{cases}
$$

Now, $\widetilde{T}(z,t)$ is odd in $z$ for all $t>0$, so it automatically vanishes at $z=0$. In turn, 

$$
\left(\widetilde{T}\vert_{z\geq0}\right)\bigg\vert_{z=0} = T(z,t)\vert_{z=0} = 0. 
$$

In other words, we can solve \eqref{eqn:kelvin_problem} by solving 

$$
    \begin{equation}
        \label{eqn:kelvin_problem_extended}
        \left\{
            \begin{aligned}
                \widetilde{T}_{t} - \kappa \widetilde{T}_{zz} &= 0
                \\
                \widetilde{T}\big\vert_{t=0}(z) &= \begin{cases} \phantom{-} T_{0} \quad & z\geq 0
                \\
                -T_{0} \quad & z< 0
                \end{cases} 
            \end{aligned}
        \right.
    \end{equation}
$$

and then defining $T(z,t)\doteq \widetilde{T}(z,t)\vert_{z\geq0}$. Using the representation formula for the heat equation (see section 7.2 of {% cite Choksi2022 %}), we have  

$$
    \begin{align*}
        \widetilde{T}(z,t) &= T_{0} \int_{0}^{\infty}S(z-z',t) \ \mathrm{d} z' - T_{0}\int_{-\infty}^{0}S(z-z',t) \ \mathrm{d} z'
        \\
        &= T_{0} \int_{0}^{\infty} S(z-z',t)-S(z+z',t) \ \mathrm{d} z',
    \end{align*}
$$

where we recall that 

$$
S(z,t) = \frac{1}{\sqrt{4\pi \kappa t}} e^{-\frac{z^2}{4\kappa t}}
$$

gives the heat kernel over the entire real line. In other words, we know that the heat kernel for $[0,\infty)$ is 

$$
    \begin{equation}
        \label{eqn:half_line_heat_kernel}
        G(z,z',t)\doteq S(z-z',t)-S(z+z',t) . 
    \end{equation}
$$

Therefore, we have found that the solution of \eqref{eqn:kelvin_problem} is 

$$
    \begin{equation}
        T(z,t) = T_{0}\int_{0}^{\infty} G(z,z',t) \ \mathrm{d} z'. 
    \end{equation}
$$

By doing some simple Gaussian integrals, we find that the above expression can be simplified to 

$$
    \begin{equation}
        \label{eqn:kelvin_soln}
        T(z,t) = T_{0}\mathrm{erf}\left(\frac{z}{\sqrt{4\kappa t}}\right)
    \end{equation}
$$

where the error function is defined by the definite integral

$$
\mathrm{erf}(y) = \frac{2}{\sqrt{\pi}}\int_{0}^{y} e^{-w^2} \ \mathrm{d} w. 
$$

Here is a schematic plot of $T(z,t):

<div align="center">
  <img width="350" height="350" alt="kelvin soln plot" src="assets/img/kelvin_plot.png"/>
</div>

This pictures tells us that the planetary core $z=\infty$ stays at $T_{0}$, but near the surface $z=0$ there is a buffer region where the temperature quickly transitions from $T_{0}$ to $0$. As time elapses, the size of this buffer region increases, corresponding to gradual cooling of the Earth's interior. Notice also that the solution is always smooth for $t>0$, which is a universal feature of heat flow on the real line. 

&nbsp; &nbsp; &nbsp; &nbsp; Since we know the geothermal gradient $\gamma$ at the present day, we can roughly estimate the age of the Earth by computing $T_{z}\vert_{z=0}(t)$, then solving 

$$
T_{z}\vert_{z=0}(t) = \gamma
$$

for $t$. A direct computation yields

$$
    \begin{align*}
        \gamma &= T_{z}\vert_{z=0}
        \\
        &= \frac{T_{0}}{\sqrt{\pi \kappa t}} \left[e^{-z^2/4\kappa t}\right]_{z=0}
        \\
        &= \frac{T_{0}}{\sqrt{\pi \kappa t}}. 
    \end{align*}
$$


Plugging in the values from the table above, we find 

$$
t \approx 30 \ \text{million years}. 
$$

Thus, under Kelvin's assumptions, we have roughly estimated the age of the Earth at around $30$ million years. Even in the nineteenth century, this result was known to be problematic: this estimate is at odds with geological observations, and appears to prohibit sufficient time for Darwinian evolution to develop complex life. So, some physical process is clearly not captured by the assumptions listed at the beginning of this section!

## Modelling Heat Sources within the Earth
In formulating the mathematical model \eqref{eqn:kelvin_problem}, we ignored the possibility that there is some chemical reaction within the Earth generating heat. In the later part of the nineteenth century, however, an alternative mechanism for heat generation was discovered: *radioactivity*. That is, in light of radioactive decay of rock within the Earth, we should have a source of heat in our problem. 

&nbsp; &nbsp; &nbsp; &nbsp; Accordingly, we want to pick a suitable function $f(z,t)$ describing heat creation by radioactivity within the Earth. Then, we solve the *inhomogeneous* problem 

$$
    \begin{equation}
        \label{eqn:real_problem}
        \left\{
            \begin{aligned}
                T_{t} - \kappa T_{zz} &= f(z,t)
                \\
                T\big\vert_{z=0} &= 0
                \\
                T\big\vert_{t=0} &= T_{0}
            \end{aligned},
        \right.
    \end{equation}
$$

compute the geothermal gradient, and solve for $t$ to determine the age of the Earth. 
&nbsp; &nbsp; &nbsp; &nbsp;  Before experimenting with some simple choices of source function $f(z,t)$, we describe how to solve \eqref{eqn:real_problem} for a general $f(z,t)$ using our knowledge of the half-line heat kernel $G(z,z',t)$. Our strategy is to use **Duhamel's principle**. This consists of the following recipe:

1.  Write $w=T-T^{h}$, where $T^{h}$ solves the homogeneous problem
    \eqref{eqn:kelvin_problem}. Then, we know that $w$ solves

    $$
    \left\{
        \begin{aligned}
            w_{t} - \kappa w_{zz} &= f(z,t)
            \\
            w\big\vert_{z=0} &= 0
            \\
            w\big\vert_{t=0} &= 0
        \end{aligned}.
    \right.
    $$

2.  Re-define our (linear!) spatial differential operator by
    $A=\kappa \partial_{z}^2$, whence the above becomes

    $$
    \begin{equation}
        \label{eqn:w_problem}
        \left\{
            \begin{aligned}
                w_{t} - Aw &= f(z,t)
                \\
                w\big\vert_{z=0} &= 0
                \\
                w\big\vert_{t=0} &= 0
            \end{aligned}.
        \right.
    \end{equation}
    $$

3.  Since $A$ is really just a matrix on function space, we may use the
    **method of integrating factors** from ODE theory to formally write
    
    $$
    \partial_{t}\left(e^{-tA}w\right) = -e^{-tA}Aw+e^{-tA}w_{t} = e^{-tA}f.
    $$

    Integrating the above with respect to $t$ tells us that the solution
    to \eqref{eqn:w_problem} is given by 
    $$
    \begin{aligned}
        w(z,t) &= e^{tA}w\vert_{t=0}+\int_{0}^{t}e^{(t-s)A}f(z,s) \ \mathrm{d} s \nonumber 
        \\
        &= \int_{0}^{t}e^{(t-s)A}f(z,s) \ \mathrm{d} s . 
    \end{aligned}
    $$

4.  Remembering the definition of $w$ and that $T^{h}$ is given by \eqref{eqn:kelvin_soln}, we know that 

    $$
        \begin{equation}
            \label{eqn:T_general}
            T(z,t) = T_{0} \ \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa t}}\right) + \int_{0}^{t}e^{(t-s)A}f(z,s) \ \mathrm{d} s.
        \end{equation}
    $$

Therefore, to completely solve \eqref{eqn:real_problem}, it remains to sensibly define $e^{tA}$ for $A=\kappa \partial_{z}^2$. We proceed by analogy with the matrix exponential: if $A$ is an $N\times N$ real matrix and $x(t)$ is an $\mathbb{R}^N$-valued function of time, then the ODE

$$
\begin{equation}
\left\{
    \begin{aligned}
        \frac{\mathrm{d} x}{\mathrm{d} t} - Ax(t) &= 0
        \\
        x(0) &= x_{0}
    \end{aligned}
\right.
\end{equation}
$$

is solved by

$$
    x(t) = e^{tA}x_{0}, \quad e^{tA} = \sum_{k=0}^{\infty} \frac{1}{k!}\left(tA\right)^{k}. 
$$

Now, when $A$ is a linear *differential operator* instead of simply a matrix, we can imagine the PDE 

$$
    v_{t}(z,t)-Av(z,t) =0
$$

as an ODE in "function space": in other words, $t\mapsto v(z,t)$ is a curve in a real vector space of suitably nice functions of $z$.  By analogy, then, it makes sense to define $e^{tA}v_{0}(z)$ for a suitably nice function $v_{0}(z)$ as the solution to the initial-value problem

$$
    \left\{
        \begin{aligned}
            v_{t}(z,t)-Av(z,t) &=0
            \\
            v|_{t=0} &= v_{0}
        \end{aligned}.
    \right.
$$

When we want to impose a homogeneous boundary condition on $e^{tA}v_{0}(z)$ (as in our geochronology problem), we should instead define the exponential as the solution to the initial-boundary-value problem

$$
    \left\{
        \begin{aligned}
            v_{t}(z,t)-Av(z,t) &=0
            \\
            v|_{z=0} &= 0 
            \\
            v|_{t=0} &= v_{0}
            \\
        \end{aligned}.
    \right.
$$

In the case $A=\kappa \partial_{z}^2$, then, we know by earlier discussion that 

$$
    \begin{equation}
        e^{tA}v_{0}(z) = v(z,t) = \int_{0}^{\infty} G(z,z',t) \ v_{0}(z') \ \mathrm{d} z',
    \end{equation}
$$

where $G(z,z',t)$ is defined by \eqref{eqn:half_line_heat_kernel}. We conclude that \eqref{eqn:T_general} can be written more precisely as 

$$
\begin{equation}
    \label{eqn:T_specific}
    T(z,t) = T_{0} \ \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa t}}\right) +\int_{0}^{t} \int_{0}^{\infty} G(z,z',t-s) \ f(z', s) \ \mathrm{d} z' \ \mathrm{d} s. 
\end{equation}
$$


{% bibliography --cited %}
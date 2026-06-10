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

In geophysics, $\gamma$ is called the **geothermal gradient**. See Table 1 below for a summary of all parameter values in the model.  

| Parameter  | Name | Value | Units |  
| :-----: | :-----: | :-----: | :-----: |
| $T_0$ | Initial Temperature of Earth | $1200$ | $\phantom{0}^{\circ}\mathrm{C}$ |
| $\kappa$  | Thermal Conductivity of Earth | $6\times 10^{-3}$ | $\text{cm}^2/\text{s}$ |
| $\gamma$ | Geothermal Gradient | $3\times 10^{-4}$ | $\phantom{0}^{\circ}\text{C}/\text{cm}$ |

**Table 1:** Parameters for the homogeneous problem. 

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


Plugging in the values from table 1, we find 

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


> Our definition of the exponential of a differential operator as the solution to an ODE in a space other than $\mathbb{R}^n$ generalizes nicely to other contexts. For instance, the exponential on a general Lie group is also defined by solving an ODE, with an element of the corresponding Lie algebra serving as initial data; see for example chapter 7 of \cite{Lee2013}. 
{: .prompt-info }

> Viewing a PDE as an ODE on a function space is the starting point for the **semigroup approach** to solving initial-value problems. In this context, we say that the family of integral operators $$ \left\{e^{\kappa t\partial_{z}^2} \ \vert \ t\geq 0\right\} $$ 
form the **semigroup generated by $\kappa \partial_{z}^2$**. This set formally has the structure of an abelian semigroup (abelian group without inverses) upon noticing that
$$
e^{\kappa t_{1}\partial_{z}^2}e^{\kappa t_{2}\partial_{z}^2} = e^{\kappa t_{2}\partial_{z}^2}e^{\kappa t_{1}\partial_{z}^2} = e^{\kappa (t_{1}+t_{2})\partial_{z}^2}. 
$$
We cannot upgrade this to a group structure, since the heat equation is not well-posed for $t\in (-\infty,0]$. For a good introduction to semigroups, see chapter 6 of {% cite MH1994 %} (be prepared: semigroup theory requires a lot of functional analysis!). 
{: .prompt-info }

&nbsp; &nbsp; &nbsp; &nbsp; Now that we have a representation formula for a solution of \eqref{eqn:real_problem} given by \eqref{eqn:T_specific}, we can start playing around with particular choices of source function $f(z,t)$. For simplicity, today we only consider time-independent sources $f=f(z)$. Let's begin with the simplest choice 

$$
    f(z) \equiv B\kappa
$$

with $B$ a constant: according to \cite{TS1990}, a suitable empirical value for $B$ is $B=\left(3.25 \times 10^{-10}\right)^{\circ}\text{C}/\text{cm}^2$. Such a choice of $f(z)$ corresponds to a uniform distribution of radioactive rock through the whole interior of the Earth. Using our earlier computations of $\int_{0}^{\infty} G(z,z',t) \ \mathrm{d} z',$ we find

$$
    \begin{align*}
        T(z,t) = T_{0} \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa t}}\right) + B\kappa\int_{0}^{t} \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa (t-s)}}\right) \ \mathrm{d} s .
    \end{align*}
$$

Then, we form the equation for the geothermal gradient: 
$$
    \begin{align*}
        \gamma &= T_{z}|_{z=0}
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{B\kappa}{\sqrt{\pi\kappa}}\int_{0}^{t}\frac{\mathrm{diff} s}{\sqrt{t-s}} 
        \\
        &= \frac{1}{\sqrt{\pi\kappa}}\left(T_{0}t^{-1/2}+B\kappa t^{1/2}\int_{0}^{1}\tau^{-1/2}\ \mathrm{diff} \tau\right)
        \\
        &= \frac{1}{\sqrt{\pi\kappa}}\left(T_{0}t^{-1/2}+2B\kappa t^{1/2}\right).
    \end{align*}
$$

Accordingly, $t$ must satisfy the quadratic equation

$$
0 = 4T_{0}\kappa^2B^2t^2+\left(4B\kappa T_{0}-a^2\right)t+1.
$$

However, this equation clearly has no real roots! Accordingly, $f\equiv \text{constant}$ cannot provide a realistic model of radioactive heat sources inside the Earth. We have arrived at the scientific conclusion that *the concentration of radioactive rock within our planet varies with distance from the core*. For any alternative discussion of why $f\equiv \text{constant}$ is not a good model, see {% cite TS1990 %}. 

&nbsp; &nbsp; &nbsp; &nbsp; Since we know the source function $f$ cannot be constant, let us take the simplest non-constant $f(z)$ possible: 

$$
    f(z) = \begin{cases} 
    B\kappa \quad & 0\leq  z\leq H
    \\
    0 \quad & z>H.
    \end{cases}
$$

This corresponds to all radioactive material being concentrated in a layer of thickness $H$ below the surface of the Earth. Again, we can take $B=\left(3.25 \times 10^{-10}\right)^{\circ}\text{C}/\text{cm}^2$. For the layer thickness $H$, we assume $H$ is on the order of the mean thickness of the Earth's crust (about $1/640$ of the Earth's radius): a suitable concrete value is $10$ km; see table 2 below for a listing of all parameter values:

| Parameter  | Name | Value | Units |  
| :-----: | :-----: | :-----: | :-----: |
| $T_0$ | Initial Temperature of Earth | $1200$ | $\phantom{0}^{\circ}\mathrm{C}$ |
| $\kappa$  | Thermal Conductivity of Earth | $6\times 10^{-3}$ | $\text{cm}^2/\text{s}$ |
| $\gamma$ | Geothermal Gradient | $3\times 10^{-4}$ | $\phantom{0}^{\circ}\text{C}/\text{cm}$ |
| $B$ | Radioactive Heat Source Strength | $3.25 \times 10^{-10}$ | $\phantom{0}^{\circ}\text{C}/\text{cm}^2$|
| $H$ | Thickness of Radioactive Layer | $10$ | $\text{km}$ |

**Table 2:** Parameters for the inhomogeneous models. 

Using \eqref{eqn:T_specific}, we can write the temperature distribution corresponding to this $f$ as

$$
    \begin{equation}
        T(z,t) = T_{0}\mathrm{erf}\left(\frac{z}{\sqrt{4\kappa t}}\right) + B\kappa\int_{0}^{t}\int_{0}^{H} G(z,z',t-s) \ \mathrm{diff} z' \ \mathrm{diff} s. 
    \end{equation}
$$

The equation for the geothermal gradient then reads

$$
    \begin{align*}
        \gamma &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + B\kappa\int_{0}^{t}\int_{0}^{H} \partial_{z}G(z,z',t)|_{z=0} \ \mathrm{d} z' \ \mathrm{d} s
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{B\kappa}{2\sqrt{\pi}}\int_{0}^{t}\int_{0}^{H} \left(\kappa (t-s)\right)^{-3/2} z' e^{-(z')^2/4\kappa (t-s)} \ \mathrm{d} z' \ \mathrm{d} s
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} - \frac{B\kappa}{2\sqrt{\pi}}\int_{0}^{t}\int_{0}^{H} \left(\kappa (t-s)\right)^{-3/2} \ \left[2\kappa (t-s)\right] \ \frac{\partial}{\partial z'}\left[e^{-(z')^2/4\kappa (t-s)} \right] \ \mathrm{d} z' \ \mathrm{d} s
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} - \frac{B\kappa}{\sqrt{\pi}}\int_{0}^{t} \left(\kappa (t-s)\right)^{-1/2} \left[e^{-(z')^2/4\kappa t} \right]_{z=0}^{H} \ \mathrm{d} s
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{B\kappa}{\sqrt{\pi}}\int_{0}^{t} \left(\kappa (t-s)\right)^{-1/2} \left[1-e^{-H^2/4\kappa (t-s)} \right] \ \mathrm{d} s
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} - \frac{B\kappa}{\sqrt{\pi}}\int_{t}^{0} \left(\kappa \tau\right)^{-1/2} \left[1-e^{-H^2/4\kappa\tau} \right] \ \mathrm{d} \tau
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} +\frac{B\kappa}{\sqrt{\pi}}\int_{0}^{t} \left(\kappa \tau\right)^{-1/2} \left[1-e^{-H^2/4\kappa\tau} \right] \ \mathrm{d} \tau
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{2B\sqrt{\kappa}}{\sqrt{\pi}}t^{1/2} -\frac{B\sqrt{\kappa}}{\sqrt{\pi}} \int_{0}^{t} \tau^{-1/2} \ e^{-H^2/4\kappa \tau} \ \mathrm{d} \tau.
    \end{align*}
$$

If we change variables according to

$$
\sigma^{2} = \frac{H^2}{4\kappa\tau}
$$

then

$$
2\sigma \  \mathrm{d} \sigma = -\frac{H^2}{4\kappa\tau^2} \mathrm{d} s = \frac{\sigma^{2}}{\tau} \ \mathrm{d} \tau
$$

hence 

$$
\mathrm{d} \tau = \frac{2\tau}{\sigma} \ \mathrm{d} \sigma. 
$$

We then find 

$$
    \begin{align*}
        \gamma &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{2B\sqrt{\kappa}}{\sqrt{\pi}}t^{1/2} - \frac{BH}{\sqrt{\pi}}\int_{\frac{H}{2\sqrt{\kappa t}}}^{\infty} \frac{e^{-\sigma^2}}{\sigma^2}  \ \mathrm{d} \sigma.
    \end{align*}
$$

Next, an integration by parts shows

$$
    \int_{\zeta}^{\infty} \frac{e^{-\sigma^2}}{\sigma^2} \ \mathrm{d} \sigma = \frac{e^{-\zeta^2}}{\zeta^2} -2\int_{\zeta}^{\infty} e^{-\sigma^2} \ \mathrm{d} \sigma =  \frac{e^{-\zeta^2}}{\zeta
    }+\sqrt{\pi}\left(\mathrm{erf}\left(\zeta\right) -1\right).
$$

Using this formula, we conclude

$$
    \begin{align}
        \gamma &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{2B\sqrt{\kappa}}{\sqrt{\pi}}t^{1/2} - \frac{BH}{\sqrt{\pi}}\left[\frac{e^{-\zeta^2}}{\zeta}+\sqrt{\pi}\left(\mathrm{erf}\left(\zeta\right) -1\right)\right]_{\zeta=\frac{H}{2\sqrt{\kappa t}}}
        \nonumber
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} + \frac{2B\sqrt{\kappa}}{\sqrt{\pi}}t^{1/2} - \frac{BH}{\sqrt{\pi}}\left[\frac{e^{-H^2/4\kappa t}}{\frac{H}{2\sqrt{\kappa t}}}+\sqrt{\pi}\left(\mathrm{erf}\left(\frac{H}{2\sqrt{\kappa t}}\right)-1\right)\right]
        \nonumber
        \\
        &= \frac{T_{0}}{\sqrt{\pi\kappa t}} +\frac{2B\sqrt{\kappa}}{\sqrt{\pi}}t^{1/2}\left(1-e^{-H^2/4\kappa t}\right) -BH\left(\mathrm{erf}\left(\frac{H}{2\sqrt{\kappa t}}\right)-1\right).
        \label{eqn:gamma_expression_full}
    \end{align}
$$

This expression can be vastly simplified using asymptotic approximations valid as $t\rightarrow\infty$. First, define 

$$
    \begin{equation}
        \tau = t^{-\frac12}.
    \end{equation}
$$

As $t\rightarrow\infty$, $\tau\rightarrow 0^{+}$. In terms of $\tau$, \eqref{eqn:gamma_expression_full} reads

$$
    \begin{equation}
        \label{eqn:radioactive_layer_geo_grad_tau}
        \gamma = \frac{T_{0}}{\sqrt{\pi\kappa }}\tau +\frac{2B\sqrt{\kappa}}{\sqrt{\pi}}\tau^{-1}\left(1-e^{-\tau^2H^2/4\kappa }\right) -BH\left(\mathrm{erf}\left(\frac{H}{2\sqrt{\kappa}}\tau\right)-1\right).
    \end{equation}
$$

When $\tau \ll 1$, we may use Taylor expansion to discover

$$
    \begin{align}
        1-e^{-\tau^2H^2/4\kappa } &= -\frac{H^2}{4\kappa} \tau^2 + \mathcal{O}\left(\tau^4\right)
        \\
        \mathrm{erf}\left(\frac{H}{2\sqrt{\kappa}}\tau\right)-1 &= -1 + \frac{H}{\sqrt{\kappa\pi}}\tau + \mathcal{O}\left(\tau^{3}\right)
    \end{align}
$$

Plugging these expansions into \eqref{eqn:radioactive_layer_geo_grad_tau} and simplifying, we obtain 

$$
    \begin{equation*}
        \gamma = BH + \frac{T_0 -\frac32 BH^2}{\sqrt{\pi\kappa}} \tau + \mathcal{O}\left(\tau^3\right), \quad \tau \rightarrow 0^{+}. 
    \end{equation*}
$$

In terms of $t$, this becomes 

$$
    \begin{equation}
        \gamma = BH + \frac{T_0 -\frac32 BH^2}{\sqrt{\pi\kappa}} t^{-\frac12} + \mathcal{O}\left(t^{-\frac32}\right), \quad t \rightarrow +\infty. 
    \end{equation}
$$

Discarding the $\mathcal{O}\left(t^{-\frac32}\right)$ terms and solving for $t$ gives 

$$
    \begin{equation}
        \label{eqn:t_final_radioactive}
        t \approx \left[\frac{T_0 -\frac32 BH^2}{\sqrt{\pi\kappa}\left(\gamma-BH\right)}\right]^{\frac12}.
    \end{equation}
$$

If we insert parameter values from table 2, the above formula yields 

$$
t\approx 6.3 \times 10^8 \ \text{years} = \mathcal{O}\left(10^9\right) \ \text{years}.
$$

This at least gives a suitable magnitude for the age of the Earth! Accordingly, assuming all radioactive rock is concentrated in a small layer near the Earth's surface gives a solid conceptual model for the actual age of our planet. 


>  Whenever your units are a bit funny (as they are in this problem), I recommend using [pint](https://pint.readthedocs.io/en/stable/) (or [MetPy](https://unidata.github.io/MetPy/latest/index.html)) to perform computations: this is how I evaluated the right-hand side of \eqref{eqn:t_final_radioactive}. These software packages can really save you from wasting time making silly unit conversion mistakes! 
{: .prompt-info }

{% bibliography --cited %}
---
math: true
title: "Introduction to Equatorial Waves in the Stratosphere"
date: 2026-05-15
categories: [Applied Math, Waves, Atmospheric Dynamics]
tags: [atmospheric dynamics, middle atmosphere, fluid mechanics, waves]     # TAG names should always be lowercase
---

## Introduction 
In this post, I provide a brief introduction to linear waves propagating in the equatorial stratosphere. The discussion here follows the classical monograph of Andrews, Holton, and Leovy {% cite AHL1987 %} very closely. For an alternative reference, see chapter 17 of the **second edition** of Vallis' textbook {% cite Vallis2017 %}. For coverage of equatorial waves in the troposphere see Majda {% cite Majda2003 %}, chapter 8 of Vallis {% cite Vallis2017 %}, or the original paper of Matsuno {% cite Matsuno1966 %}. 

## Primitive Equations in Log-Pressure Coordinates
Our first task is to re-formulate the primitive equations in a suitable
coordinate system. Let's set up some notation: explanations of any
notation not defined here may be found in any standard atmospheric
dynamics text such as {% cite Holton2004 %} or {% cite Vallis2017 %}. If we denote the
material derivative by 

$$
    \begin{equation}
    \label{eqn:Dt}
    D_t = \partial_t + u\partial_x + v\partial_y + \omega \partial_p = \partial_t + u\partial_x + v\partial_y + w \partial_z.
    \end{equation}
$$

and the *horizontal* material derivative by 

$$
    \begin{equation}
    \label{eqn:Dth}
    D_t^h = \partial_t + u\partial_x + v\partial_y,
    \end{equation}
$$

then in isobaric coordinates the primitive equations take the form (cf. {% cite Holton2004 %} Ch. 3)

$$
    \begin{alignat}{3}
        \label{eqn:isobaric_momentum}
        D_t \mathbf{v} + \mathbf{f}(y) \times \mathbf{v} &= -\nabla_p \Phi \quad &&\text{(horizontal momentum)},
        \\
        \partial_p \Phi &= - \frac{RT}{p} \quad &&\text{(hydrostatic balance)},
        \\
        \partial_x u + \partial_y v + \partial_p \omega &= 0 \quad  &&\text{(continuity)},
        \\
        \label{eqn:isobaric_thermo}
        D_t^h T - S_p \omega &= c_p^{-1}J \quad &&\text{(thermodynamic)}. 
    \end{alignat}
$$

This $5\times 5$ system of nonlinear PDE is certainly familiar to every student of atmospheric science. In the context of middle atmosphere applications, however, alternatives to isobaric coordinates are a common sight. For example, one can justifiably use isentropic coordinates owing to the strong static stability characterizing the stratosphere. In these notes, we use **log-pressure** coordinates since they result in a simplification of \eqref{eqn:isobaric_thermo} where the static stability parameter $S_p$ no longer appears (see {% cite Holton2004 %}, section 8.4.1). 

&nbsp; &nbsp; &nbsp; &nbsp; The vertical log-pressure coordinate $z$ is defined with respect to a
choice of **scale height** $H>0$ and **surface pressure** $p_s>0$ as

$$
    \begin{equation}
    \label{eqn:log_pressure_coord}
    z = - H \ln\left(p/p_s\right).
    \end{equation}
$$ 

The scale height gives rise to a choice of surface temperature $T_S$ according to

$$
    \begin{equation}
    \label{eqn:scale_height}
    H = \frac{RT_s}{g};
    \end{equation}
$$

the motivation for this choice comes from solving the equation of hydrostatic balance in an isothermal atmosphere where $T\equiv T_S$.

&nbsp; &nbsp; &nbsp; &nbsp; Next, we collect some preliminary results and definitions required to convert \eqref{eqn:isobaric_momentum}-\eqref{eqn:isobaric_thermo} to log-pressure coordinates. First and most importantly, the multivariable chain rule tells us how to change vertical coordinates.

$$
    \begin{equation}
    \label{eqn:vertical_coord_derivative_change}
    \partial_p = \frac{\partial z}{\partial p}\partial_z = -\frac{H}{p} \partial_z.
    \end{equation}
$$

Also, since the primitive equations involve $D_t$, we must determine the vertical velocity in log-pressure coordinates. Using \eqref{eqn:vertical_coord_derivative_change}, 

$$
    \begin{equation}
    \label{eqn:logp_vertical_velocity}
    w = \frac{Dz}{Dt} = -\frac{H}{p}\frac{Dp}{Dt} = -\frac{H}{p}\omega.
    \end{equation}
$$

Let's define a **reference density profile** by

$$
    \begin{equation}
    \label{eqn:reference_density}
    \rho_0(z) = \rho_s \exp\left(-z/H\right)
    \end{equation}
$$ 
where $p_s = R \rho_s T_s$. A quick exercise with the hydrostatic balance equation shows that this would be the *actual* density profile if $T\equiv T_S$. Additionally, we define the **buoyancy frequency** by

$$
    \begin{equation}
    \label{eqn:logp_buoyancy_freq}
    \widetilde{N}^2 = \frac{R}{H} \left[\frac{\kappa T}{H} + \partial_z T\right].
    \end{equation}
$$

Note that this is *not* the Brunt-Väisälä buoyancy frequency $N^2$
commonly encountered in atmospheric thermodynamics. However, the two
(angular!) frequencies are related by 
$$
    \begin{equation}
    \label{eqn:bvf_to_bf}
    \widetilde{N}^2 = \left(\frac{T}{T_s}\right)^2 N^2.
    \end{equation}
$$ 

One often encounters literature where $\widetilde{N}^2$ is denoted by $N^2$, and you should be careful about this notation difference to avoid headaches when verifying calculations. Note that, for an isothermal atmosphere, $\widetilde{N}^2 = N^2$: by now, you should realize that log-pressure coordinates work well when deviations from isothermality are small.

&nbsp; &nbsp; &nbsp; &nbsp; With all of these preliminaries taken care of, we can write the primitive equations in log-pressure coordinates as 

$$
    \begin{alignat}{3}
        \label{eqn:log_pressure_momentum}
        D_t \mathbf{v} + \mathbf{f}(y)\times \mathbf{v} &= -\nabla_z \Phi \quad &&\text{(horizontal momentum)},
        \\
        \label{eqn:log_pressure_hydrostatic_balance}
        \partial_z\Phi &= \frac{RT}{H}  \quad &&\text{(hydrostatic balance)},
        \\
        \label{eqn:log_pressure_continuity}
        \partial_x u + \partial_y v + \rho_0^{-1}\partial_z \left(\rho_0 w\right) &= 0 \quad &&\text{(continuity)}, 
        \\
        \label{eqn:log_pressure_thermo}
        D_t^h T + \frac{\widetilde{N}^2 H}{R} w &= \frac{\kappa}{R}J \quad &&\text{(thermodynamic)}. 
    \end{alignat}
$$

How can we check these formulas? In light of \eqref{eqn:vertical_coord_derivative_change} and \eqref{eqn:logp_vertical_velocity}, the only nontrivial part of the conversion is transforming \eqref{eqn:isobaric_thermo} to \eqref{eqn:log_pressure_thermo}. To see why \eqref{eqn:log_pressure_thermo} is correct, first note that $\kappa = R/c_p$ by definition. Next, we expand $S_p \omega$ to discover

$$
    \begin{alignat*}{3}
        -S_p \omega &= S_p\frac{P}{H} w \quad &&\text{by \eqref{eqn:logp_vertical_velocity}}
        \\
        &= \left[\frac{RT}{c_p p} - \partial_p T\right]\frac{P}{H}w \quad && \text{by definition of } S_p
        \\
        &= \left[\frac{RT}{c_p p} +\frac{H}{P} \partial_z T\right]\frac{P}{H}w \quad && \text{by \eqref{eqn:vertical_coord_derivative_change}} 
        \\
        &= \left[\frac{\kappa T}{H} +\partial_z T\right]w \quad &&\text{by definition of } \kappa
        \\
        &= \frac{\widetilde{N}^2 H}{R} w \quad &&\text{by definition of } \widetilde{N}^2. 
    \end{alignat*}
$$

Owing to the simple relationship \eqref{eqn:bvf_to_bf}, \eqref{eqn:log_pressure_thermo} is much easier to understand than its isobaric coordinate counterpart \eqref{eqn:isobaric_thermo}. 
\par For the wave analysis presented below, we elect to eliminate $T, w$ from the above system to get a reduced $3\times 3$ system for $u,v, \Phi$. If we use \eqref{eqn:log_pressure_thermo} to write $w$ in terms of $T$ and then use \eqref{eqn:log_pressure_hydrostatic_balance} to write $T$ in terms of $\Phi$, the primitive equations become 

$$
    \begin{alignat}{2}
        \label{eqn:log_pressure_momentum_again}
        D_t \mathbf{v} + \mathbf{f}(y)\times \mathbf{v} &= -\nabla_z \Phi,
        \\
        \label{eqn:log_pressure_evol_eqn_for_phi}
        \partial_x u + \partial_y v + \rho_0^{-1}\partial_z \left[\frac{\rho_0}{\widetilde{N}^2}\left(\frac{\kappa}{H}J - D_t^h\partial_z\Phi\right) \right] &= 0. 
    \end{alignat}
$$

## Wave Solutions

Now, we put the primitive equations to work to predict the existence of
a zoo of atmospheric waves. For simplicity, we make the following
assumptions:

1.  the motion approximately occurs in an equatorial $\beta$-plane so
    $$\label{eqn:beta-plane-approx}
                \mathbf{f}(y) = \beta y \ \widehat{z};$$

2.  diabatic heating can be ignored so $J=0$;

3.  the buoyancy frequency $\widetilde{N}^2$ is approximately constant;

4.  the solution to \eqref{eqn:log_pressure_momentum_again}-\eqref{eqn:log_pressure_evol_eqn_for_phi} is a small perturbation of the trivial basic state
    
    $$
    \overline{u} = \overline{v} = \overline{\Phi} = 0.
    $$

Linearizing \eqref{eqn:log_pressure_momentum_again}-\eqref{eqn:log_pressure_evol_eqn_for_phi} under these assumptions, the perturbations $u, v, \Phi$ (we have avoided the customary primes for cleanliness) obey 

$$
    \begin{align*}
        \partial_t u - \beta y v+ \partial_x \Phi &= 0, 
        \\
        \partial_t v + \beta y u + \partial_y\Phi &= 0,
        \\
         \partial_x u + \partial_y v - \rho_0^{-1}\partial_z \left[\frac{\rho_0}{\widetilde{N}^2}\partial_{z}\partial_{t}\Phi \right] &= 0.
    \end{align*}
$$

Using the definition of the reference density $\rho_0$ \eqref{eqn:reference_density}, the above can be written more cleanly as 

$$
    \begin{align}
        \partial_t u - \beta y v+ \partial_x \Phi &= 0, \label{eqn:lin_x_mom}
        \\
        \partial_t v + \beta y u + \partial_y\Phi &= 0, \label{eqn:lin_y_mom}
        \\
         \partial_x u + \partial_y v - \frac{1}{\widetilde{N}^2}\left(\partial_{z}^2\partial_t -\frac{1}{H} \partial_z\partial_{t}\right)\Phi &= 0.  \label{eqn:lin_evol_eqn_phi}
    \end{align}
$$

While \eqref{eqn:lin_x_mom} - \eqref{eqn:lin_evol_eqn_phi} is *linear*, the $\beta$-effect gives rise to a few variable coefficients. 

> In {% cite Majda2003 %} chapter 9, Majda shows how equatorial waves confined to the troposphere can be modelled using a system of linear *hyperbolic* conservation laws. This allows for an elegant change-of-variables under which the governing equations naturally include the **ladder operators** familiar from the theory of the quantum harmonic oscillator (see the discussion below). For our present work involving waves propagating in the stratosphere, however, we find that\eqref{eqn:lin_x_mom} - \eqref{eqn:lin_evol_eqn_phi} is not hyperbolic, so Majda's trick cannot be applied. I think it would be fun to obtain an analogue of \eqref{eqn:lin_x_mom} - \eqref{eqn:lin_evol_eqn_phi} in a non-hydrostatic atmosphere and see if there is enough "hyperbolic structure" to make the ladder operators pop out naturally! 
{: .prompt-info }

&nbsp; &nbsp; &nbsp; &nbsp; To solve \eqref{eqn:lin_x_mom} - \eqref{eqn:lin_evol_eqn_phi} while accommodating the $y$-dependent coefficients, we introduce the ansatz below following {% cite AHL1987 %}, section 4.7: 

$$
    \begin{equation}
        \label{eqn:equatorial_wave_ansatz}
        \left[u, v, \Phi\right]^{\text{T}} = e^{r z} \Re \left\{e^{i\left(kx+mz-\omega t\right)}\left[\widehat{u}(y), \widehat{v}(y), \widehat{\Phi}(y)\right]^{\text{T}}\right\},
    \end{equation}
$$

where T denotes transposition and $r>0$ is an altitudinal growth rate. \eqref{eqn:equatorial_wave_ansatz} represents a wave propagating in both the zonal and vertical directions and growing with increasing altitude; this growth is physically motivated by higher altitudes having a lower density (that is, less inertia). The wave's amplitude also depends on the meridional coordinate $y$. Plugging \eqref{eqn:equatorial_wave_ansatz} into \eqref{eqn:lin_x_mom} - \eqref{eqn:lin_evol_eqn_phi} and being careful about $e^{rz}$, we obtain an ODE system where $y$ is the manipulated variable:

$$
    \begin{align*}
        -i\omega \widehat{u} -\beta y \widehat{v} + ik \widehat{\Phi} &= 0,
        \\
        -i\omega \widehat{v} + \beta y \widehat{u} +\partial_y \widehat{\Phi} &= 0,
        \\
        ik\widehat{u} + \partial_{y}\widehat{v} + \frac{i\omega}{\widetilde{N}^2}\left(r + im\right)\left(\frac{1}{H} - r-im\right) \widehat{\Phi} &= 0.
    \end{align*}
$$

If we choose 

$$
    r = \frac{1}{2H},
$$

then 

$$
    \left(r + im\right)\left(\frac{1}{H} - r-im\right) = \left\vert\frac{1}{2H} + im\right\vert^2 
$$

and our ODE system simplifies to

$$
    \begin{align}
        -i\omega \widehat{u} -\beta y \widehat{v} + ik \widehat{\Phi} &= 0, \label{eqn:wave_odes_1st_one}
        \\
        -i\omega \widehat{v} + \beta y \widehat{u} +\partial_y \widehat{\Phi} &= 0,
        \\
        ik\widehat{u} + \partial_{y}\widehat{v} - \frac{i\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \widehat{\Phi} &= 0. \label{eqn:wave_odes_3rd_one}
    \end{align}
$$

If we did not choose $r=\left(2H\right)^{-1}$, we would end up having to split our system (and, ultimately, our dispersion relation) into real and imaginary parts. This requires a fair amount of non-instructive arithmetic: I am happy to choose the canonical growth rate and press on with the physics. 

&nbsp; &nbsp; &nbsp; &nbsp; At this point, Andrews and friends {% cite AHL1987 %} make the additional assumption 

$$
    \left\vert m\right\vert \gg (2H)^{-1}.
$$

That is, they restrict to the case where many vertical wavelengths can be observed between the lower boundary and two scale heights. Said differently, $\left\vert m\right\vert$ must be much larger than the wave's altitudinal amplitude growth rate. This means \eqref{eqn:wave_odes_3rd_one} can be replaced with 

$$
    ik\widehat{u} + \partial_{y}\widehat{v} - \frac{i\omega m^2}{\widetilde{N}^2}\widehat{\Phi}=0.
$$

To add some value, in these notes we'll keep the $2H$ terms around most of the time; according to {% cite AHL1987 %}, this is useful for upper-stratospheric settings (and therefore the semiannual oscillation of the tropical stratopause).

&nbsp; &nbsp; &nbsp; &nbsp; Next, we reduce \eqref{eqn:wave_odes_1st_one}-\eqref{eqn:wave_odes_3rd_one} to a $2\times 2$ system by eliminating $\widehat{u}$. This gives

$$
    \begin{align}
        -i\omega \widehat{v} + \frac{\beta y}{i\omega}\left(-\beta y \widehat{v} + ik\widehat{\Phi}\right) +\partial_y \widehat{\Phi} &= 0,
        \label{eqn:first_reduced_wave_ode}
        \\
        \frac{k}{\omega}\left(-\beta y \widehat{v} + ik\widehat{\Phi}\right) + \partial_{y}\widehat{v} - \frac{i\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \widehat{\Phi} &= 0.  \label{eqn:second_reduced_wave_ode}
    \end{align}
$$

Clearly, $\widehat{v} = \widehat{\Phi} =0$ solves the above system. Are there any interesting solutions when only one of the variables is identically zero? 

### Case 1: $\widehat{\Phi} =0$

In this case, \eqref{eqn:wave_odes_1st_one}-\eqref{eqn:wave_odes_3rd_one} becomes 

$$
    \begin{align*}
        \left[-\omega + \frac{\beta^2 y^2}{\omega} \right]\widehat{v} &= 0,
        \\
        \frac{k}{\omega}\left(-\beta y \widehat{v}\right) + \partial_{y}\widehat{v} &= 0. 
    \end{align*}
$$

The first equation cannot be satisfied by any constant $\omega$ unless $\widehat{v} \equiv 0$. So, there are no nontrivial solutions in this case. 

### Case 2: $\widehat{v} =0$

Here, \eqref{eqn:wave_odes_1st_one}-\eqref{eqn:wave_odes_3rd_one} simplifies to 

$$
    \begin{align*}
        \frac{\beta y}{\omega}\left(k\widehat{\Phi}\right) +\partial_y \widehat{\Phi} &= 0,
        \\
        \left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \widehat{\Phi} &= 0. 
    \end{align*}
$$

The second of these equations gives the dispersion relation 

$$
    \begin{equation}
        \label{eqn:kelvin_dispersion_relation}
        \omega^2 = \frac{k^2\widetilde{N}^2}{m^2 + \left(2H\right)^{-2}}.
    \end{equation}
$$

The first equation is a separable ODE, which is easily solved to yield

$$
    \begin{equation}
        \label{eqn:kelvin_structure_function}
        \widehat{\Phi}(y) = \widehat{\Phi}_0 \exp\left[-\frac{\beta k y^2}{2\omega}\right]
    \end{equation}
$$

where $\widehat{\Phi}_0 \in \mathbb{C}$ is an integration constant.

### Case 2: $\widehat{\Phi} \neq0$, $\widehat{v} \neq 0$

Isolating $\partial_y \widehat{v}$ in \eqref{eqn:second_reduced_wave_ode}, we have 

$$
    \begin{equation}
        \label{eqn:funny_polarization}
        -\partial_y \widehat{v}  = i\left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \widehat{\Phi} - \frac{k\beta y}{\omega} \widehat{v}.
    \end{equation}
$$

Differentiating both sides with respect to $y$ and using \eqref{eqn:first_reduced_wave_ode}, we discover 

$$
    \begin{align*}
        -\partial_y^2 \widehat{v}  &= i\left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \partial_y\widehat{\Phi} - \partial_y\left\{\frac{k\beta y}{\omega} \widehat{v}\right\}
        \\
        &= i\left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \left\{i\omega \widehat{v} - \frac{\beta y}{i\omega}\left[-\beta y \widehat{v} + ik\widehat{\Phi}\right]\right\}- \partial_y\left\{\frac{k\beta y}{\omega} \widehat{v}\right\}
        \\
        &= i\left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \left\{i\left[\omega - \frac{\beta^2y^2}{\omega} \right] \widehat{v} - \frac{\beta k y}{\omega}\widehat{\Phi}\right\}- \partial_y\left\{\frac{k\beta y}{\omega} \widehat{v}\right\}.
    \end{align*}
$$

Using \eqref{eqn:second_reduced_wave_ode} again to express $\widehat{\Phi}$ in terms of $\widehat{v}, \partial_y \widehat{v}$, we have 

$$
    \begin{align*}
        -\partial_y^2 \widehat{v}  &= -\left\{\frac{k^2}{\omega} - \frac{\omega}{\widetilde{N}^2}\left[m^2 + \left(2H\right)^{-2}\right] \right\} \left\{\omega - \frac{\beta^2y^2}{\omega}\right\}\widehat{v}- \frac{\beta ky}{\omega}\left[\frac{k\beta y}{\omega}\widehat{v} -\partial_y \widehat{v}\right]- \partial_y\left\{\frac{\beta k y}{\omega} \widehat{v}\right\}
        \\
        &= -\left\{k^2 - \frac{\omega^2}{\widetilde{N}^2}\left[m^2 +(2H)^{-2}\right] + \frac{\beta k}{\omega}\right\}\widehat{v} - \frac{\beta^2y^2}{\widetilde{N}^2}\left[m^2+(2H)^{-2}\right] \widehat{v}.
    \end{align*}
$$

We conclude that $\widehat{v}$ obeys the linear ODE 

$$
    \begin{equation}
        \label{eqn:dimensionful_vhat_ode_yucky}
        0 = \left\{\partial_y^2 +\left[\frac{\omega^2}{\widetilde{N}^2}\left[m^2+(2H)^{-2}\right] - k^2 -\frac{k\beta}{\omega}\right] - \frac{\beta^2}{\widetilde{N}^2}\left[m^2+(2H)^{-2}\right] y^2\right\}\widehat{v}.
    \end{equation}
$$

To lighten our notational burden, let's define 

$$
    \begin{equation}
        \langle m\rangle \doteq \sqrt{m^2 + (2H)^{-2}}.
    \end{equation} 
$$

Think of this as a "smoother" version of $\vert m\vert$: we in fact have 

$$
    \langle m\rangle\approx \vert m \vert \quad \text{if} |m|\gg 2H.
$$ 

With this in mind, \eqref{eqn:dimensionful_vhat_ode_yucky} becomes 

$$
    \begin{equation}
        \label{eqn:dimensionful_vhat_ode}
        0 = \left\{\partial_y^2 +\left[\frac{\langle m\rangle^2\omega^2}{\widetilde{N}^2} - k^2 -\frac{k\beta}{\omega}\right] - \frac{\langle m\rangle^2\beta^2}{\widetilde{N}^2} y^2\right\}\widehat{v}.
    \end{equation}
$$

Let's tidy this equation by introducing a linear change of variables 

$$
    y\mapsto \eta = \alpha y
$$

where $\alpha$ is a constant with units of $\text{length}^{-1}$ to be chosen later. The chain rule gives 

$$
    \partial_y = \alpha \partial_\eta
$$

so \eqref{eqn:dimensionful_vhat_ode} becomes 

$$
    0 = \left\{\partial_\eta^2 + \alpha^{-2}\left[\frac{\langle m\rangle^2\omega^2}{\widetilde{N}^2} - k^2 -\frac{k\beta}{\omega}\right] - \alpha^{-4}\left(\frac{\beta\langle m\rangle}{\widetilde{N}}\right)^2\eta^2\right\}\widehat{v}(\eta).
$$

To simplify the term involving $\eta^2$, we should choose 

$$
    \alpha = \left(\frac{\langle m\rangle\beta}{\widetilde{N}}\right)^{\frac12},
$$

which indeed has the correct units. This means 

$$
    \begin{equation}
        \label{eqn:eta_defn}
        \eta =  \left(\frac{\langle m\rangle\beta}{\widetilde{N}}\right)^{\frac12} y.
    \end{equation}
$$

We also compress our notation using 

$$
    \begin{equation}
        \label{eqn:lambda_defn}
        2\lambda + 1 \doteq \frac{\widetilde{N}}{\beta \langle m\rangle}\left\{\frac{\langle m\rangle^2\omega^2}{\widetilde{N}^2} -k^2 -\frac{k\beta}{\omega}\right\},
    \end{equation}
$$

whence \eqref{eqn:dimensionful_vhat_ode} becomes simply 

$$
    \begin{equation}
        \label{eqn:qho_eqn}
        \left\{\partial_\eta^2 + 2\lambda+1 - \eta^2\right\}\widehat{v}(\eta).
    \end{equation}
$$

&nbsp; &nbsp; &nbsp; &nbsp; We investigate \eqref{eqn:qho_eqn} by following the treatment of quantum-mechanical simple harmonic motion in Griffiths' excellent book {% cite Griffiths2018 %}. Let's get some insight into this equation by rewriting it as 

$$
    \partial_\eta^2 \widehat{v} = -\left[2\lambda+1 -\eta^2\right]\widehat{v}.
$$

Therefore, $\widehat{v}(\eta)$ will be 

1.  exponentially decaying if $$ \vert\eta\vert^2 > 2\lambda + 1 $$, or

2.  oscillatory if $$ \vert\eta\vert^2 < 2\lambda + 1 $$. 

Strictly speaking, the first case also permits $\widehat{v}$ to grow exponentially in $\eta$, but to ensure a finite-energy (equatorially-localized) solution we ignore this possibility. Therefore, solutions of \eqref{eqn:qho_eqn} should look roughly like trigonometric functions near $\eta = 0$, but gradually they deform to ensure rapid decay as $\vert\eta\vert\rightarrow\infty$. Further, for $\vert\eta\vert \gg 2\lambda + 1$ we approximately have 

$$
\partial_\eta^2\widehat{v} \approx \eta^2 \widehat{v} \Rightarrow \widehat{v} \approx \text{constant} \times e^{-\eta^2/2} 
$$

via a rough dominant-balance argument. This motivates us to look for solutions of \eqref{eqn:dimensionful_vhat_ode} in the form 

$$
\widehat{v}(\eta) = e^{-\eta^2/2}\psi(\eta).
$$

Grinding out the calculations, this means that $\psi(\eta)$ obeys the **Hermite equation**

$$
    \begin{equation}
        \label{eqn:hermite}
        0 = \psi''(\eta) - 2\eta \psi'(\eta) + 2\lambda \psi(\eta). 
    \end{equation}
$$

The method of power series may be used to show that the only solutions to \eqref{eqn:hermite} with sub-exponential growth (read: corresponding to equatorially-localized $\widehat{v}(\eta)$) require 

$$
    \lambda = n \in \left\{0,1,2,...\right\}.
$$

Further, these solutions are scalar multiples of the **Hermite polynomials** defined by the recurrence relation 

$$
    \begin{align}
        H_0(\eta) &= 1,
        \\
        H_1(\eta) &= 2\eta,
        \\
        H_{n+1}(\eta) &= 2\eta H_n(\eta) - 2n H_{n-1}(\eta). 
    \end{align}
$$

Thus we have an infinite family of solutions to \eqref{eqn:qho_eqn} in the form 

$$
    \begin{equation}
        \label{eqn:qho_solns}
        \widehat{v}(\eta) = \text{constant} \times D_n(\eta), \quad D_n(\eta)\doteq  e^{-\eta^2/2}H_n(\eta) \quad n=0,1,2,... \ .
    \end{equation}
$$

The functions $D_n(\eta)$ in \eqref{eqn:qho_solns} are sometimes called **parabolic cylinder functions**, up to a choice of normalization. Here is a sketch of the first few $D_n$:

<div align="center">
  <img width="625" height="625" alt="pcf_plot" src="assets/img/pcf.png" />
</div>

This plot substantiates our decay/oscillation transition intuition from earlier! Additionally, the illustration shows that that $D_n$ goes through $(n+1)/2$ "wavelengths" before it begins to decay. Note also that the constraint of integer $\lambda$ implies our dispersion relation is parameterized: \eqref{eqn:lambda_defn} implies 

$$
    \begin{equation}
        \label{eqn:eq_wave_disperion_relation}
        \frac{\omega^2\langle m\rangle^2}{\widetilde{N}^2}-k^2 -\frac{k\beta}{\omega} = (2n+1)\frac{\beta \langle m\rangle}{\widetilde{N}}, \quad n= 0,1,2,... \quad . 
    \end{equation}
$$

## Wave Spectrum 

Here, we analyze the wave solutions computed in the previous section in greater detail, paying attention to the various branches of the dispersion relation. 

### Kelvin 

These waves are characterized by having a vanishing meridional velocity disturbance $\widehat{v}=0$. Therefore, their dispersion relation is given by \eqref{eqn:kelvin_dispersion_relation}, which is very similar to that of midlatitude gravity waves {% cite Vallis 2017 %}, section 17.2.2. Since we are primarily interested in upward-propagating waves, we typically focus only on the branch 

$$
    \begin{equation}
        \label{eqn:kelvin_upward_branch}
        \omega = -\frac{k\widetilde{N}}{\langle m\rangle} = -\frac{k\widetilde{N}}{m\sqrt{1 +\vert 2Hm \vert^{-2}}}
    \end{equation}
$$

since this corresponds to an upward group velocity in the limit $\vert 2H m\vert 
\rightarrow + \infty$. The geopotential for Kelvin waves obeys 

$$
    \Phi(x,y,z,t) = \exp\left[\frac{z}{H}-\frac{\beta ky^2}{2\omega}\right] \Re\left\{\widehat{\Phi_0} \ e^{i(kx+mz-\omega t)}\right\}
$$

where $k,m,\omega$ are constrained by \eqref{eqn:kelvin_upward_branch}. Define the zonal phase speed by 

$$
    \begin{equation}
        \label{eqn:kelvin_zonal_phase_speed}
        c = \omega/k.
    \end{equation}
$$

If $c<0$ then the Kelvin wave geopotential $\Phi$ is unphysical and not equatorially localized. If $c=0$, $\Phi$ is not even defined! Therefore, for Kelvin waves, we must have $c>0$, corresponding to \emph{eastward} zonal phase propagation. In light of \eqref{eqn:kelvin_upward_branch}, this means $m<0$ is necessary for upward propagation. Finally, notice how rotational effects only appear in the meridional structure function: $\beta$ and $c$ control the extent of equatorial localization. 

### Yanai/Mixed

These are the simplest waves with a nonzero meridional velocity perturbation. They correspond to the $n=0$ branch of the dispersion relation \eqref{eqn:eq_wave_disperion_relation}. By recognizing a difference of squares, we can write the dispersion relation along this branch as 

$$
    \begin{equation}
        \label{eqn:yanai_simplified_dispersion}
        \omega\left(\langle m\rangle\omega + \widetilde{N}k\right)\left(\langle m\rangle\omega -\widetilde{N}k\right) = \beta\widetilde{N}\left(\langle m\rangle\omega + \widetilde{N}k\right)
    \end{equation}
$$

To simplify this expression further, we return to \eqref{eqn:funny_polarization}. If 

$$
    \langle m\rangle\omega + \widetilde{N}k = 0
$$

so that \eqref{eqn:yanai_simplified_dispersion} is trivial, then in particular the zonal phase speed $c$ obeys $c<0$. However, \eqref{eqn:funny_polarization} then reduces to 

$$
    \partial_y \widehat{v} = \frac{\beta y}{c} \widehat{v} \Rightarrow \widehat{v} = \exp\left[\frac{\beta y^2}{2c}\right],
$$

which is unphysical. We conclude $\langle m\rangle\omega + \widetilde{N}k = 0$, whence \eqref{eqn:yanai_simplified_dispersion} becomes 

$$
    \begin{equation}
        \label{eqn:yanai_dispersion_final}
        \langle m\rangle = \frac{\widetilde{N}}{c^2}\left(c + \frac{\beta}{k^2}\right).
    \end{equation}
$$

Since $m\neq 0$ is required for vertical propagation, Yanai waves obey 

$$
    \begin{equation}
        \label{eqn:yanai_are_super_rossby}
        c > -\frac{\beta}{k^2}.
    \end{equation}
$$

Accordingly, these waves have a zonal phase speed strictly greater than that of pure Rossby waves. This is why they are not just equatorially-modified Rossby waves! 

### Poincaré/Rossby

Our discussion of these modes is very brief: see \cite{AHL1987} or \cite{Vallis2017} for additional details. These auxiliary waves appear when $n>0$ in \eqref{eqn:eq_wave_disperion_relation}. Notice that this dispersion relation is \emph{quadratic} in $\langle m\rangle$, hence the quadratic formula gives

$$
    \begin{equation}
        \langle m \rangle = 
        \begin{cases}
        &\widetilde{N}\beta \omega^{-2} \left\{n+\frac12 + \left[\left(n+\frac12\right)^2 + \frac{\omega k}{\beta}\left(1+\frac{\omega k}{\beta}\right) \right]^{\frac12}\right\}
        \quad \text{(Poincaré branches)},
        \\
        &\widetilde{N}\beta \omega^{-2} \left\{n+\frac12 - \left[\left(n+\frac12\right)^2 +\frac{\omega k}{\beta}\left(1+\frac{\omega k}{\beta}\right)\right]^{\frac12}\right\}
        \quad \text{(Rossby branches)}.
        \end{cases}
    \end{equation}
$$

In the case where $\vert 2Hm\vert \gg 1$, this becomes

$$
    \begin{equation}
        m = 
        \begin{cases}
        -\mathrm{sgn}(\omega)\hspace{-0.35cm}&\widetilde{N}\beta \omega^{-2} \left\{n+\frac12 + \left[\left(n+\frac12\right)^2 + \frac{\omega k}{\beta}\left(1+\frac{\omega k}{\beta}\right) \right]^{\frac12}\right\}
        \quad \text{(Poincaré branches)},
        \\
        &\widetilde{N}\beta \omega^{-2} \left\{n+\frac12 - \left[\left(n+\frac12\right)^2 +\frac{\omega k}{\beta}\left(1+\frac{\omega k}{\beta}\right)\right]^{\frac12}\right\}
        \quad \text{(Rossby branches)}.
        \end{cases}
    \end{equation}
$$

after choosing signs to guarantee a positive vertical component of the group velocity. 

### Nondimensional Dispersion Relations in the Case $\vert 2Hm\vert \gg 1$

Suppose the zonal wavenumber $k$ is fixed and $\vert 2Hm\vert$ is very large. We may then define a dimensionless $\omega$ by 

$$
    \begin{equation}
        \label{eqn:omegahat}
        \widehat{\omega} = \frac{\omega k}{\beta}
    \end{equation}
$$

and a dimensionless, $H$-modified $m$ by 

$$
    \begin{equation}
        \label{eqn:mhat}
        \widehat{m} = \frac{m\beta}{\widetilde{N}k^2}.
    \end{equation}
$$

With these simplifications in mind, the various branches of the dispersion relation reduce to

$$
    \begin{align}
        \widehat{m}_{\text{Kelvin}} &= -\frac{1}{\widehat{\omega}},
        \\
        \left\vert \widehat{m}_{\text{Yanai}}\right\vert &= \frac{\widehat{\omega} + 1}{\widehat{\omega}^2},
        \\
        \widehat{m}_{\text{Poincaré}} &= -\frac{\mathrm{sgn}\widehat{\omega}}{\widehat{\omega}^2} \left\{n+\frac12 + \left[
            \left(n+\frac12\right)^2 + \widehat{\omega}\left(1+\widehat{\omega}\right)
        \right]^{\frac12}
        \right\},
        \\
        \widehat{m}_{\text{Rossby}} &= \frac{1}{\widehat{\omega}^2} \left\{n+\frac12 - \left[
            \left(n+\frac12\right)^2 + \widehat{\omega}\left(1+\widehat{\omega}\right)
        \right]^{\frac12}
        \right\}.
    \end{align}
$$

The figure below, adapted from figure 4.23 in {% cite AHL1987 %}, summarizes these dispersion relations visually.

<div align="center">
  <img width="625" height="625" alt="eq_waves_dispersion_relation" src="assets/img/eq_waves_dispersion.png" />
</div>

For an analogue of this plot valid in the (tropospherically relevant) case of *no vertical propagation*, see figure 8.6 in {% cite Vallis2017 %}. 

### Propagation Conditions

To end this basic coverage of equatorial waves in the stratosphere, we introduce conditions under which Kelvin and Yanai waves can propagate in a nontrivial background flow. This is a baby step towards understanding the Holton-Lindzen-Plumb theory of the quasi-biennial oscillation or QBO \cite{AHL1987, Baldwin2001}: I hope to return to the subject of QBO modelling in a future post. For simplicity, let's assume $\vert 2Hm\vert\gg 1$ for this last little part of our discussion. 

&nbsp; &nbsp; &nbsp; &nbsp; Now, all the above developments trivially generalize to a situation with a constant zonal background flow 

$$
    \mathbf{v} = \overline{u} \ \mathbf{\widehat{x}}
$$

since, by boosting to a frame of reference moving zonally with speed $\overline{u}$, the frequency transforms as

$$
    \omega \mapsto \omega -\overline{u}k.
$$

Accordingly, for upward-propagating Kelvin waves, \eqref{eqn:kelvin_upward_branch} implies that the zonal phase speed $c_{\text{Kelvin}}$ obeys

$$
    c_{\text{Kelvin}} -\overline{u} = \frac{\omega}{k} - \overline{u} = -\frac{\widetilde{N}}{m}.
$$

However, for Kelvin waves, we already know that upward propagation requires $m<0$ (this is a consequence of their *eastward* zonal phase speed). Consequently, a Kelvin wave in a nontrivial zonal background flow must obey

$$
    \begin{equation}
        \overline{u} < c_{\text{Kelvin}}.
    \end{equation}
$$

Said differently, Kelvin waves are inhibited within a strong eastward (positive) mean flow. Further, in the "critical limit" $c_{\text{Kelvin}}\rightarrow \overline{u}^{+}$, the dispersion relation requires $m\rightarrow -\infty$. Therefore, the vertical component of the wave's group velocity tends to zero in this limit.  

&nbsp; &nbsp; &nbsp; &nbsp;  We may play a similar game with westward-moving Yanai waves. In this case, a bit of elbow grease and \eqref{eqn:yanai_dispersion_final} tell us that the vertical component of the group velocity is 

$$
    c_g^{z} = -\frac{\left(\omega-\overline{u}k\right)^2}{\widetilde{N}\left[2\beta + \left(\omega-\overline{u}k\right)k\right]},
$$

However, from dispersion plot above we find that $c_{g}^z>0$ is required, so the above tells us

$$
    -\left(c_{\text{Yanai}} - \overline{u}\right) > 0.
$$

Re-arranging gives 

$$
    \begin{equation}
        \overline{u} > c_{\text{Yanai}},
    \end{equation}
$$

so a strong westward flow inhibits Yanai wave propagation. This is the opposite condition to that encountered for Kelvin waves! Since Kelvin waves always travel eastward and we have assumed a westward Yanai wave in this discussion, we have found that "like opposes like" when it comes to determining how a background flow opposes wave propagation. Finally, notice that, similar to the Kelvin case, $c_g^z$ vanishes for Yanai waves as $c_{\text{Yanai}}\rightarrow \overline{u}^{-}$. Thus, group velocity vanishes as we approach a critical phase speed, irrespective of whether our waves are Kelvin or Yanai. 

{% bibliography --cited %}
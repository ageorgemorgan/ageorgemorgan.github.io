---
math: true
title: "Global Existence for a Nonlinear Wave Equation: A Cool Example from Klainerman and Nirenberg"
date: 2026-05-07
categories: [Math, Analysis of PDEs, Evolution Equations]
tags: [evolution equations]     # TAG names should always be lowercase
---

## Introduction 
Today, we'll study a very nasty-looking nonlinear wave equation that actually admits global-in-time solutions for a suitable class of nice initial data. Specifically, we look at the IVP

$$
    \begin{equation}
        \label{eqn:nlw_kn}
        \left\{
            \begin{aligned}
            v_{tt}-\Delta v  +\left(v_{t}\right)^2-|\nabla_{x} v|^2 &= 0 \quad \forall \ (x,t)\in \mathbb{R}^3\times \mathbb{R},
            \\
            v|_{t=0} &= f(x)\quad \forall \ x\in \mathbb{R}^3, 
            \\
            v_{t}|_{t=0} &= g(x) \quad \forall \ x\in \mathbb{R}^3. 
            \end{aligned} 
            \right.
    \end{equation}
$$

Following {% cite Klainerman1980 %} we'll show that, if $f,g$ are compactly supported and "small enough" in a suitable norm, then a $v(x,t)\in C^{2}\left(\mathbb{R}^3\times \mathbb{R}\right)$ solving \eqref{eqn:nlw_kn} exists. This requires a clever transformation taking our our nonlinear PDE to a *linear* wave equation that may be analyzed using the Kirchhoff formula {% cite Choksi2022 %} (p.130). I originally learned about this problem from the exercises in chapter 6 of {% cite Craig2018 %}.

## Function Spaces
Before beginning in earnest, we need to define norms to quantify what we mean by "small enough" initial data. 

### Definition 
- We denote the set of continuous, bounded, real-valued functions on $$\mathbb{R}^3$$ by $$C_{b}\left(\mathbb{R}^3\right)$$. This becomes a normed real vector space when we define $$\left\|f\right\|_{C_{b}\left(\mathbb{R}^3\right)} \doteq \sup_{x\in \mathbb{R}^3}\left|f(x)\right|$$.

- If $$\mathbf{V}\colon \mathbb{R}^3 \rightarrow\mathbb{R}^3$$ is a vector field with entries in $$C_{b}\left(\mathbb{R}^3\right)$$, we define $$\left\|\mathbf{V}\right\|_{C_{b}\left(\mathbb{R}^3\right)} \doteq \sup_{x\in \mathbb{R}^3}\left|\mathbf{V}(x)\right|$$ where $$\left|\mathbf{V}(x)\right|$$ is the Euclidean norm of $$\mathbf{V}(x)$$. 

- We define $$C^{1}_{b}\left(\mathbb{R}^3\right)$$ to be the subspace of $$C_{b}\left(\mathbb{R}^3\right)$$ consisting of those bounded continuous functions $$f(x)$$ whose gradients $$\nabla f(x)$$ have entries in $$C_{b}\left(\mathbb{R}^3\right)$$. We put another norm on this space by defining $$\left\|f\right\|_{C^{1}_{b}\left(\mathbb{R}^3\right)} \doteq  \left\|f\right\|_{C_{b}\left(\mathbb{R}^3\right)}+ \left\|\nabla f\right\|_{C_{b}\left(\mathbb{R}^3\right)}$$. 

## A Decay Estimate for the Linear Wave Equation 
We'll get comfortable with these new norms by proving the following crucial estimate for solutions of the linear wave equation in $\mathbb{R}^3$; this estimate is one of the two main steps in constructing global-in-time solutions of \eqref{eqn:nlw_kn}. 

### Proposition [Decay Rate for Wave Equation in $\mathbb{R}^3$]

Suppose $$F(x) \in C^{1}_{b}\left(\mathbb{R}^3\right)$$ and $$G(x) \in C_{b}\left(\mathbb{R}^3\right)$$ are supported in the ball of radius $$R>1$$ centred at the origin (denoted $$B(0, R)$$). Consider the initial value problem 
$$ 
    \begin{equation} \label{eqn:wave_ivp}
        \left\{
            \begin{aligned}
            u_{tt}-\Delta u &= 0 \quad \forall \ (x,t)\in \mathbb{R}^3\times \mathbb{R},
            \\
            u|_{t=0} &= F(x)\quad \forall \ x\in \mathbb{R}^3, 
            \\
            u_{t}|_{t=0} &= G(x) \quad \forall \ x\in \mathbb{R}^3. 
            \end{aligned} 
            \right.
    \end{equation} 
$$ 
This problem admits a unique global-in-time solution $$u(x,t) \in C^{2}\left(\mathbb{R}^3\times \mathbb{R}\right)$$ satisfying the following decay property: for all $$t\neq 0$$, 
$$
    \begin{equation}
        \label{eqn:wave_eqn_decay}
        \left\|u(\cdot,t)\right\|_{C_{b}\left(\mathbb{R}^3\right)} \leq \frac{R^2}{|t|}\left(\left\|F\right\|_{C^{1}_{b}\left(\mathbb{R}^3\right)} + \left\|G\right\|_{C_{b}\left(\mathbb{R}^3\right)}\right).
    \end{equation} 
$$ 
In particular, $$u(x,t)\rightarrow 0$$ as $$|t|\rightarrow \infty$$, uniformly in $$x$$.

### Proof
To prove this result, first note that existence of $u(x,t)$ follows from Kirchhoff's formula: 
$$
    \begin{equation}
        \label{eqn:kirchhoff}
        u(x,t) = \frac{1}{4\pi t^2}\oint_{\partial B(x,|t|)} F(y) + \nabla F(y)\cdot(y-x) + |t|G(y) \ \mathrm{d} S(y). 
    \end{equation}
$$
Uniqueness of the Kirchhoff solution then follows from a simple energy argument. So, the only hard part here is to obtain the temporal decay estimate \eqref{eqn:wave_eqn_decay}. We split into two cases depending on $|t|$.

#### <b>Case 1: $|t|\geq 1$.</b> 

By hypothesis and the Kirchhoff formula, $u(x,t)$ depends only on the values of $F$ and $G$ on    $\partial B(x,|t|) \cap B(0,R)$: $F,G$ vanish at all other points in $\partial B(x,|t|)$. Now, notice that the largest sphere that can be inscribed inside the closure $\overline{B(0,R)}$ is simply $\partial B(0,R)$. Therefore, 
$$ 
    \text{Area}\left(\partial B(x,|t|) \cap B(0,R)\right)\leq  \text{Area}\left(\partial B(0,R)\right)=4\pi R^2. 
$$
With this observation in mind, we have  
$$ 
    \begin{align*}
        |u(x,t)| &\leq \frac{1}{4\pi t^2} \oint_{\partial B(x,|t|) \cap B(0,R)} |F(y)| + |\nabla F(y)| \ |y-x| + |t| \ |G(y)| \ \mathrm{d} S(y)
        \\
        &\leq \frac{1}{4\pi t^2} \ \left(\left\|F\right\|_{C_{b}}+|t| \ \left\|\nabla F\right\|_{C_{b}} + |t|\ \left\|G\right\|_{C_{b}}\right)\oint_{\partial B(x,|t|) \cap B(0,R)}\mathrm{d} S(y)
        \\
        &\leq \frac{R^2}{t^2}\left(\left\|F\right\|_{C_{b}}+|t|\ \left\|\nabla F\right\|_{C_{b}} + |t|\ \left\|G\right\|_{C_{b}}\right)
        \\
        &\leq \frac{R^2}{|t|}\left(\left\|F\right\|_{C_{b}}+\left\|\nabla F\right\|_{C_{b}} + \left\|G\right\|_{C_{b}}\right)
        \\
        &=\frac{R^2}{|t|}\left(\left\|F\right\|_{C^{1}_{b}}+ \left\|G\right\|_{C_{b}}\right), 
    \end{align*} 
$$ 

so this case is all done. 

#### <b>Case 2: $|t|< 1$.</b> 

As in the previous case, we apply the Kirchhoff formula to estimate 
$$ 
    \begin{align}
        |u(x,t)| &\leq \frac{1}{4\pi t^2} \oint_{\partial B(x,|t|)} |F(y)| + |\nabla F(y)| \ |y-x| + |t| \ |G(y)| \ \mathrm{d} S(y) \nonumber
        \\
        &\leq \left\|F\right\|_{C_{b}}\frac{1}{4\pi |t|^2} \oint_{\partial B(x,|t|)} \mathrm{d} S(y) + \left(\left\|\nabla F\right\|_{C_{b}}+\left\| G\right\|_{C_{b}}\right) \frac{1}{4\pi |t|} \oint_{\partial B(x,|t|)} \mathrm{d} S(y) \nonumber
        \\
        &\leq \left\|F\right\|_{C_{b}} + \left(\left\|\nabla F\right\|_{C_{b}}+\left\| G\right\|_{C_{b}}\right) |t| \nonumber 
        \\
        &< \left\|F\right\|_{C^{1}_{b}}+ \left\|G\right\|_{C_{b}}, \label{eqn:important_intermediate}
    \end{align}
$$ 
where we have used $|t|<1$ to obtain the last line. We conclude by observing that, for $|t|<1$ and $R>1$, we automatically have $1< R^2/|t|$. 

Putting both these cases together, the proof is finished. $\square$ 

## Back to the Nonlinear Problem 

At long last, we're ready to show that the big bad IVP \eqref{eqn:nlw_kn} always has at least one solution for small enough ICs. 

### Theorem [Klainerman-Nirenberg {% cite Klainerman1980 %}, pp. 45-46]

Suppose $$f(x) \in C^{1}_{b}\left(\mathbb{R}^3\right)$$ and $$g(x) \in C_{b}\left(\mathbb{R}^3\right)$$ are supported in the ball of radius $$R>1$$ centred at the origin. Then, \eqref{eqn:nlw_kn} admits a global-in-time solution $$v(x,t) \in C^{2}\left(\mathbb{R}^3\times \mathbb{R}\right)$$ provided that $$\left\|f\right\|_{C^{1}_{b}\left(\mathbb{R}^3\right)}$$  and $$\left\|g\right\|_{C_{b}\left(\mathbb{R}^3\right)}$$ are sufficiently small. 

### Proof 
The proof is delightfully simple. First, we make a change of responding variable $v\mapsto u$ according to
$$
    \begin{equation}
    \label{eqn:v_to_u}
        u(x,t) = e^{v(x,t)} -1. 
    \end{equation}
$$
The inverse transformation is 
$$
    \begin{equation}
    \label{eqn:u_to_v}
        v(x,t) = \log\left(u(x,t)+1\right). 
    \end{equation}
$$
Note that this inverse is only well-defined if $u>-1$. 

&nbsp; &nbsp; &nbsp; &nbsp; Let's pretend we have found a $v(x,t)$ solving \eqref{eqn:nlw_kn}. What IVP $u(x,t)$ does satisfy? A direct computation shows that 
$$
    \begin{align*}
        u_{t} &= v_{t}e^{v}
        \\
        u_{x_{i}} &= v_{x_{i}}e^{v}, \quad i=1,2,3
    \end{align*}
$$
so 
$$
    \begin{align*}
        u_{tt} &= \left[v_{tt}+(v_{t})^2\right]e^{v}
        \\
        u_{x_{i}} &=  \left[v_{x_{i}x_{i}}+(v_{x_{i}})^2\right]e^{v}, \quad i=1,2,3. 
    \end{align*}
$$
Using the nonlinear wave equation satisfied by $v$ then gives 
$$
    \begin{align*}
        0 &= \left[u_{tt}e^{-v}-(v_{t})^2\right] -\left((\Delta u) e^{-v}- |\nabla u|^2\right)  +\left(v_{t}\right)^2-|\nabla_{x} v|^2
        \\
        &= e^{-v}\left[u_{tt}-\Delta u\right].
    \end{align*}
$$
Therefore, $u$ satisfies the linear wave equation 
$$
    \begin{equation}
        \label{eqn:kn_u_is_lw}
        u_{tt}-\Delta u =0 \quad \forall \ (x,t) \in \mathbb{R}^3\times \mathbb{R}. 
    \end{equation}
$$
Now, we check what initial conditions $u$ satisfies. By inspection of the ICs in \eqref{eqn:nlw_kn}, we see that 
$$
    \begin{align}
        u|_{t=0} &= F(x) \doteq e^{f(x)}-1, \label{eqn:ICs_nlwa}
        \\
        u_{t}|_{t=0} &= G(x) \doteq g(x)e^{f(x)} \label{eqn:ICs_nlwb}. 
    \end{align}
$$
Since $f, g$ are supported in $B(0,R)$, we find that $F, G$ are supported in this ball as well. In summary, we have shown that if $v(x,t)$ solving \eqref{eqn:nlw_kn} exists then its corresponding $u(x,t)$ satisfies \eqref{eqn:wave_ivp} with $F, G$ supported inside $B(0,R)$. By reversing the computations performed above, we can likewise show that if 

- $$ u(x,t) $$ solves \eqref{eqn:wave_ivp} with the particular $$F, G$$ defined in \eqref{eqn:ICs_nlwa}, \eqref{eqn:ICs_nlwb} and

- $$ \|u(\cdot, t)\|_{C_{b}\left(\mathbb{R}^3\right)}<1 $$,

then $v$ defined by \eqref{eqn:u_to_v} solves \eqref{eqn:nlw_kn}. This gives us a recipe for constructing $v(x,t)$ explicitly by solving the linear wave equation and then bounding its solution!

&nbsp; &nbsp; &nbsp; &nbsp; So, let's define $u$ to be the solution of \eqref{eqn:wave_ivp} guaranteed to exist (and decay in time) by our earlier proposition. We now split into two cases depending on the size of $|t|$. 

#### <b>Case 1: $|t|\geq 1$.</b> 

In this case we see straight from \eqref{eqn:wave_eqn_decay} that $$\|u(\cdot, t)\|_{C_{b}\left(\mathbb{R}^3\right)}<1$$ provided
$$
    \begin{equation}
        \label{eqn:smallness_condition}
        \left\|F\right\|_{C^{1}_{b}\left(\mathbb{R}^3\right)} + \left\|G\right\|_{C_{b}\left(\mathbb{R}^3\right)} < \frac{1}{R^2}. 
    \end{equation}
$$

#### <b>Case 2: $|t|< 1$.</b>

Here we can't use \eqref{eqn:wave_eqn_decay} directly. Instead, we go back to the intermediary result \eqref{eqn:important_intermediate} in the proof of proposition \ref{prop:wave_eqn_decay}, which tells us that $$\|u(\cdot, t)\|_{C_{b}\left(\mathbb{R}^3\right)}<1$$ provided
$$
    \left\|F\right\|_{C^{1}_{b}\left(\mathbb{R}^3\right)} + \left\|G\right\|_{C_{b}\left(\mathbb{R}^3\right)} < 1. 
$$

Putting both these cases together and using that $R>1$, we find that $\|u(\cdot, t)\|_{C_{b}\left(\mathbb{R}^3\right)}<1$ for all $t\neq 0$ *if we impose the smallness condition \eqref{eqn:smallness_condition}*. So, provided \eqref{eqn:smallness_condition} holds, we can define $v(x,t)$ by \eqref{eqn:u_to_v} and get a global-in-time solution to \eqref{eqn:nlw_kn} as outlined above. $\square$

&nbsp; &nbsp; &nbsp; &nbsp; Notice again how the size of the initial data played into our analysis, but this time for the better. For most interesting nonlinear evolution equations, small initial data are unlikely to give rise to blowup, and they also tend to keep useful norms of the initial data "nice". For this reason, you'll often go to modern PDE talks where the speaker focuses entirely on small initial data. 

{% bibliography --cited %}
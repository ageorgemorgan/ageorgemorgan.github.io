---
math: true
title: "Analytic Issues with Laplace's Equation Part 2: Hadamard's Example of Ill-Posedness"
date: 2025-06-09
categories: [Math, Analysis of PDEs, Elliptic Theory]
tags: [elliptic equations, regularity theory]     # TAG names should always be lowercase
---

## Introduction 
In [an earlier post](https://ageorgemorgan.github.io/posts/laplace-issues/) (henceforth called simply "Part 1"), we learned that the geometry of a domain $\Omega$ affects the regularity of the solution to the Dirichlet problem for Laplace's equation on $\Omega$. In this follow-up post, we discuss a classical example due to Hadamard showing that boundary value problems for the Laplace equation can actually tbe <i>ill-posed</i>: while a unique solution exists, arbitrarily small changes to the boundary conditions lead to enormous changes in the solution! As with Part 1, complete details can be understood using nothing more sophisticated than separation-of-variables.  

## Statement of Hadamard's Problem
Consider the following boundary value problem on the entire plane:
$$
    \begin{equation}
        \label{eqn:BVP}
        \left\{
            \begin{aligned}
                \Delta u &= 0 \quad \forall \ (x,y)\in \mathbb{R}^2
                \\
                u|_{y=0} &= 0
                \\
                u_{y}|_{y=0} &= \frac{1}{n}\sin(nx).
                \end{aligned} 
        \right.
    \end{equation}
$$
This means we're imposing wave-equation-like initial conditions on Laplace's equation ($y$ is acting like a time variable). Thus we could justifiably call the above problem a Cauchy (initial value) problem rather than a boundary value problem. 

&nbsp;&nbsp;&nbsp;&nbsp; Recall that \eqref{eqn:BVP} is said to be <b>well-posed</b> if the following conditions are met (for more discussion see for instance {% cite Choksi2022 %}, p.10):
<ol>
    <li>
        a solution $u(x,y)$ exists,    
    </li>
    <li>
        $u(x,y)$ is <i>unique</i>, and 
    </li>
    <li>
        $u(x,y)$ depends continuously on the boundary data. 
    </li>
</ol>

## Solution via Separation of Variables
Assume 
$$
    u(x,y) = X(x)Y(y). 
$$
Plugging into Laplace's equation gives 

$$
    0 = X''Y + XY'' \Rightarrow \frac{X''(x)}{X(x)} = \frac{-Y''(y)}{Y(y)}. 
$$

The above expression can only hold if there exists a separation constant $\mu\in\mathbb{R}$ so that 

$$
    \begin{equation}
         \frac{X''(x)}{X(x)} = \frac{-Y''(y)}{Y(y)} = \mu. 
    \end{equation}
$$

Thus we have converted the PDE into two ODEs:
$$
            \begin{align}
                X'' &= \phantom{-}\mu X, \label{eqn:ODE_x}
                \\
                Y'' &= -\mu Y. \label{eqn:ODE_y}
            \end{align}
$$
We must now consider two cases depending on the sign of the separation constant $\mu$. 

### Case 1: $\mu>0$
Right away, we know that $\mu=\lambda^2$ for $\lambda>0$. Then, the general solutions of \eqref{eqn:ODE_x}, \eqref{eqn:ODE_y} are 

$$
    \begin{align*}
        X(x) &= A_{1}\cosh\left(\lambda x\right) + A_{2}\sinh\left(\lambda x\right),
        \\
        Y(y) &= B_{1}\cos\left(\lambda y\right) + B_{2}\sin\left(\lambda y\right).
    \end{align*}
$$

Now, $u|_{y=0}=0$ immediately implies $B_{1}=0$. For the other condition, we compute 

$$
    \begin{align*}
        u_{y}|_{y=0} &= X(x)Y'(0) 
        \\
        &= \left[A_{1}\cosh\left(\lambda x\right) + A_{2}\sinh\left(\lambda x\right)\right]\lambda B_{2}
        \\
        &= \frac{1}{n} \sin\left(nx\right). 
    \end{align*}
$$

However, $\cosh(\lambda x), \sin(\lambda x)$, and $\sin(nx)$ are linearly independent, so there is no way the above expression can hold. We conclude that this case ($\mu>0$) is impossible. 

### Case 2: $\mu \leq 0$
Write $\mu=-\lambda^2$ for $\lambda\geq0$. The ODEs in \eqref{eqn:ODE_x}, \eqref{eqn:ODE_y} have general solutions 

$$
    \begin{align*}
            X(x) &= A_{1}\cos\left(\lambda x\right) + A_{2}\sin\left(\lambda x\right),
            \\
            Y(y) &= B_{1}\cosh\left(\lambda y\right) + B_{2}\sinh\left(\lambda y\right).
     \end{align*}
$$

$u|_{y=0}=0$ implies $B_{1}=0$, and the other condition yields 

$$
    \begin{align*}
        u_{y}|_{y=0} &= X(x)Y'(0) 
        \\
        &= \left[A_{1}\cos\left(\lambda x\right) + A_{2}\sin\left(\lambda x\right)\right]\lambda B_{2}
        \\
        &= \frac{1}{n} \sin\left(nx\right). 
    \end{align*}
$$

Using linear independence again, the above holds provided

$$   
    \begin{align*}
        A_{1}&=0,
        \\
        \lambda A_{2}B_{2} &= \frac{1}{n}, \ \text{and}
        \\
        \lambda &= n . 
    \end{align*}
$$

Consequently, the solution here is 

$$
    \begin{align}
    u(x,y) =X(x)Y(y) &= A_{2}B_{2}\sin(\lambda x)\sinh(\lambda y) \nonumber
    \\ & = \frac{1}{n^2} \sin(n x)\sinh(n y) \label{eqn:PDE_soln}.
    \end{align}
$$

So, the solution to \eqref{eqn:BVP} must be given by \eqref{eqn:PDE_soln}. By construction, this solution is unique. 

## The Emergence of Ill-Posedness
To see that our work implies \eqref{eqn:BVP} is ill-posed, we investigate how our solution behaves as $n\rightarrow \infty$. For notational clarity, let us define 
$$
    \begin{equation}
    u_{n}(x,y) \doteq \frac{1}{n^2} \sin(n x)\sinh(n y).
    \end{equation}
$$

Fix any $(x,y)\in \mathbb{R}^2$ such that $y \neq 0$ and 

$$
    \frac{x}{\pi} \notin \mathbb{Q}. 
$$

If you want a concrete example, pick $(x,y)=(\pi^2, 1)$. In particular, this choice of $x$ means 

$$
    |\sin(nx)|>0 \quad \forall \ n\geq 1. 
$$

Then, we find 

$$
    \begin{align*}
        \lim_{n\rightarrow\infty}|u_{n}(x,y)| &= \lim_{n\rightarrow \infty} n^{-2}\left|\sinh(ny)\right||\sin(nx)|
        \\
        &= \frac12 \lim_{n\rightarrow \infty} n^{-2} \ e^{n|y|} \ |\sin(nx)|
        \\
        &= +\infty,
    \end{align*}
$$

since $e^{n|y|}$ blows up faster than $n^{2}$ (for instance, by L'Hopital's rule) and sine is bounded. We conclude that

$$
    \lim_{n\rightarrow\infty}|u_{n}(x,y)| = +\infty
$$

on the complement of 

$$
    \mathcal{S} = \left\{y= 0 \ \text{or} \ \frac{x}{\pi}\in \mathbb{Q}\right\}. 
$$

In the jargon of measure theory, $\mathcal{S}$ is said to be a <b>null set</b>: it has Lebesgue measure zero and therefore covers zero area (think if $\mathcal{S}$ as a collection of isolated dust particles lying on top of a piece of paper). So, as $n\rightarrow \infty$, $u_{n}(x,y)$ blows up almost everywhere: we cannot even make sense of the ``pointwise limit" as an element of $L^{\infty}\left(\mathbb{R}^2\right)$, the space of Lebesgue-measurable functions that are finite almost everywhere!

&nbsp;&nbsp;&nbsp;&nbsp; Why does this blowup imply ill-posedness of \eqref{eqn:BVP}? First notice that the Cauchy data obeys

$$
    \frac{1}{n} \sin(nx) \rightarrow 0 \ \text{uniformly as} \ n\rightarrow \infty. 
$$

Accordingly, the limiting case of \eqref{eqn:BVP} is 

$$
    \begin{equation}
    \label{eqn:BVP_inf}
      \left\{
        \begin{aligned}
        \Delta u &= 0 \quad \forall \ (x,y)\in \mathbb{R}^2
        \\
        u|_{y=0} &= 0
        \\
        u_{y}|_{y=0} &= 0
        \end{aligned} ,
        \right.
    \end{equation}
$$

which admits the trivial solution $u\equiv 0$. However, $u_{n}(x,y)$ obviously cannot converge to $$u\equiv 0$$ in any reasonable sense. Accordingly, there is no way the data-to-solution correspondence can be continuous. Therefore, \eqref{eqn:PDE_soln} is not well-posed!
 
&nbsp;&nbsp;&nbsp;&nbsp; In practice, it's useful to think of the continuity requirement in the well-posedness definition as a very weak version of <b>stability</b>. In this problem, we have effectively shown that the trivial solution $u\equiv 0$ of \eqref{eqn:BVP_inf} is unstable: if we pick $n\gg1$, then the Cauchy data in \eqref{eqn:BVP} is very close to zero (pointwise), but the solution $u_{n}$ is very far from $0$ (pointwise). 

{% bibliography --cited %}

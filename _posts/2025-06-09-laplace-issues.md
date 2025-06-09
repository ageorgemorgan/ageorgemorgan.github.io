---
math: true
title: "Analytic Issues with Laplace's Equation Part 1: Loss of Regularity"
date: 2025-06-09 +0500 
categories: [Math, Analysis of PDEs, Elliptic Theory]
tags: [elliptic equations, regularity theory]     # TAG names should always be lowercase
---

## Introduction 
Throughout your mathematical life, you'll often hear Laplace's equation described as one of the nicest PDEs out there: this equation posses an elegant maximum principle and mean-value property. Additionally, there are beautiful connections between harmonic functions on a domain $$\Omega\subseteq \mathbb{R}^2$$ and *holomorphic* functions on the same domain $$\Omega$$ viewed as a subset of $$\mathbb{C}$$. However, just like us humans, no PDE is perfect. Even Laplace's equation can exhibit some pathological behaviour. 

&nbsp;&nbsp;&nbsp;&nbsp; In this mini-series of posts, we'll explore some analytical problems that can occur when dealing with Laplace's equation. Today we discover that, if we do not pose Laplace's equation on a smooth domain, then the derivatives of our solution may behave quite badly near the domain's boundary. Interestingly, this nasty example can be completely understood using a simple separation-of-variables computation! 

## Loss of Regularity in the Wedge Problem
Fix a particular angle $\beta \in (0,2\pi)$ and a radius $a>0$. Let us define a (nonsmooth) domain $\Omega\subseteq \mathbb{R}^2$ using polar coordinates:
$$
  \begin{equation}
    \Omega \doteq \left\{0< r< a, \ 0<\theta<\beta\right\}.
  \end{equation}
$$
$$\Omega$$ therefore represents a wedge when $\beta\leq \pi$ or a "Pac-Man" shape when $\beta>\pi$. 

We then consider the boundary value problem of finding $$u\colon \overline{\Omega}\rightarrow\mathbb{R}$$ such that, for a given continuous function $$g(\theta),$$ 
$$
  \begin{equation}
    \label{eqn:wedge_BVP}
    \left\{
    \begin{aligned}
      \Delta u &= 0 \quad \forall \ (r,\theta)\in \Omega
      \\
      u|_{\theta=0} &=0
      \\
      u|_{\theta=\beta} &= 0
      \\
      u|_{r=a} &= g(\theta)
      \\
      \left|u(r,\theta)\right| &<\infty \quad \forall \ (r,\theta)\in \overline{\Omega}.
    \end{aligned} 
    \right.
  \end{equation}
$$
Solutions to this problem can be ill-behaved depending on the angle $$\beta$$ subtended by the wedge, as the question below demonstrates. 

### Question 
a) Solve \eqref{eqn:wedge_BVP} using polar coordinates. 

b) Show that, if $$\beta> \pi$$, then $u_{r}$ is not necessarily bounded on $$\overline{\Omega}$$ for arbitrary $$g(\theta)$$.

c) Put a restriction on $$g(\theta)$$ guaranteeing that $$u_{r}$$ is bounded on $$\overline{\Omega}$$.

### Solution

a) This first part of the solution mostly follows {% cite Strauss2008 %}, section 6.3. We employ separation of variables, looking first for solutions of the form 
$$
  \begin{equation}
    u(r,\theta) = R(r)\Theta(\theta). 
  \end{equation}
$$
Since the polar coordinate Laplacian takes the form 
$$
  \begin{equation}
    \Delta u = u_{rr}+r^{-1}u_{r}+r^{-2}u_{\theta\theta}. 
  \end{equation}
$$
we find that 
$$
\begin{equation}\label{eqn:separated}
-\frac{\Theta''}{\Theta} = \frac{r^2R''+rR'}{R} = \lambda^2
\end{equation}
$$
for a constant $$\lambda$$.

&nbsp;&nbsp;&nbsp;&nbsp; We focus on the angular equation first. $$\Theta$$ obeys the BVP
$$
\begin{equation}
\label{eqn:angular_BVP}
  \left\{
    \begin{aligned}
    \Theta'' &= -\lambda^2\Theta,
    \\
    \Theta|_{\theta=0} &= 0
    \\
    \Theta|_{\theta=\beta} &= 0
    \end{aligned} 
    \right.
\end{equation}
$$
Immediately, we know the solution to this problem is 
$$
  \begin{equation}
    \Theta_{n}(\theta) = \text{constant}\times\sin\left(\frac{n\pi}{\beta}\theta\right), \quad \lambda=\lambda_{n} = \frac{n\pi}{\beta}, \ n=1,2,3,...
  \end{equation}
$$
&nbsp;&nbsp;&nbsp;&nbsp; Now that we've found the eigenvalues, we can turn to the radial equation: 
$$
\begin{equation}
\label{eqn:radial_BVP}
  \left\{
    \begin{aligned}
    r^2R''+rR' +\lambda_{n}^2R &= 0,
    \\
    |R| &< \infty.
    \end{aligned} 
    \right.
\end{equation}
$$
We look for solutions in the form $$R=r^{a}$$. Plugging into the above equation, we find that $$|a|=|\lambda_{n}|$$. To ensure $$|u|<\infty$$, we only take the positive values of $a$. This yields
$$
  \begin{equation}
    R_{n}(r) = \text{constant}\times r^{\frac{n\pi}{\beta}}. 
  \end{equation}
$$
Then, using the principle of superposition, we can express our solution $$u(r,\theta)$$ as a Fourier series with $$r$$-dependent coefficients: 
$$
\begin{equation}
\label{eqn:wedge_FS}
    u(r,\theta) = \sum_{n=1}^{\infty} A_{n} \ r^{\frac{n\pi}{\beta}} \ \sin\left(\frac{n\pi}{\beta}\theta\right). 
\end{equation}
$$
&nbsp;&nbsp;&nbsp;&nbsp; Next, we impose the boundary condition at $r=a$ to determine the expansion coefficients $$A_{n}$$:
$$
  \begin{equation}
    g(\theta) = u(a,\theta) = \sum_{n=1}^{\infty} A_{n} \ a^{\frac{n\pi}{\beta}} \ \sin\left(\frac{n\pi}{\beta}\theta\right). 
  \end{equation}
$$
Using trigonometric orthogonality relations, we know then that the $$n^{\text{th}}$$ Fourier sine series coefficient of $$g$$ on $$[0,\beta]$$ is 
$$ 
  \begin{equation} 
    A_{n} \ a^{\frac{n\pi}{\beta}} = \frac{2}{\pi}\int_{0}^{\beta} g(\theta) \ \sin\left(\frac{n\pi}{\beta}\theta\right)\ \mathrm{d} \theta.
  \end{equation}
$$
We conclude that
$$ 
\begin{equation}
    \label{eqn:wedge_fourier_coeffs}
    A_{n} = \frac{2a^{-n\pi/\beta}}{\pi}\int_{0}^{\beta} g(\theta) \ \sin\left(\frac{n\pi}{\beta}\theta\right)\ \mathrm{d} \theta.
\end{equation}
$$
Combining \eqref{eqn:wedge_FS} and \eqref{eqn:wedge_fourier_coeffs} gives our complete solution to \eqref{eqn:wedge_BVP}.

b) Using \eqref{eqn:wedge_FS} (and glossing over the technical details of the term-by-term differentiation) we find that 
$$
  \begin{equation}
    u_{r} = \sum_{n=1}^{\infty} A_{n} \ \frac{n\pi}{\beta} \ r^{\frac{n\pi}{\beta}-1} \ \sin\left(\frac{n\pi}{\beta}\theta\right). 
  \end{equation}
$$
The $$n=1$$ term includes $$r^{\frac{\pi}{\beta}-1}$$, which blows up as $$r\rightarrow 0^{+}$$ when $$\beta>\pi$$. Therefore, for such $$\beta$$, $$u_{r}$$ may be unbounded on $$\overline{\Omega}$$. 

c) We know that the $$n=1$$ term is what causes $$u_{r}$$ to become unbounded near $$r=0$$. So, by picking $$g(\theta)$$ so that $$A_{1}=0$$, this troublesome term doesn't appear. Concretely, we may use \eqref{eqn:wedge_fourier_coeffs} to express this constraint as 
$$
  \begin{equation}
    \label{eqn:finiteness_constraint}
    0=\int_{0}^{\beta} g(\theta) \ \sin\left(\frac{\pi}{\beta}\theta\right)\ \mathrm{d} \theta.
  \end{equation}
$$
If \eqref{eqn:finiteness_constraint} holds, then $$u_{r}$$ extends to a continuous, bounded function on $$\overline{\Omega}$$.


&nbsp;&nbsp;&nbsp;&nbsp; So, we have found that solutions to Laplace's equation on a nonsmooth domain $$\Omega$$ do not necessarily extend to smooth functions on $$\overline{\Omega}$$. In particular, the sharpness in $$u$$ (severe lack of boundedness in $$u_{r}$$) is inherited from the sharpness of the wedge (the sharp point at $$r=0$$). In the language of analysis, we say that our solution **loses regularity** when $$\beta$$ exceeds the threshold value $$\beta=\pi$$: think of "regular" as meaning "very very good derivatives". 

> When we demand that the underlying domain $$\Omega$$ and the boundary data $$g(\theta)$$ are smooth, it turns out that we can say quite a bit about smoothness of the corresponding harmonic function using the principle of **elliptic regularity**: see {% cite Evans2010 %} for details. 
{: .prompt-info }

{% bibliography --cited %}

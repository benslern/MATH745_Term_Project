# MATH745_Term_Project

For my term project I will compute and analyse solutions of the Generalized Constantine-Lax-Majda equation. The Constantine-Lax-Majda (CLM) equation was proposed as a one-dimensional model for the three-dimensional vorticity equation and has the form

$$\omega_t - u_x\omega = 0, \qquad u_x = H\omega$$

where $\omega$ is the vorticity of the fluid, $u$ is the fluid velocity. For this analysis the equations will be examined on the periodic domain $x\in[-\pi,\pi]$. The Hilbert transform of the vorticity is defined as

$$H\omega(x,t) = \frac{1}{2\pi}\int_{-\pi}^{\pi} \omega(y,t)\cot\left(\frac{x-y}{t}\right) \,dy,$$

which has a singularity at $x=y$ and so must be evaluated using the Cauchy Principle Value. The CLM equation was later expanded on by De Gregorio to include the convection term $u\omega_x$. The De Gregorio equation has the form

$$\omega_t + u\omega_x - u_x\omega = 0, \qquad u_x = H\omega.$$

The De Gregorio equation was later generalized by Okamoto et al to include the real parameter $a$. The result is the Generalized Constantin-Lax-Majda equation (GCLM)

$$\omega_t + au\omega_x - u_x\omega = 0, \qquad u_x = H\omega.$$

The parameter $a$ is a real value that determines the relative strength of the convection term. When $a=0$ the GCLM is equal to the CLM, and when $a=1$ it is equal to the De Gregorio equation.

The CLM equation was proposed in 1985 and has been shown to exhibit finite time blow-up. While the De Gregorio equation was proposed shortly after in 1989, questions still remain as to whether or not the equation permits finite time blow-up. The goal of this analysis will be to examine the affect the parameter $a$ has on the finite time blow-up of the GCLM equation. I will explore solutions of the GCLM for values of $a\in[-1,1]$. 

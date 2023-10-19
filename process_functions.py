
# Dissolution_process.py

"""
Computes the dissolution (temporal) profile of a bubble.

Freely inspired from :
* [Venerus, JNNFM (2015)](https://www.sciencedirect.com/science/article/pii/S037702571400202X)
* [Fyrillas _et al_. Langmuir (2000)](https://pubs.acs.org/doi/abs/10.1021/la990784y)
* [Kloek, Meinders and van Vliet, JCIS (2001)](https://www.sciencedirect.com/science/article/pii/S0021979701974545)
* [Macosko, 1994](https://www.wiley.com/en-sg/Rheology:+Principles,+Measurements,+and+Applications-p-9780471185758) for the notations

I consider :
* a neo-Hookean elastic medium, with ${\bf T} = -p {\bf I} + G {\bf B}$ under yielding, and a purely plastic medium above it 
* a yielding criterion that is von Mises, as used in Venerus 2015 
* a yielded region of size $S$, as defined by Venerus 2015 at Equation (20).
* plastic stresses above yielding in the form of $2 \sqrt{3} \sigma_Y \ln (R/S)$ with $S$ the size of the yielded region (Eq. 22 from Venerus 2015)
* a surface tension that is constant with the bubble radius $R(t)$, i.e. no surfactants or anything

"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

### BUBBLE DISSOLUTION CLASS ########

class bubble_params:
    """
    A "bubble_params" class to quickly compute relevant quantities
    Some values are "baked in", e.g. ambient pressure, Henry's law constants
    or diffusion coefficient. Check the values out.

    CONSTRUCTOR : 
    - bub = bubble_params() 

    NOTES : 
    - Temperature is set to 295 by default
    - Surface tension set to 0.072 by default
    
    """
    def __init__(self, R=100e-6, G=0, yieldstress=np.inf, surftens=0.072, surfelast=0.0, T=295.0, f=1.0, 
                 p0=1.013e5, ignore_transient=True, initial_stretch=False):
        self.r0 = R              # Initial radius
        self.R0 = R
        self.T = T             # Temperature
        self.surftens= surftens     # Surface tension of the material (sigma)
        self.surfelast = surfelast # Surface elasticity of the material (= 1/sigma d(sigma)/dR if we follow gibbs)
        self.yieldstress = yieldstress # Yield stress of the material
        self.G = G               # Shear modulus of the gel (in linear regime)
        self.f = f                  # Undersaturation
        self.p0 = p0            # Ambient pressure

        self.ignore_transient = ignore_transient # Ignore the diffusion layer buildup at beginning of dissolution
        self.initial_stretch = initial_stretch # Enforces that the initial state of the bubble corresponds to a "stretched" surrounding material 
                                               # instead of it being at rest

        self.fit_message = ''   # Will be used later on to say if fit to experimental data was successful
        self.kHn2 = 6.1e-6
        self.kHo2 = 1.3e-5       # Henry's constant of O2 in water
        self.Mn2 = 28e-3
        self.Mo2 = 32e-3         # Molar mass (SI) of O2
        self.D = 2.1e-9            # Diffusion coefficient of air in water
        self.Rg = 8.314          # Ideal gas constant    
        self.xo2 = 0.21
        self.xn2 = 0.79

        self.time_ODEfit = np.array([])
        self.r_ODEfit = np.array([])

        # Things that will be computer later on ...
        self.kHeff = (self.xo2*self.kHo2*self.Mn2 + self.xn2*self.kHn2*self.Mo2)/(self.xo2*self.Mo2 + self.xn2*self.Mn2)

    def tcarac(self, r0=None):
        if r0 == None:
            r0 = self.r0
        return 4/self.kHeff*r0**2/(3*self.D*self.Rg*self.T)

    def Ca(self, surftens=None, r0=None):
        if surftens == None:
            surftens = self.surftens
        if r0 == None:
            r0 = self.r0
        return r0*self.p0/(2*surftens)

# 
def model_dissolution_Venerus(params):
    """ Models the dissolution of a bubble in a YS medium __including elastic and plastic contributions
    to strain__ (due to geometry), surface tension terms but ignoring viscous effects in the process
    since I arbitrarily judge that these are small as the strain rates in bubble dissolution are mostly of
    order 1e-3 1/s. """

    # Dimensionless groups
    dimless = params.kHeff*params.Rg*params.T
    f = params.f
    
    # Compute the "functions" that return the size of the yielded region,
    # the pressure inside the bubble and its derivative with respect to R/R0
    # for an input R (and, of course, rheological parameters and the initial bubble size R0) 
    # once and for all (it is a bit time consuming)
    S, _ = S_fun(params)
    Pterm = P_fun(params, S)
    dPdRterm = dPdR_fun(params, S)

    def _dissolution_ode_venerus_eps(x, eps:float):
        """ This returns what is equal to -dot{eps} in the 
        equation governing the evolution of the bubble radius as 
        a function of time in the model of Venerus (but borrowing
        multiple conventions from Kloek and Fyrillas)
        """

        # Functions and their current values
        P = Pterm(eps)   
        dPdR = dPdRterm(eps)
        deps = -6*dimless*(P - f)/(3*P - eps*dPdR)*(x/eps + np.pi**(-0.5)) 
        return deps

    def _dissolution_ode_venerus_eps2(x:float, eps2:float):
        """ This is basically the differential equation for epsilon² as a function of x,
        I have tried to use that to regularize even more the problem of bubble dissolution,
        but in the end it does not change anything. You can use it in `solve_ivp` below and adapt
        the `r_ODEfit` to give R(t) again to check for yourself. """
        # Functions and their current values
        P = Pterm(eps2)   
        dPdR = dPdRterm(eps2) # Should it be different ?
        if eps2 > 0:
            deps2 = -12*dimless*(P - f)/(3*P - 2*eps2*dPdR)*(x + (eps2/np.pi)**(-0.5)) 
        else: 
            deps2 = 0*eps2
        return deps2
    
    def _dissolved(x:float, eps2:float):
        """ An event function examining whether the bubble is dissolved"""
        return eps2
    _dissolved.direction = -1
    _dissolved.terminal = True

    fit_up_to = params.tcarac()*100
    t_eval = np.linspace(0, fit_up_to, 3000)
    x_eval = np.sqrt(t_eval*params.D)/params.r0

    solution = solve_ivp(_dissolution_ode_venerus_eps, (0, np.max(x_eval)), [1], 
                         method='RK45', t_eval=x_eval, atol=1e-7, rtol=1e-4, events=_dissolved) # Solution of the type R = f(sqrt(t))

    t_dissolved = np.squeeze(solution.t_events)
    if not t_dissolved:
        t_dissolved = np.nan

    params.time_ODEfit = np.append(solution.t, t_dissolved)**2*params.r0**2/(params.D)
    params.r_ODEfit = np.append(solution.y[0], 0)*params.r0
    params.x_ODEfit = solution.t
    params.y_ODEfit = solution.y[0]
    params.tdiss_ODEfit = solution.t_events[0]**2*params.r0**2/(2*params.D)
    params.mess_ODEfit = solution.message
    
    return params


def S_fun(params): 
    """ Computes the extent of the yielded region, for a given R_0 and R
        (see Venerus, JNNFM 2015)

    ARGS
    -----
    * params : `bubble_params` object with attributes R0, G, yieldstress
    * Rs [list or ndarray] : the radii at which you want to compute S_fun

    RETURNS
    -------
    * S_fun [ndarray] : the extent of the yielded region for the Rs considered
    NOTE : if there is no solution of the equation for S_fun, the code automatically 
    sets S_fun = R, which means the plastic term is set to 0. 
    """

    G = params.G
    sigmaY = params.yieldstress

    if G == 0:
        S_fun = lambda R : R
        dS_fun = lambda R : 1
        return S_fun, dS_fun

    Rs = np.concatenate((np.logspace(-5,-0.5,500), \
                         np.linspace(0.317,3.16,10000), \
                         np.logspace(0.5,3,500)))
    Ss = np.zeros_like(Rs)
    vfun = lambda R, S_fun : ((1 - (R**3 - 1**3)/S_fun**3)**(+4/3) - 1)**2 \
                     + 2*((1 - (R**3 - 1**3)/S_fun**3)**(-2/3) - 1)**2 \
                     - 2*(sigmaY/G)**2

    for no, R in enumerate(Rs):
        Ss[no] = fsolve(lambda S_fun : vfun(R, S_fun), x0=R)[0]

    condition = (np.abs(vfun(Rs, Ss)) > 1e-6) | (Ss < Rs)
    Ss[condition] = Rs[condition]

    diffS = np.diff(Ss)/np.diff(Rs)
    Rs_diffS = 0.5*(Rs[:-1] + Rs[1:])

    S_fun = interp1d(x=Rs, y=Ss, kind='linear', bounds_error=False, fill_value='extrapolate') 
    dS_fun = interp1d(x=Rs_diffS, y=diffS, kind='linear', bounds_error=False, fill_value='extrapolate')

    return S_fun, dS_fun

def G_fun(params, S_fun, dimensionless=True):
    """ G_FUN : ELASTIC part of the matrix resistance to 
    a bubble dissolution. 

    ARGS
    ----
    * params : the parameters of your bubble (needs to include attributes `yieldstress`and `p0` for 
    the dimensionless case)
    * S_fun : the interpolator (obtained from `S_fun`) giving the extent of the yielded region
    as a function of R for a given R0
    * dimensionless [default True] : do you want a dimensionless result or a result in Pa ?

    RETURNS
    -----
    * G_fun : the interpolator giving the value of the ELASTIC stress if you give it a radius 
    (e.g. if you call `G_fun(1e-3)`, it will output you a number)

    NB : dimensionless means we plot the G/p0 term = f(R/R0)
    non dimensionless terms means we plot G term = f(R/R0) still """

    Rs = np.concatenate((np.logspace(-5,-0.5,500), \
                         np.linspace(0.317,3.16,10000), \
                         np.logspace(0.5,3,500)))
    Svals = S_fun(Rs)
    Gvals = 5/2 - 1/2*(1 - ((Rs**3 - 1**3)/(Svals**3)))**(1/3)*(5 - ((Rs**3 - 1**3)/Svals**3))

    if not dimensionless: 
        Gvals *= params.G
    else:
        Gvals *= params.G/params.p0
    
    Gfun = interp1d(x=Rs, y=Gvals, kind='linear', bounds_error=False, fill_value='extrapolate')
    return Gfun

def YS_fun(params, S_fun, dimensionless=True):
    """ YS_FUN : ELASTIC part of the matrix resistance to 
    a bubble dissolution. 

    ARGS
    ----
    * params : the parameters of your bubble (needs to include attributes `yieldstress`and `p0` for 
    the dimensionless case)
    * S_fun : the interpolator (obtained from `S_fun`) giving the extent of the yielded region
    as a function of R for a given R0
    * dimensionless [default True] : do you want a dimensionless result or a result in Pa ?

    RETURNS
    -----
    * YS_fun : the interpolator giving the value of the PLASTIC stress if you give it a radius 
    (e.g. if you call `YS_fun(1e-3)`, it will output you a number)

    NB : dimensionless means we plot the YS/p0 term as a function of (R/R0)
    non dimensionless terms means we plot YS term  as a function of (R/R0) still """

    Rs = np.concatenate((np.logspace(-5,-0.5,500), \
                         np.linspace(0.317,3.16,10000), \
                         np.logspace(0.5,3,500)))
    
    Svals = S_fun(Rs)
    YSvals = 2*np.sqrt(3)*np.log(Svals/Rs)
    YSvals[Rs < 1] = -YSvals[Rs < 1]

    if not dimensionless:
        YSvals *= params.yieldstress
    else: 
        YSvals *= params.yieldstress/params.p0

    YSfun = interp1d(x=Rs, y=YSvals, kind='linear', bounds_error=False, fill_value='extrapolate')

    return YSfun

def P_fun(params, S_fun):
    """ P_FUN : Generates a function computing the pressure inside a bubble 
    for a given R and initial conditions (given by a `params` object)

    ARGS
    ----
    * params : the parameters of your bubble (needs to include attributes `yieldstress`and `p0` for 
    the dimensionless case)
    * S_fun : the interpolator (obtained from `S_fun`) giving the extent of the yielded region
    as a function of R for a given R0

    RETURNS
    -----
    * P_fun : the interpolator giving the value of the pressure inside the bubble ; it will, e.g. give 
    you some value if you call `P_fun(0.3)`
    
    NOTE: for R/R0 < 1e-5, it might start to break since the interpolator extrapolates linearly out of bounds,
    but that would mean R < 1e-7 mm unless your bubble is initially _really_ big, so it shouldn't be an issue.
    """

    delta_s = 2*params.surftens/params.R0/params.p0  # Dimensionless surface tension
    Gterm = G_fun(params, S_fun)
    YSterm = YS_fun(params, S_fun)
    return lambda eps: 1 + Gterm(eps) + YSterm(eps) + delta_s/eps # Do not forget the p0 term ... !

def dPdR_fun(params, S_fun):
    """ dPdR_FUN : Generates a function computing the differential of p(R)
    for a given R and initial conditions (given by a `params` object)
    (actually, it computes dP(R/R0) /d(R/R0), as all radii in the contributions 
    to p are computed as functions of R/R0). We stupidly compute it from P_fun by evaluating 
    it at plenty of places ; then we differentiate the results and build an interpolator to allow
    the user to compute dP/dR at any value. 

    NOTE: for R/R0 < 1e-5, it might start to break since the interpolator extrapolates linearly out of bounds,
    but that would mean R < 1e-7 mm unless your bubble is initially _really_ big, so it shouldn't be an issue.

    ARGS
    ----
    * params : the parameters of your bubble (needs to include attributes `yieldstress`and `p0` for 
    the dimensionless case)
    * S_fun : the interpolator (obtained from `S_fun`) giving the extent of the yielded region
    as a function of R for a given R0

    RETURNS
    -----
    * dPdR_fun : the interpolator giving the value of dP/dR inside the bubble ; it will, e.g. give 
    you some value if you call `dPdR_fun(0.3)`
    """

    Pterm = P_fun(params, S_fun)
    Rs = np.concatenate((np.logspace(-5,-0.5,500), \
                         np.linspace(0.317,3.16,10000), \
                         np.logspace(0.5,3,500)))
    
    Rs_diff = 0.5*(Rs[1:] + Rs[:-1])
    dPvals = np.diff(Pterm(Rs))/np.diff(Rs)

    return interp1d(x=Rs_diff, y=dPvals, kind='linear', bounds_error=False, fill_value='extrapolate')


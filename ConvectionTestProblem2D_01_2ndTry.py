from fipy import CellVariable, Gmsh2D, TransientTerm, DiffusionTerm, Viewer, serialComm, FaceVariable, ImplicitSourceTerm
from fipy.tools import numerix

from IPython import embed
import sys
import logging
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

from oct2py import octave

#head
def get_solution_along_ray(soln,R_inner,R_outer,cellSize,rayAngle,N_pointsAlongRay):
    '''
    define a ray from origin outward at angle rayAngle.
    obtain solution values on this ray via interpolation.
    avoid nan's by staying within the convex hull of mesh points
    '''

    assert rayAngle>(cellSize/R_inner), 'rayAngle has to be large enough for ray to end up within convex hull of grid pts'
    assert rayAngle<(np.pi-cellSize/R_inner), 'rayAngle has to be smaller than pi so that ray ends up within convex hull of grid pts'    
    x,y=soln.mesh.cellCenters
    z=soln.value

    #change to rotated coordinate system
    xp=np.cos(rayAngle)*x+np.sin(rayAngle)*y
    yp=-np.sin(rayAngle)*x+np.cos(rayAngle)*y

    #want solution along ray yp=0
    #stay within the convex hull of grid points
    s_for_interp=np.linspace(R_inner+cellSize,R_outer-cellSize,N_pointsAlongRay)

    xp_Plaid,yp_Plaid=np.meshgrid(s_for_interp,np.linspace(0,0,1))

    soln_Plaid=scipy.interpolate.griddata((xp, yp), z, (xp_Plaid, yp_Plaid), method='linear')  #SETTING

    solnAlongRay=soln_Plaid.ravel()
    
    return s_for_interp,solnAlongRay
    
#head
logLevel='INFO'
logFilename='convectionTestProblem2D_01.log'
logging.basicConfig(filename=logFilename,filemode="w",level=logLevel,format='%(asctime)s - %(levelname)s - %(message)s')
del logLevel,logFilename

R_outer=2e-6  #SETTING  must match what is in .geo file
R_inner=1e-6  #SETTING  must match what is in .geo file
cellSize=.05*(R_outer-R_inner)  #SETTING  must match what is in .geo file

#run gmsh gui to produce this mesh (.msh file)
#gmsh infiniteCylinder01.geo
filename='infiniteCylinder01.msh'
mesh=Gmsh2D(filename,communicator=serialComm)
del filename

T_initial=425.08  #deg K
T_infinity=293.15  #deg K  #SETTING

var=CellVariable(mesh=mesh,value=T_initial-T_infinity)  #let var be T-T_infinity
    
rho=6980.  #kg/m^3
cp=227.  #J/kg/K
k=59.6  #W/m/K
D_thermal=k/rho/cp

X_faces,Y_faces=mesh.faceCenters

surfaceFaces=(Y_faces>0) & ((X_faces**2+Y_faces**2)**.5 > R_outer-cellSize/10.)

#convectionCoeff=200.  #W/m^2/K
Bi_desired=10.  #SETTING
convectionCoeff=Bi_desired*k/R_inner
logging.info('convection coefficient is %.2E' % convectionCoeff)

Bi=convectionCoeff*R_inner/k  #Biot number
logging.info('Biot number is %.2E' % Bi)

#ref: https://github.com/usnistgov/fipy/blob/develop/documentation/USAGE.rst#applying-robin-boundary-conditions
#warning: avoid confusion between convectionCoeff in that fipy documentation, which refers to terms involving "a", and convectionCoeff here, which refers to a heat transfer convection coefficient at a boundary
Gamma0=D_thermal
Gamma = FaceVariable(mesh=mesh, value=Gamma0)
mask=surfaceFaces
Gamma.setValue(0., where=mask)
dPf = FaceVariable(mesh=mesh, value=mesh._faceToCellDistanceRatio * mesh.cellDistanceVectors)
Af = FaceVariable(mesh=mesh, value=mesh._faceAreas)
#RobinCoeff = (mask * Gamma0 * Af / (dPf.dot(a) + b)).divergence  #a is zero in our case
b=1.
RobinCoeff = (mask * Gamma0 * Af / b).divergence  #a is zero in our case
#eq = (TransientTerm() == DiffusionTerm(coeff=Gamma) + RobinCoeff * g - ImplicitSourceTerm(coeff=RobinCoeff * mesh.faceNormals.dot(a)))   #a is zero in our case
# g in this formulation is -convectionCoeff/k*var, where var=T-T_infinity
eq = (TransientTerm() == DiffusionTerm(coeff=Gamma) + ImplicitSourceTerm(RobinCoeff * -convectionCoeff/k))

# embed()
# sys.exit()

#either show the fipy viewer or make a T(s,t) plot; doing both at once is too much of a hassle to debug
showViewer=False  #SETTING
if showViewer:
    #viewer=Viewer(vars=var,datamin=T_infinity,datamax=T_initial)
    viewer=Viewer(vars=var)

    viewer.plot()

dt_explicit=cellSize**2/2./D_thermal

tFinal_nd=2.  #SETTING  dimensionless final temperature
tFinal=R_inner**2*tFinal_nd/D_thermal

# steps=401
# dt=tFinal/steps

dt=.9*dt_explicit
steps=int(tFinal/dt)

val=dt/dt_explicit
print 'dt over explicit limit dt is %.2E' % val
logging.info('dt over explicit limit dt is %.2E' % val)
del val

N_plots_desired=7  #SETTING
plotEvery=steps/N_plots_desired  #SETTING

showAnalytical=True

N_pointsAlongRay=10  #SETTING
rayAngle=np.pi/4.  #SETTING

if not showViewer:
    plt.ioff()
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel('r, m')
    ax.set_ylabel('temperature, K')
    ax.set_title('T(r) at selected times; black is FiPy and red is analytical solution')
    
    #plot the initial condition
    s_ray,var_ray=get_solution_along_ray(soln=var,R_inner=R_inner,R_outer=R_outer,cellSize=cellSize,rayAngle=rayAngle,N_pointsAlongRay=N_pointsAlongRay)
    ax.plot(s_ray,var_ray+T_infinity,'k.')

A=15  #SETTING for analytical solution

for step in range(steps):

    eq.solve(var=var,dt=dt)

    if (step>0) and (step % plotEvery==0):
        if showViewer:
            viewer.plot()
            print 'step is %d' % step
            #time.sleep(1.5)
            raw_input('hit enter to continue')
            
        else:
            s_ray,var_ray=get_solution_along_ray(soln=var,R_inner=R_inner,R_outer=R_outer,cellSize=cellSize,rayAngle=rayAngle,N_pointsAlongRay=N_pointsAlongRay)
            ax.plot(s_ray,var_ray+T_infinity,'k.')

            print 'step is %d; min T along ray is %.2E, max T along ray is %.2E' % (step,var_ray.min()+T_infinity,var_ray.max()+T_infinity)

        if showAnalytical and (not showViewer):
            t=dt*(step+1)
            tv=t*D_thermal/R_inner**2  #dimensionless

            Td_a=[]  #dimensionless temperature
            for s in s_ray:
                rv=s/R_inner
                Rprime=R_outer/R_inner
                #Td_a.append(octave.fdR23B01T0(rv=rv,tv=tv,R=Rprime,Bi=Bi,A=A))  #this throws an error
                Td_a.append(octave.fdR23B01T0(rv,tv,Rprime,Bi,A))                

            Td_a=np.array(Td_a)
            T_a=Td_a*(T_infinity-T_initial)+T_initial  #dimensional temperature
            ax.plot(s_ray,T_a,'r.')

if not showViewer:
    plt.show()

#raw_input('hit enter to continue')

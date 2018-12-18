# Convection2DFiPyExample01
This repo contains a series of efforts at numerically solving a simple 2D heat conduction problem having a simple convection boundary condition applied to a curved boundary, using **FiPy**.

The intention is to figure out how to apply such a boundary condition in FiPy; until now, it was not obvious or clear to me how this is done.

[FiPy](https://www.ctcms.nist.gov/fipy/) is a Python-based solver for partial differential equations, created and maintained by a [team](https://www.ctcms.nist.gov/fipy/documentation/CREDITS.html) at NIST.

The following is a sketch of the transient heat conduction problem:

[problem sketch](reportImages/2018-11-16_22-15-18Report.jpg)

The matching analytical solution is obtained from the [UNL Exact Analytical Conduction Toolbox](http://exact.unl.edu/exact/contents/display.php?eqtype=Heat%20Equation,%201D%20Hollow%20Cylinder&&name=Hollow%20cylinder%20with%20heating%20through%20convection%20at%20outer%20surface%20and%20insulated%20at%20inner%20surface) in the form of MATLAB code, which I have put in the form of scripts for use with [GNU Octave](https://www.gnu.org/software/octave/).

This repo contains a series of efforts to obtain agreement between that analytical solution, and a script that solves the problem via FiPy.  I interacted with the team at NIST using their [FiPy mailing list](https://www.ctcms.nist.gov/fipy/documentation/MAIL.html).  To find the relevant thread from late 2018, search the [mailing list archive](https://www.mail-archive.com/fipy@nist.gov/) for "how to apply convection boundary condition in 2D?".  Eventually it appears that the FiPy team succeeded in modifying their code as required.  Part of their work is found in [their forked repo](https://github.com/guyer/Convection2DFiPyExample01).  The resulting Python script is [here](ConvectionTestProblem2D_01_20181215Version.py).  This [report](Report20181216.pdf) shows good agreement between FiPy and the analytical solution for a particular set of problem input parameters, after several mesh refinements.

## License

GNU General Public License v3.0

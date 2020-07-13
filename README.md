## StabilizedRK

A collection of C++ classes that implement ODE solution methods designed to facilitate the efficient evolution of systems of ODE's to steady state.  The foundation for the methods are a family of "stabilized" Runge-Kutta methods, methods that are first order accurate, but whose region of absolute stability contains a large portion of the negative real axis. The method of determiniation of the coefficients of this family is described in the report 

Christopher R. Anderson and Christopher Elion, "Accelerated Solutions of Nonlinear Equations Using Stabilized Runge-Kutta Methods", CAM report 04-26, April 2004 (ftp://ftp.math.ucla.edu/pub/camreport/cam04-26.pdf

In addition to classes for constructing the coefficients of the stabilized Runge-Kutta methods there are classes that implement different adaptive timestepping techniques to facilitate rapid convergence to steady state. 


### Prerequisites
C++11
### Versioning
Release : 1.0.1
### Date 
July 13, 2020 
### Authors
Chris Anderson
### License
GPLv3  For a copy of the GNU General Public License see <http://www.gnu.org/licenses/>.
### Acknowledgements



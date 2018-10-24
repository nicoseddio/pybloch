import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import math
import cmath
import copy
import warnings
warnings.filterwarnings("ignore", message="Casting complex values to real discards the imaginary part")

class BlochSphere:
    def __init__(self, name = ""):
        #housekeeping
        self.title = "Bloch Sphere"
        if(len(name)>0):
            self.title = "Bloch Sphere: "+name
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.label_dist_mult = 1.1

        #sphere visuals
        # - outer circles
        circleX = plt.Circle((0,0),1,fill=False,edgecolor="b",alpha=0.1)
        circleY = plt.Circle((0,0),1,fill=False,edgecolor="b",alpha=0.1)
        circleZ = plt.Circle((0,0),1,fill=False,edgecolor="b",alpha=0.1)
        self.ax.add_patch(circleX)
        art3d.pathpatch_2d_to_3d(circleX, z=0, zdir="x")
        self.ax.add_patch(circleY)
        art3d.pathpatch_2d_to_3d(circleY, z=0, zdir="y")
        self.ax.add_patch(circleZ)
        art3d.pathpatch_2d_to_3d(circleZ, z=0, zdir="z")
        # - axis
        self.ax.plot(xs=[-1,1.5],ys=[0,0],color="lightgray", linestyle='--')
        self.ax.plot(xs=[0,0],ys=[-1,1.5],color="lightgray", linestyle='--')
        self.ax.plot(xs=[0,0],ys=[0,0],zs=[-1,1.5],color="lightgray", linestyle='--')
        # - axis labels
        self.ax.text(1.5,0,0,"x+")
        self.ax.text(0,1.5,0,"y+")
        self.ax.text(0,0,1.5,"z+")
        # - 'background'
        bgu = np.linspace(0, 2 * np.pi, 21)
        bgv = np.linspace(0, np.pi, 21)
        bgx = np.outer(np.cos(bgu), np.sin(bgv))
        bgy = np.outer(np.sin(bgu), np.sin(bgv))
        bgz = np.outer(np.ones(np.size(bgu)), np.cos(bgv))
        self.ax.plot_surface(bgx, bgy, bgz, color='b', alpha=0.05)

        #window setup
        self.fig.canvas.set_window_title(self.title)
        self.ax.set_aspect("equal")
        self.ax.set_xlim3d(-1, 1)
        self.ax.set_ylim3d(-1, 1)
        self.ax.set_zlim3d(-1, 1)
        self.ax.set_xticks([-1,0,1],minor=False)
        self.ax.set_yticks([-1,0,1],minor=False)
        self.ax.set_zticks([-1,0,1],minor=False)
    
    def add_pt_by_theta_phi(self,theta,phi,ptName=""):
        '''
        theta  - angle in radians from +z axis
        phi    - angle in radians from +x axis
        ptName - str name of this point (optional)

        returns a list with coords and alpha/beta:
        [ [a,b] , [x,y,z] ]
        '''
        # valid_ab = self.check_alpha_beta(theta,phi)
        # if(valid_ab == False):
        #     warningMessage = "This theta and phi look like an alpha and beta:\n"+str(theta)+", "+str(phi)+"\nDid you mean to add the point with alpha and beta?"
        #     warnings.warn(warningMessage)

        if(len(ptName)<1):
            ptName = str(theta)+", "+str(phi)
        endx = cmath.sin(theta) * cmath.cos(phi)
        endy = cmath.sin(theta) * cmath.sin(phi)
        endz = cmath.cos(theta)


        self.ax.plot(xs=[0,endx],ys=[0,endy],zs=[0,endz],
            label=ptName)
        self.ax.text(endx * self.label_dist_mult,
                     endy * self.label_dist_mult,
                     endz * self.label_dist_mult,
                     ptName)
        return [
            [ self.alpha_fr_theta(theta),
              self.beta_fr_theta_phi(theta,phi)],
            [endx,endy,endz]
        ]
    
    def add_pt_by_alpha_beta(self, alpha, beta, ptName=""):
        '''
        alpha  - probability of achieving |0>
        beta   - probability of achieving |1>
        ptName - str name of this point (optional)

        a^2 + b^2 = 1

        returns a list with coords and theta/phi:
        [ [theta,phi] , [x,y,z] ]
        '''
        # valid_ab = self.check_alpha_beta(alpha,beta)
        # if(valid_ab == False):
        #     warningMessage = "This alpha and beta don't match:\n"+str(alpha)+", "+str(beta)+"\nDid you mean to add the point with theta and phi?"
        #     warnings.warn(warningMessage)
        theta = self.theta_fr_alpha(alpha)
        phi   = self.phi_fr_beta_theta(beta,theta)
        pt    = self.add_pt_by_theta_phi(theta,phi,ptName)
        return [
            [ theta, phi ],
            pt[1]
        ]
    
    def show_sphere(self, showGrid=False, title=""):
        if(showGrid):
            self.ax.set_axis_on()
        else:
            self.ax.set_axis_off()
        if(len(title)>0):
            self.fig.canvas.set_window_title(title)
        
        plt.show()

        #cleanup?
        # self.fig.canvas.set_window_title(self.title)
    
    # def reset_sphere(self):

    
    #maths
    def theta_fr_alpha(self, alpha):
        return 2.0*cmath.acos(alpha)
    def alpha_fr_theta(self, theta):
        return cmath.cos(theta/2.0)
    def beta_fr_theta_phi(self, theta, phi):
        out = 0.0
        a = cmath.sin(theta/2.0)
        b = complex(0.0,theta)
        if(b==0):
            warningMessage = "Natural Log of zero attempted: "+str(a)+" * ln("+str(b)+"). Returning zero."
            warnings.warn(warningMessage)
        else:
            out = a * cmath.log(b)
        return out
    def phi_fr_beta_theta(self, beta, theta):
        out = 0
        b = beta/cmath.sin(theta/2.0)
        if(b==0):
            warningMessage = "Natural Log of zero attempted: -1 * ln("+str(b)+"). Returning zero."
            warnings.warn(warningMessage)
        else:
            out = (-1.0)*complex(0, cmath.log(b))
        return out
    def theta_fr_beta_phi(self, beta, phi):
        out = 0
        b = (-1.0)*complex(0.0, phi)
        if(b==0):
            warningMessage = "Natural Log of zero attempted: 2*arcsin( "+str(beta)+" * ln("+str(b)+") ). Returning zero."
            warnings.warn(warningMessage)
        else:
            out = 2.0*cmath.asin(beta*cmath.log(b))
        return out
    def beta_fr_alpha(self,alpha):
        return cmath.sqrt(1.0-math.pow(alpha, 2.0))
    def alpha_fr_beta(self,beta):
        return self.beta_fr_alpha(beta)
    def generate_alpha_beta_fr_theta(self,theta):
        return [cmath.sin(theta), cmath.cos(theta)]
    def check_alpha_beta(self,alpha,beta,tolerance=0.1):
        return (alpha**2)+(beta**2)-1 <= tolerance


if  __name__ == "__main__":
    b = BlochSphere()
    # b.add_pt_by_theta_phi(0.0,1.5)
    # b.add_pt_by_theta_phi(0.0,1.0)
    # b.add_pt_by_theta_phi(0.0,2.0)
    b.add_pt_by_alpha_beta(1.0/3.0,b.beta_fr_alpha(1.0/3.0),"oig")
    b.add_pt_by_theta_phi(4.3,2.4,"Hike")
    b.add_pt_by_theta_phi(1.4,3.7)
    b.add_pt_by_theta_phi(5.3,4.9,"Benji")
    b.show_sphere()


#sphere maths

#radius = sqrt(x^2 + y^2 + z^2)
#polar = arccos(z/radius)         aka theta
#azimuthal = atan2(y, x)          aka phi

#x = radius * sin(polar) * cos(azimuthal)
#y = radius * sin(polar) * sin(azimuthal)
#z = radius * cos(polar)


#latitude = polar - 90°
#longitude = azimuthal

#polar = latitude + 90°
#azimuthal = longitude

#radius = sqrt(x^2 + y^2 + z^2)
#latitude = arcsin(z/radius)
#longitude = atan2(y, x)

#x = radius * cos(latitude) * cos(longitude)
#y = radius * cos(latitude) * sin(longitude)
#z = radius * sin(latitude)

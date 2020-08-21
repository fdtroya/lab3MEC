from scipy.optimize import fsolve
import math
import numpy as np
import matplotlib.pyplot as plt



resolution=1000
step=2*math.pi/resolution
r1=0.212725
r2=0.0508
r3=0.1825498
r4=0.244475
theta4=2.3911
theta2_dot=52.3599

#variables
theta2=0
theta2s=[]


#incog
theta3=0
theta3_dot=0
theta3_doubledot=0
theta3_doubledots=[]

theta1=0    
theta1_dot=0    
theta1_doubledot=0
theta1_doubledots=[]

acg2x=[]
acg2y=[]

acg1x=[]
acg1y=[]

acg3x=[]
acg3y=[]



def displacement(p):
    theta3, theta1 = p
    return (r2*math.cos(theta2)-r4*math.cos(theta4)-r3*math.cos(theta3)+r1*math.cos(theta1), r2*math.sin(theta2)-r4*math.sin(theta4)-r3*math.sin(theta3)+r1*math.sin(theta1))

def velocity():
    sol=np.array([-r2*math.sin(theta2)*theta2_dot,r2*math.cos(theta2)*theta2_dot])
    equa=np.array([[-r3*math.sin(theta3),r1*math.sin(theta1)],[r3*math.cos(theta3),-r1*math.cos(theta1)]])
    X=np.linalg.inv(equa).dot(sol)
    return (X[0],X[1])

def acceleration():
    result1=-r2*math.cos(theta2)*theta1_dot**2 +r3*math.cos(theta3)*theta3_dot**2-r1*math.cos(theta1)*theta1_dot**2
    result2=-r2*math.sin(theta2)*theta1_dot**2 +r3*math.sin(theta3)*theta3_dot**2-r1*math.sin(theta1)*theta1_dot**2
    sol=np.array([result1,result2])
    equa=np.array([[-r3*math.sin(theta3),r1*math.sin(theta1)],[r3*math.cos(theta3),-r1*math.cos(theta1)]])
    X=np.linalg.inv(equa).dot(sol)
    return (X[0],X[1])
    


for i in range(resolution):
    theta2=theta2+step
    theta3, theta1 =  fsolve(displacement, (theta3,theta1))
    theta3_dot,theta1_dot=velocity()
    theta3_doubledot,theta1_doubledot=acceleration()
    theta2s.append(theta2)
    theta3_doubledots.append(theta3_doubledot)
    theta1_doubledots.append(theta1_doubledot)
    acg2x_o=-(r2/2)*(math.cos(theta2)*theta2_dot**2)
    acg2y_o=-(r2/2)*(math.sin(theta2)*theta2_dot**2)
    acg1x_o=2*acg2x_o-(r1/2)*(math.cos(theta1)*theta1_dot**2+theta1_doubledot*math.sin(theta1))
    acg1y_o=2*acg2y_o+(r1/2)*(-math.sin(theta1)*theta1_dot**2+theta1_doubledot*math.cos(theta1))
    acg3x_o=(-r3/2)*(math.cos(theta3)*theta3_dot**2+theta3_doubledot*math.sin(theta3))
    acg3y_o=(r3/2)*(-math.sin(theta3)*theta3_dot**2+theta3_doubledot*math.cos(theta3))

    acg2x.append(acg2x_o)
    acg2y.append(acg2y_o)
    acg1x.append(acg1x_o)
    acg1y.append(acg1y_o)
    acg3x.append(acg3x_o)
    acg3y.append(acg3y_o)
    



    

plt.xlabel("theta 2 [rad]")
plt.plot(theta2s,theta1_doubledots)
plt.ylabel("alpha 1 [rad/s^2]")
plt.show()
plt.plot(theta2s,theta3_doubledots)
plt.ylabel("alpha 3 [rad/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg2x)
plt.ylabel("acg2x [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg2y)
plt.ylabel("acg2y [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg1x)
plt.ylabel("acg1x [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg1y)
plt.ylabel("acg1y [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg3x)
plt.ylabel("acg3x [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()
plt.plot(theta2s,acg3y)
plt.ylabel("acg3y [m/s^2]")
plt.xlabel("theta 2 [rad]")
plt.show()






















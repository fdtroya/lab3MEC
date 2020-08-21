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
    

plt.xlabel("theta 2 [rad]")
plt.plot(theta2s,theta1_doubledots)
plt.ylabel("alpha 1 [rad/s**2]")
plt.show()
plt.plot(theta2s,theta3_doubledots)
plt.ylabel("alpha 3 [rad/s**2]")
plt.show()
























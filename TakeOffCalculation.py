import numpy as np
g = 9.81
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
g = 9.81

def MaximalTakeOffWeight(crew_weight, payload_weight, fuel_weight_fraction, empty_weight):
    return (crew_weight + payload_weight + empty_weight)/(1 - fuel_weight_fraction)

def BreguetRange(propeller_efficiency, lift_by_drag, start_weight, end_weight, power_sfc):
    return (propeller_efficiency/(power_sfc*g))*lift_by_drag*np.log(start_weight/end_weight)

def _A_TakeOffDistance(thrust, weight, friction_coeff):
    return g*((thrust/weight) - friction_coeff)

def _B_TakeOffDistance(density, weight, wing_area, drag_coeff, lift_coeff, friction_coeff, a=0):
    return (g/weight)*(0.5*density*wing_area*(drag_coeff - friction_coeff*lift_coeff) + a)

def TakeOffDistance(density, weight, wing_area, drag_coeff, lift_coeff, friction_coeff, thrust, velocity, a=0):
    A = _A_TakeOffDistance(thrust, weight, friction_coeff)
    B = _B_TakeOffDistance(density, weight, wing_area, drag_coeff, lift_coeff, friction_coeff)

    return (1/(2*B))*np.log((A)/(A - B*velocity*velocity))

def LoiterEndurance(lift_by_drag, power_sfc, velocity, propeller_efficiency, start_weight, end_weight):
    return ((lift_by_drag*propeller_efficiency)/(power_sfc*velocity))*np.log(start_weight/end_weight)

def TakeOffWeightFraction(fuel_consumption_rate=None): #To be changed according to fuel consumption rate of the additional thruster
    return 1

def ClimbWeightFraction(fuel_consumption_rate=None):
    return 1

def CruiseWeightFraction(breguet_range, propeller_efficiency, lift_by_drag, power_sfc):
    return np.exp(-(breguet_range*(power_sfc*g))/(lift_by_drag*propeller_efficiency))

def LoiterWeightFraction(lift_by_drag, power_sfc, propeller_efficiency, endurance, velocity):
    return np.exp(-(endurance*power_sfc*velocity*g)/(propeller_efficiency*lift_by_drag))

def LandingWeightFraction():
    return 0.9965

def DragCoefficient(min_drag_coeff, lift_coeff, lift_coeff_drag_min, aspect_ratio, e):
    return min_drag_coeff + (1/(np.pi*e*aspect_ratio))*((lift_coeff - lift_coeff_drag_min)**2)
def C_D(C_D_min, AR, e, C_L_D_min, C_L):
    k = 1/(np.pi*e*AR)
    return C_D_min + k*(C_L - C_L_D_min)**2

def Thrust(thruster, power, velocity):
    return thruster + 2*(power/velocity)

def A(thrust, weight, mu, thruster, theta):
    return g*((thrust/weight) + (thruster/weight)*np.cos(theta) - mu*(1 - (thruster/weight)*np.cos(theta)))

def B(weight, rho, S, dragCoeff, C_L, mu):
    return (g/weight)*(0.5*rho*S*(dragCoeff - mu*C_L))

def v_takeoff(cl_max,w,thrust,theta,rho,S):
    net_weight=w-thrust*np.cos(theta)
    return np.sqrt(2*net_weight/(cl_max*rho*S))

def nozzle_angle_delta_cl(angle):
    deg=angle*180/np.pi
    return 0.000061308596165739*deg-0.00103

def nozzle_angle_delta_cd(angle):
    deg=angle*180/np.pi
    return 0.0000523586889301176*deg-0.00495

def takeoffDistance(weight, rho, S, C_L_max, mu, AR, C_D_min, e, C_L_D_min, powerPerEngine, takeOffVelocity=None,thruster=0, theta=0):
    dragCoeff = C_D(C_D_min, AR, e, C_L_D_min, C_L_max)
    takeOffVelocity=v_takeoff(C_L_max,weight,thruster,theta,rho,S)
    #print(takeOffVelocity)
    totalThrust = Thrust(thruster, powerPerEngine, takeOffVelocity)
    a = A(totalThrust, weight, mu, thruster, theta)
    b = B(weight, rho, S, dragCoeff, C_L_max, mu)
    return (1/(2*b))*(np.log(a/(a-b*takeOffVelocity*takeOffVelocity)))

def weight_to_thrust(weight):
    return (0.051532*weight+0.8666)*1000

if __name__ == "__main__":
    ideal_range = 1000000 # meters
    propeller_efficiency = 0.63
    lift_by_drag_max = 16 # Used During Cruise
    power_sfc = 0.325/(10*10*10*60*60) # kg/Ws
    loiter_time = 15*60 # s
    loiter_speed = 56.59 # m/s
    lift_by_drag_loiter_max = 0.86*lift_by_drag_max

    cruiseW = CruiseWeightFraction(ideal_range, propeller_efficiency, lift_by_drag_max, power_sfc)
    loiterW = LoiterWeightFraction(lift_by_drag_loiter_max, power_sfc, propeller_efficiency, loiter_time, loiter_speed)
    
    fuel_fraction = 1 - TakeOffWeightFraction()*ClimbWeightFraction()*loiterW*cruiseW*LandingWeightFraction()
    print(fuel_fraction)
    efficiency=0.75
    thrusterWeight = 61.2
    angles=np.linspace(0,np.pi/2,90)
    
    engine_weights=np.linspace(10,150,200)
    print(engine_weights.shape)
    vals=np.zeros((angles.shape[0],engine_weights.shape[0]))
    cnt_i=0
    cnt_j=0
    for i in angles:
        cnt_j=0
        for j in engine_weights:
            vals[cnt_i][cnt_j]=takeoffDistance((MaximalTakeOffWeight(170,1620,fuel_fraction,3900+j))*g,0.76*1.225,32,1.7+0.21+nozzle_angle_delta_cl(i),0.22,9,0.029+nozzle_angle_delta_cd(i),0.63,0.17,efficiency*581646,49.33, weight_to_thrust(j),i)
            cnt_j=cnt_j+1
        cnt_i+=1
            #vals.append(takeoffDistance((MaximalTakeOffWeight(170,1620,fuel_fraction,3900+j))*g,0.76*1.225,32,1.7+0.21+nozzle_angle_delta_cl(i),0.22,9,0.029+nozzle_angle_delta_cd(i),0.63,0.17,efficiency*581646,49.33, weight_to_thrust(j),i))
    #plt.plot(angles,vals)
    #plt.show()
    fig=plt.figure(figsize=(20,20))
    ax=fig.add_subplot(111,projection='3d')
    x,y=np.meshgrid(angles,engine_weights)

    print(np.unravel_index(np.argmin(vals), vals.shape))

    ax.plot_surface(x,y,vals.T, cmap=cm.jet)
    ax.set_xlabel("Angles (Radian)")
    ax.set_ylabel("Engine Weights (KG)")
    ax.set_zlabel("Take Off Distance")
    np.ma
    print(vals[:, 100])
    plt.show()
    print(vals)
    print(vals.shape)
    min_w=500
    weights=[]
    config_angles=[]
    angle=0
    min_takeoff_config=0
    for i in range(90):
        for j in range(200):
            if(vals[i][j]<375 and angles[i]>np.pi/4):
                weights.append(engine_weights[j])
                config_angles.append(angles[i])
                if(engine_weights[j]<min_w):
                    min_w=engine_weights[j]
                    angle=angles[i]
                    min_takeoff_config=vals[i][j]
    print(min(weights))
    print(weight_to_thrust(min(weights)))
    print(takeoffDistance((MaximalTakeOffWeight(170,1620,fuel_fraction,3900+61.2))*g,0.76*1.225,32,1.7+0.21+nozzle_angle_delta_cl(np.pi/4),0.22,9,0.029+nozzle_angle_delta_cd(np.pi/4),0.63,0.17,efficiency*581646,49.33, 5700,np.pi/4))
    #print(takeoffDistance((MaximalTakeOffWeight(170,1620,fuel_fraction,3900+73))*g,0.76*1.225,32,1.7+0.21+nozzle_angle_delta_cl(np.pi/3),0.22,9,0.029+nozzle_angle_delta_cd(np.pi/3),0.63,0.17,efficiency*581646,49.33, 3250,np.pi/3))
    print(min_w,angle,min_takeoff_config)

    #print(takeoffDistance((6400 + thrusterWeight)*g,0.76*1.225,32,1.7,0.22,9,0.029,0.63,0.17,efficiency*581646,49.33, 12680, 0))
    # 440 5700; 188 16000;
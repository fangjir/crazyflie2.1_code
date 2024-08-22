#@author:fangjie
#@verson:
import asyncio
import rospy
import math
import time
import csv
from tkinter import N
import xml.etree.cElementTree as ET
from threading import Thread
from matplotlib import markers
from threading import Timer


from geometry_msgs.msg import PoseStamped
from geometry_msgs.msg import TwistStamped
from geometry_msgs.msg import AccelStamped


import cflib.crtp
from cflib.crazyflie import Crazyflie
from cflib.crazyflie.log import LogConfig
from cflib.crazyflie.mem import MemoryElement

from cflib.utils import uri_helper

from cflib.crtp.crtpstack import CRTPPacket
from controller import PID
from controller import PID1
from controller import ESO,ESO1,FTDO,dobc

from numpy        import arcsin, rad2deg
from numpy        import sqrt,arctan
from numpy        import cos, sin, deg2rad

from numpy.linalg import norm


# Guidance
import numpy as np
import numpy as np
import matplotlib.pyplot as plt

from minsnap_traj import minimum_snap_traj,minimum_snap_traj_p2p,get_traj


uri = uri_helper.uri_from_env(default='radio://0/100/2M/E7E7E7E701')

# The name of the rigid body in QTM that represents the Crazyflie
rigid_body_name = 'cf1'

# True: send position and orientation; False: send position only
send_full_pose = False


orientation_std_dev = 8.0e-3


global pos_cb , vel_cb, acc_cb, euler_cb, angvel_cb
global const_thrust, h_int, thrusts_, acc_real,z_list,pos_ref,desired_vel,dt,posex,posey
pos_cb = [0,0,0]
pos_ref = [0,0,0]
vel = [0,0,0]
euler_cb = [0,0,0]
angvel_cb = [0,0,0]
vel_ref = [0,0,0]
posex = 0
posey = 0

#RLS
acc_real = 0.0
rho = 1.0
# p = 1000000.0
gamma = 0
acc2thr_ = 3846.0
thrust_current = 0


z_list = []
x_list =[]
y_list= []
d_list = []
q_list = []
r_list = []
d_hat1 = 0
d_hat2 = 0
d_hat21 =0
d_hat3 = 0
x1_hat = 0
x2_hat = 0

x1 = 0
x2 = 0
x3 = 0
dt = 0.02
t = 0
i=0
acc1=0
acc2=0
acc3=0

run_enabled = 1


thrust_limit = (10001, 60000)
roll_limit   = (-10.0, 10.0)
pitch_limit  = (-10.0, 10.0)
# controller1 = PID1(9,6,3,10,dt)# 6 5 7
# controller2 = PID1(10,6.2,3.5,10,dt)# 

# controller3 = PID(14,7,5.4,dt) # 6 7 3.5 6 7 4
controller1 = PID1(4.5,0,5,10,dt)# k2,k1 2.1 1.0,2.1 2.7    2.44 2.4    2,5
controller2 = PID1(4.5,0,5,10,dt)  # 4.5 5 current use

controller3 = PID(5,0,3,dt) #one 6 7 3.5  second 6 7 4  third  14 7 5.4


#ESO 
omega1 = 5
omega2 = 4 #x
omega3 = 4 #y
d = 0
acc = [0,0,9.81]
a1 = [0,0,0]
# a = 0
command = [0,0,0,0]


class QtmWrapper(Thread):
    def __init__(self, body_name):
        Thread.__init__(self)
        self.sub1 = rospy.Subscriber("/vrpn_client_node/MCServer/3/pose",PoseStamped,self.pos_callback)
        self.sub2 = rospy.Subscriber("/vrpn_client_node/MCServer/3/velocity",TwistStamped,self.vel_callback)
        #self.sub3 = rospy.Subscriber("/vrpn_client_node/MCServer/3/acceleration",AccelStamped,self.acc_callback)
        self.on_pose = None
        self.first = 1
        self.i = 1
        self.j = 1
        self.z = 1
        self.pos_init = [0,0,0]
        self.pos_mocap = [0,0,0]
        self.vel_mocap = [0,0,0]
        self.acc_cur = [0,0,0]
        self.phi = 0
        self.theta = 0
        self.pesi = 0
        self.pesi_d = 0
        self.start()

    def close(self):
        self.sub1.unregister()
        self.sub2.unregister()

    def pos_callback(self,msg):
            if self.first == 1:
                self.pos_init = np.array([(msg.pose.position.y),-1*(msg.pose.position.x),(msg.pose.position.z)])
                #self.pos_init[2] = 0.3
                self.pesi_d = (msg.pose.orientation.z)
                self.first = 0
            self.pos_mocap = np.array([(msg.pose.position.y-self.pos_init[0])/1000,(-msg.pose.position.x-self.pos_init[1])/1000,(msg.pose.position.z-self.pos_init[2])/1000])
            self.pesi = (msg.pose.orientation.z)
            self.phi = (msg.pose.orientation.x)
            self.theta = (msg.pose.orientation.y) 
                #if self.on_pose:
                    #self.on_pose([msg.pose.position.x/1000,msg.pose.position.y/1000,msg.pose.position.z/1000,msg.pose.orientation ])

    
    def vel_callback(self,msg):
        
            self.vel_mocap = np.array([(msg.twist.linear.y/1000),(-msg.twist.linear.x/1000),(msg.twist.linear.z/1000)])
            
    def acc_callback(self,msg):
        if self.z==1:
            self.acc_cur = np.array([msg.accel.linear.x/1000,msg.accel.linear.y/1000,msg.accel.linear.z/1000 ])
            self.z=0
        else:
            self.z = self.z + 1



def reset_estimator(cf):
    cf.param.set_value('kalman.resetEstimation', '1')
    time.sleep(0.1)
    cf.param.set_value('kalman.resetEstimation', '0')
    cf.param.set_value('mypara.controller_switch', '0')
    time.sleep(0.1)
    # time.sleep(1)
    #wait_for_position_estimator(cf)



def stop_logconf():
    log_conf = LogConfig(name='Position', period_in_ms=100)
    log_conf.stop()

    log_conf = LogConfig(name='Stabilizer', period_in_ms=100)
    log_conf.stop()



def adjust_orientation_sensitivity(cf):
    cf.param.set_value('locSrv.extQuatStdDev', orientation_std_dev)


def activate_kalman_estimator(cf):
    cf.param.set_value('stabilizer.estimator', '2')

    # Set the std deviation for the quaternion data pushed into the
    # kalman filter. The default value seems to be a bit too low.
    cf.param.set_value('locSrv.extQuatStdDev', 0.06)

def position_callback_euler(timestamp, data, logconf):
    global euler_cb
    # euler_cb[0] = data['stateEstimate.roll']
    # euler_cb[1] = data['stateEstimate.pitch']
    euler_cb[2] = data['stateEstimate.yaw']

def start_position_printing_euler(cf):
    log_conf = LogConfig(name='Stabilizer', period_in_ms=10)
    # log_conf.add_variable('stateEstimate.roll', 'FP16')
    # log_conf.add_variable('stateEstimate.pitch', 'FP16')
    log_conf.add_variable('stateEstimate.yaw', 'FP16')
    # log_conf.add_variable('acc.x', 'FP16')
    # log_conf.add_variable('acc.y', 'FP16')
    # log_conf.add_variable('acc.z', 'FP16')
    cf.log.add_config(log_conf)
    log_conf.data_received_cb.add_callback(position_callback_euler)
    log_conf.start()
    

def stop(cf):
    cf.close_link()
    print("start to close link")

def run(event):
    global run_enabled
    if run_enabled == 1:
        cf.commander.send_setpoint(0.0,0.0,0.0,0)
        run_enabled = 0
    calc_control_input(cf,qtm_wrapper)



def calc_control_input(cf,qt):  #have 17*0.02 time delay
    global t,x,acc,command,d_hat1,d_hat2,d_hat3,d_hat21,x1 ,x2,a1,a,posex,posey,x1,x2,x3,euler_cb,acc1,acc2,i,p,a,acc3
    time_finish = time.time() 
    # z = qt.pos_mocap[2]
    # x = qt.pos_mocap[0]
    # y = qt.pos_mocap[1]
    # pos_ref[0]=p[0][i]
    # # print(pos_ref[0])
    # pos_ref[1]=p[1][i]
    # pos_ref[2]=p[2][i]
    # acc1=a[0][i]
    # acc2=a[1][i]
    # acc3=a[2][i]
    # i=i+1
    y_list.append(qt.pos_mocap[1]-pos_ref[1]) 
    r_list.append(qt.pos_mocap[0]-pos_ref[0])
    x_list.append(qt.pos_mocap[2]-pos_ref[2])#qt.pos_mocap[2]-pos_ref[2]
    if time_finish - time_start >= 5.0:
        pos_ref[0]=p[0][i]
        # print(pos_ref[0])
        pos_ref[1]=p[1][i]
        pos_ref[2]=p[2][i]
        acc1=a[0][i]
        acc2=a[1][i]
        acc3=a[2][i]
        i=i+1
    #     t = time_finish - time_start - 6.0
        # pos_ref[0] = pos_ref[0] + 0.002
        # pos_ref[1] = pos_ref[1] + 0.001
    #     pos_ref[0]=-0.3+0.3*math.cos(0.5*t)
    #     pos_ref[1]= 0.3*math.sin(0.5*t)
    #     acc1=-0.075*math.cos(0.5*t)
    #     acc2=-0.075*math.sin(0.5*t)
    #     # if time_finish-time_start >4.0  and time_finish-time_start <14:
    #     #   pos_ref[2]=pos_ref[2]+0.001
    #     # if time_finish-time_start>14:
    # #     #     pos_ref[2]=pos_ref[2]-0.001
            
    cosyaw = cos(qt.pesi_d)
    sinyaw = sin(qt.pesi_d)
    set_body_x = pos_ref[0]*cosyaw+pos_ref[1]*sinyaw
    set_body_y = -pos_ref[0]*sinyaw+pos_ref[1]*cosyaw
    set_acc1=acc1*cosyaw+acc2*sinyaw
    set_acc2=-acc1*sinyaw+acc2*cosyaw

    state_body_x = qt.pos_mocap[0]*cosyaw +qt.pos_mocap[1]*sinyaw
    state_body_y = -qt.pos_mocap[0]*sinyaw+qt.pos_mocap[1]*cosyaw
 #################   
    #update control value
    #desired_acc,deisred_vel,act_vel,desired_pos,pos

    acc[0] = (controller1.update(set_acc1,0,0,set_body_x,state_body_x))#x,pitch
    acc[1] = (controller2.update(set_acc2,0,0,set_body_y,state_body_y))#y,roll
    acc[2] =  controller3.update(acc3,0,0,pos_ref[2],qt.pos_mocap[2])+9.81
    
    #######
    # d_hat3 = estimate1.update(qt.vel_mocap[2],acc[2]-9.81-d_hat3)  #z
    d_hat1 = estimate2.update2(state_body_x,acc[0]-d_hat1)
    d_hat2 = estimate3.update2(state_body_y,acc[1]-d_hat2)
    d_hat3 = estimate1.update2(qt.pos_mocap[2],acc[2]-9.81-d_hat3)

    x1 = estimate2.get_x1()
    x2 = estimate3.get_x1()
    d_list.append(pos_ref[0])
    #d_hat1->z
    #d_hat2->x
    #d_hat3->y
################################ 
    a1[0] = acc[0]-d_hat1
    a1[1] = acc[1]-d_hat2
    a1[2] = acc[2]-d_hat3
    a2 = norm(a1)
    yaw = deg2rad(euler_cb[2])
    q_list.append(pos_ref[1])
    pitch_d = arctan((a1[0]*cos(yaw)+a1[1]*sin(yaw))/a2)
    roll_d = arctan(cos(pitch_d)*(a1[0]*sin(yaw)-a1[1]*cos(yaw))/a2)
    # pitch_d = arctan((a1[0]*cos(0)+a1[1]*sin(0))/a1[2])
    # roll_d = arctan(cos(pitch_d)*(a1[0]*sin(0)-a1[1]*cos(0))/a1[2])

    roll_d  = rad2deg(roll_d)
    pitch_d = rad2deg(pitch_d)
    
    thrust_r = 3004*(a2)
###############################

    thrust_r = np.clip(thrust_r,*thrust_limit) #thrust
    if qt.pos_mocap[2] >0.8 :
       cf.commander.send_stop_setpoint() 
       timer1.shutdown()
    cf.commander.send_setpoint(roll_d,pitch_d,0,int(thrust_r)) #-pitch,,,,x-> pitch,y->roll
    #cf.commander.send_setpoint(0,0,0,int(10002))
                            #roll,pitch  
    ###############
    #record
    z_list.append(pos_ref[2])

    # print(i)

if __name__ == '__main__':
    cflib.crtp.init_drivers()
    rospy.init_node('mocap',anonymous=True)
    qtm_wrapper = QtmWrapper(rigid_body_name)

    cf =Crazyflie(rw_cache='./cache')
    cf.open_link(uri)

    adjust_orientation_sensitivity(cf)
    activate_kalman_estimator(cf)
    reset_estimator(cf)
    start_position_printing_euler(cf)

    pos_ref[2] = 0.4
    posex = pos_ref[0]
    posey = pos_ref[1]

    # estimate1 = ESO(0,3*omega1,3*omega1*omega1,omega1*omega1*omega1,dt) #z
    # estimate2 = ESO(0,3*omega2,3*omega2*omega2,omega2*omega2*omega2,dt) #x
    # estimate3 = ESO(0,3*omega3,3*omega3*omega3,omega3*omega3*omega3,dt) #y
    estimate1 = ESO1(0,0.01,0.00001,100,dt) #z
    estimate2 = ESO1(0,0.0001,0.00001,1000,dt) #x
    estimate3 = ESO1(0,0.001,0.00001,1000,dt) #y  0.001  0.1 10000

    print("start to print velocity")
    
    way_points =np.array([[0,0,0.4,0],[0.3,0.3,0.55,0],[0.5,0.5,0.65,0]])
    print("trajectory point:",way_points)
    time_set = np.array([0,5,10.02])
    n_order = 5
    n_obj = 4
    v_i = [0,0,0,0]
    a_i = [0,0,0,0]
    v_e = [0,0,0,0]
    a_e = [0,0,0,0]
    sample_rate = 50
    Matrix_x, Matrix_y, Matrix_z = minimum_snap_traj_p2p(way_points, time_set, n_order, n_obj, v_i, a_i, v_e, a_e)
    p, a,t_list= get_traj(Matrix_x, Matrix_y, Matrix_z, time_set, sample_rate)
    # print(p[0])
    time_start = time.time()
    timer1 = rospy.Timer(rospy.Duration(dt), run)
    time.sleep(14.98)
    #cf.commander.send_setpoint(0.0,0.0,0.0,0)
    timer1.shutdown()

    force = 35000
    while force>0:
        cf.commander.send_setpoint(0, 0, 0, force)
        time.sleep(0.1)
        force-=500
        if force<30000:
            break
    cf.commander.send_setpoint(0.0,0.0,0.0,0)
    stop(cf)
    time.sleep(1)
    qtm_wrapper.close()
    ### plots
    
    with open('errory.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(y_list)
    
    with open('errorx.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(r_list)

    with open('yref.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(q_list)
    with open('xref.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(d_list)
    with open('errorz.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(x_list)
    with open('zref.csv','w',newline='') as myflie:
        wr = csv.writer(myflie,quoting=csv.QUOTE_ALL)
        wr.writerow(z_list)
    fig = plt.figure()

    ax1 = fig.add_subplot(611)# z axis
    ax1.plot(z_list) 
    ax1.set_ylabel('zref')
    plt.grid()

    ax2 = fig.add_subplot(612)  #x estimate error
    ax2.plot(x_list) 
    ax2.set_ylabel('errorz')
    plt.grid()

    ax3 = fig.add_subplot(613)  #  error of y
    ax3.plot(y_list) 
    ax3.set_ylabel('e_of_y')
    plt.grid()

    ax4 = fig.add_subplot(614)  # error of x
    ax4.plot(r_list) 
    ax4.set_ylabel('e_of_x')
    plt.grid()

    ax5 = fig.add_subplot(615)  #  y estimate  error
    ax5.plot(d_list) 
    ax5.set_ylabel('pos_ref[0]')
    plt.grid()

    ax6 = fig.add_subplot(616)  #  yaw from body
    ax6.plot(q_list) 
    ax6.set_ylabel('pos_ref[1]')
    plt.grid()

    plt.show() 

    #estimate1.plot()

    # fig = plt.figure()
    # plt.scatter(t_list[1:], p[0], marker = 'x', color = 'blue', s = 2, label = '0')
    # plt.show()

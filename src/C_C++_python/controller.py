'''PID class'''
# first simple PID
import math
import numpy as np
import matplotlib.pyplot as plt


class PID:
        def __init__(self,kp,ki,kd,dt):
            self.kp = kp 
            self.ki = ki 
            self.kd = kd 
            self.kf = 1
            self.error = 0
            self.pre_error = 0
            self.dt = dt 
            self.outp = 0
            self.outi = 0
            self.deriv = 0
            self.out  = 0
            self.integ = 0
            self.ff = 0
            self.flag1 = 1
            self.flag2 = 1
            self.fc = 10.0
            self.Ts = dt
            self.pi = 3.14159
            self.rc = 1.0/(2.0* self.pi * self.fc)
            self.alpha = self.Ts/(self.Ts + self.rc)
            self.last_deriv = 0
            self.g = 9.81

            self.alpha1 = 5/7
            self.alpha2 = 2* self.alpha1/(1+self.alpha1)
        def update(self,desired_acc,desired_vel,act_vel,desired,measure):
            self.out = 0
            self.error = desired - measure 
            #self.outp = self.kp * self.error 
            #self.outp = self.kp*math.pow(abs(self.error),self.alpha1)*self.flag1
            self.outp = self.kp*self.error
            self.integ = self.integ+self.dt*self.error
            self.outi = self.ki *self.integ
            
            self.deriv = self.kd*((self.error - self.pre_error)/self.dt)  #deriv           
            self.deriv = self.last_deriv + (self.deriv - self.last_deriv) * self.alpha
            self.last_deriv = self.deriv
            self.pre_error = self.error
            self.out = desired_acc+(self.outp + self.outi + self.deriv)
            return self.out
        
        def get_p(self):
            return self.outp

        def get_i(self):
             return self.outi
        
        def get_d(self):
             return self.deriv

class PID1:
        def __init__(self,kp,ki,kd,fc,dt):
            self.kp = kp 
            self.ki = ki 
            self.kd = kd 
            self.error = 0
            self.pre_error = 0
            self.dt = dt 
            self.outp = 0
            self.outi = 0
            self.deriv = 0
            self.out  = 0
            self.integ = 0
            self.fc = fc
            self.Ts = dt
            self.pi = 3.14159
            self.rc = 1.0/(2.0* self.pi * self.fc)
            self.alpha = self.Ts/(self.Ts + self.rc)
            self.last_deriv = 0
            self.g = 9.81

            self.alpha1 = 5/7
            self.alpha2 = 2* self.alpha1/(1+self.alpha1)
        def update(self,desired_acc,desired_vel,act_vel,desired,measure):
            self.out = 0
            self.error = desired - measure 
            #self.outp = self.kp * self.error 
            #self.outp = self.kp*math.pow(abs(self.error),self.alpha1)*self.flag1
            self.integ = self.integ+self.dt*self.error
            self.outp = self.kp*self.error
            self.outi = self.ki *self.integ
            
            self.deriv = self.kd*((self.error - self.pre_error)/self.dt)  #deriv          
            #self.deriv = self.last_deriv + (self.deriv - self.last_deriv) * self.alpha
            self.last_deriv = self.deriv
            self.pre_error = self.error
            self.out = desired_acc+(self.outp + self.outi + self.deriv)
            return self.out
        
        def get_p(self):
            return self.outp

        def get_i(self):
             return self.outi
        
        def get_d(self):
             return self.deriv

class ESO:
     def __init__(self,x1,peta1,peta2,peta3,dt):
          self.peta1 = peta1
          self.peta2 = peta2
          self.peta3 = peta3
          self.dt    = dt
          self.x1    = x1
          self.x2    = 0
          self.x3    = 0 
          self.z1    = 0
          self.z2    = 0
          self.kp    = 1
          self.d     = 0
          self.ki    = 0.1
          self.inte  = 0
        

     def  update(self,y,u):
          e = self.x1 - y
          self.x1 = self.x1 + self.dt*(self.x2 - self.peta1 * e)
          self.x2 = self.x2 + self.dt*(self.x3 + u - self.peta2 *e)
          self.x3 = self.x3 + self.dt*(-self.peta3*e)
          #self.x2 = self.x2 + self.dt*(self.x3 + u - self.peta2 *e)
          #self.x1 = self.x1 + self.dt*(self.x2 - self.peta1 * e)
          return self.x3
     def update2(self,y,u):
          self.z2=self.z2+self.dt*(u+self.d)
          self.z1=self.z1+self.dt*self.z2
        #   self.z2=self.z2+self.dt*(u+self.d)
          e=y-self.z1
          self.inte = self.inte+self.dt*e
          self.d = self.kp*e+self.ki*self.inte
          return self.d
     
     def get_x1(self):
          return self.x1
     
     def get_x2(self):
          return self.x2

class ESO1:
     def __init__(self,x1,q,r,p,dt):
        #   self.peta1 = peta1
        #   self.peta2 = peta2
        #   self.peta3 = peta3
          self.peta1 = 15
          self.peta2 = 75
          self.peta3 = 125
          self.dt    = dt
          self.x1    = x1
          self.x2    = 0
          self.x3    = 0 
          self.x4    = 0
          self.Q=np.array([[q,0,0],
                         [0,q,0],
                         [0,0,q]])
          self.R=r
          self.P=np.array([[p,0,0],
                         [0,p,0],
                         [0,0,p]])
          self.A=np.array([[1,self.dt,0],
                         [0,1,self.dt],
                         [0,0,1]])
          self.B=np.array([[0],[self.dt],[0]])
          self.H=np.array([[1,0,0]])
          self.x_pri=np.array([[0],
                        [0],
                        [0]])
          self.x_post=np.array([[0],
                        [0],
                        [0]])
          self.x1_hat=[]
          self.x2_hat=[]
          self.x3_hat=[]

     def  update(self,y,u):
          e = self.x1 - y
          #self.x1 = self.x1 + self.dt*(self.x2 - self.peta1 * e)
          #self.x2 = self.x2 + self.dt*(self.x3 + u - self.peta2 *e)
          self.x3 = self.x3 + self.dt*(-self.peta3*e)
          self.x2 = self.x2 + self.dt*(self.x3 + u - self.peta2 *e)
          self.x1 = self.x1 + self.dt*(self.x2 - self.peta1 * e)
          return self.x3
     
     def update2(self,y,u):
          self.x_pri=np.dot(self.A,self.x_post)+self.B*u  ###update 1
          self.P=np.dot(np.dot(self.A,self.P),(self.A.T))+self.Q  ## update 2

          k1=np.dot(self.P,(self.H.T))
          k2=np.dot(np.dot(self.H,self.P),self.H.T)+self.R
          k=np.dot(k1,1/k2)

          kalman_e=y-self.x_pri[0,0]
          self.x_post = self.x_pri +np.dot(k,kalman_e)
          self.P=np.dot((np.eye(3)-k*self.H),self.P)
        #   self.x1_hat.append(self.x_post[0,0])
        #   self.x2_hat.append(self.x_post[1,0])
        #   self.x3_hat.append(self.x_post[2,0])
          return self.x_post[2,0]
     def plot(self):
        fig = plt.figure()
        ax2 = fig.add_subplot(311)  #x estimate error
        ax2.plot(self.x1_hat) 
        ax2.set_ylabel('x1_hat')
        plt.grid()

        ax4 = fig.add_subplot(312)  # error of x
        ax4.plot(self.x2_hat) 
        ax4.set_ylabel('x2_hat')
        plt.grid()

        ax5 = fig.add_subplot(313)  #  y estimate  error
        ax5.plot(self.x3_hat) 
        ax5.set_ylabel('x3_hat')
        plt.grid()

        plt.show() 

     def get_x1(self):
          return self.x_post[0,0]
     
     def get_x2(self):
          return self.x_post[1,0]
    
class dobc:
     def __init__(self,l,dt):
          self.z=0
          self.z_pre=0
          self.l = l
          self.dt = dt
     def update(self,x2,u):
          self.z=self.z_pre + self.dt*(-self.l*self.z-self.l*(self.l*x2+u))
          self.z_pre=self.z

          return self.z+self.l*x2

class feedback:
        def __init__(self,kp,kd,dt):
            self.kp = kp #k1k2 + 1
            self.kd = kd #k1 + k2
            self.error = 0
            self.error_dot = 0
            self.dt = dt 
            self.out  = 0
            # self.ff = 0
            # self.fc = 10.0
            # self.Ts = dt
            # self.pi = 3.14159
            # self.rc = 1.0/(2.0* self.pi * self.fc)
            # self.alpha = self.Ts/(self.Ts + self.rc)

            # self.alpha1 = 5/7
            # self.alpha2 = 2* self.alpha1/(1+self.alpha1)


        def update(self,desired_acc,desired_vel,act_vel,desired,measure):
            self.out = 0
            self.error = desired - measure 
            self.error_dot = desired_vel - act_vel
            self.out = self.kp * self.error + self.kd * self.error_dot + desired_acc
            return self.out

class FTDO():
    def __init__(self,x1,lamda1,lamda2,lamda3,L,dt):
        self.x1 = 0 # 速度
        self.x2 = 0  # 扰动d
        self.x3 = 0  # 一阶扰动
        self.lamda1 = lamda1
        self.lamda2 = lamda2
        self.lamda3 = lamda3
        self.L = L
        self.dt = dt
    
    def sig(self,e):
        if e>0:
            return  1
        elif (e==0):
             return 0
        else:
            return -1
        

    def update(self,y,u):
        e1 = self.x1 - y
        e1 = self.x2 - self.lamda1*math.pow(self.L,1/3)*math.pow(abs(e1),2/3)*self.sig(e1)
        e2 = self.x2 - e1
        e2 = self.x3 - self.lamda2*math.pow(self.L,1/2)*math.pow(abs(e2),1/2)*self.sig(e2)
        e3 = self.x3 - e2
        e3 = -self.lamda3*self.L*self.sig(e3)

        #FTDO
        self.x1 = self.x1 + self.dt*(u+e1)
        self.x2 = self.x2 + self.dt*(e2)
        self.x3 = self.x3 + self.dt*(e3)

        return self.x2
    def get_x1(self):
          return self.x1
     
    def get_x2(self):
          return self.x2

class polynomial:
        def __init__(self,p):
            self.p = p
        def eval(self,t):
            assert t >= 0
            x = 0.0
            for i in range(0, len(self.p)):
                x = x * t + self.p[len(self.p) - 1 - i]
            return x

    # compute and return derivative
        def derivative(self):
            return polynomial([(i+1) * self.p[i+1] for i in range(0, len(self.p) - 1)])
        
class trajectoryoutput:
        def __init__(self):
            self.pos = []
            self.vel = []
            #self.yaw = None

class polynominal4D:
        def __init__(self,duration,px,py,pz,pyaw):
            self.duration = duration
            self.px = polynomial(px)
            self.py = polynomial(py)
            self.pz = polynomial(pz)
            self.pyaw = polynomial(pyaw)
        
        def derivative(self):
            return polynominal4D(
                self.duration,
                self.px.derivative().p,
                self.py.derivative().p,
                self.pz.derivative().p,
                self.pyaw.derivative().p
            )

        def eval(self,t):
            result = trajectoryoutput()
            result.pos = np.array([self.px.eval(t),self.py.eval(t),self.pz.eval(t)])
            derivative = self.derivative()
            result.vel = np.array([derivative.px.eval(t),derivative.py.eval(t),derivative.pz.eval(t)])

            return result

class trajectory():
        def __init__(self):
            self.polynomials = []
            self.duration = 0.0

        def n_pieces(self):
            return len(self.polynomials)

        def loadcsv(self, filename):
            data = np.loadtxt(filename, delimiter=",", skiprows= 1, usecols=range(25))
            self.polynomials = [polynominal4D(row[0], row[1:7], row[7:13], row[13:19], row[19:25]) for row in data]
            self.duration = np.sum(data[:,0])

        def eval(self,t):
            assert t >= 0
            assert t <= 20

            current_t = 0.0
            for p in self.polynomials:
                if t <= current_t + p.duration:
                    return p.eval(t)
                current_t = current_t + p.duration
class sqrt_control:
        # error is position ,v is velocity
        def __init__(self,error):
            self.error = error
            self.kp = 0.5
            self.second_ord_lim = 1 #a_max
            self.linear_dist = 0
        
        def get_velocity(self):
            self.linear_dist = self.second_ord_lim * self.kp * self.kp
            if self.error > self.linear_dist :
                return math.sqrt(2.0 * self.second_ord_lim * (self.error-self.linear_dist/2))
            elif self.error < - self.linear_dist:
                return -1 * math.sqrt(2.0 * self.second_ord_lim * (-self.error-self.linear_dist/2))
            else:
                return self.error* self.kp
            
class finite_control:
        def __init__(self,error):
            self.error = error

            self.kp = 1.1
            self.alpha= 1/3
        def sign(self):
            if self.error > 0:
                return 1
            elif self.error < 0:
                return -1
            else:
                return 0
        def sig(self):
            return self.sign()*pow(abs(self.error),self.alpha)
        def get_velocity(self):
            return self.kp*self.sig()
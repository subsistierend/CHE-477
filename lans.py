import numpy as np
import math
import fire
import matplotlib.pyplot as plt

# the setting of the simulation: step length and total time
class Setting():
    def __init__(self,time_step,total_time):
        self.time_step=time_step
        self.total_time=total_time

# this is the particle that will be simulated
class Point():
    def __init__(self,mass,initial_position,velocity,temperature,damping_coefficient):
        self.mass=mass
        self.dc=damping_coefficient
        self.pos=initial_position
        self.vel=velocity
        self.temp=temperature
        self.status=[self.pos,self.vel]
        self.acc=0 # the particle has no acceleration at T0
        
    # find the potential term of the given position
    def find_force(self,potential):
        #if the particle exceed the boundry, warp to other side
        '''
        if self.pos < potential[0,0]:
            self.pos = self.pos + potential[-1,0]
        elif self.pos > potential[0,-1]:
            self.pos = self.pos - potential[-1,0]
        '''
        for i in range(len(potential)-1):
            #for positions not in the potential profile, use the closest position instead  
            if self.pos-potential[i,0]<(potential[i+1,0]-potential[i,0])/2:
                self.force=potential[i,2]
                break
    
    #the accleration is depend on drag, random(solvent) and potential terms    
    def update_acc(self,setting,potential): 
        self.find_force(potential)
        
        if self.dc!=0:
            # F = ma = - λv + η(t) - dU/dx
            self.acc=(-self.dc*self.vel+np.random.normal(0,math.sqrt(2*self.temp*self.dc))+self.force)/self.mass
        # if damping coefficient is 0, numpy.random.normal will return an error
        elif self.dc==0:
            self.acc=self.force/self.mass
    
    #use Euler integration to update the velocity
    def update_vel(self,setting):     
        self.vel=self.vel+self.acc*setting.time_step
    
    #use Euler integration to update the position
    def update_pos(self,setting):        
        self.pos=self.pos+self.vel*setting.time_step
    
    #update the kinematics parameter of the point
    def update(self,setting,potential): 
        self.update_acc(setting,potential)
        self.update_vel(setting)
        self.update_pos(setting)
        self.status=[self.pos,self.vel]
        return self.status

def readinput():
    with open('potential.txt') as file_obj:
        contents=file_obj.readlines()
    potential=[]
    for line in contents:
        point=list(map(float,line.strip().split(' ')))
        #remove the index in the input file
        point.remove(point[0]) 
        potential.append(point)
    potential=np.array(potential)
    return potential

#take in the setting, point information and potential profile and returns a report
def report_simulation(setting,point,potential): 
    report = []
    for i in range(0,math.ceil(setting.total_time/setting.time_step)): 
        report.append(point.update(setting,potential))
    return report

#plot the result
def trajectory(report):
    #red line indicates the trajectory of the particle
    fig,ax1 = plt.subplots()
    ax1.plot(report[0], report[2], 'r', label='position')
    ax1.set_ylabel('time (s)')
    ax1.set_xlabel('position', color='r')
    ax1.tick_params('x', colors='r')
    
    #blue line gives the velocity profile of the particle 
    ax2 = ax1.twiny()
    ax2.plot(report[1],report[2], 'b', label='velocity')
    ax2.set_xlabel('velocity', color='b')
    ax2.tick_params('x', colors='b')
    
    plt.show()

#run the simulation
def run_simulation(mass,initial_position,velocity,temperature,damping_coefficient,step_time,total_time):
    
    #read the input
    #potential[0]: x [1]: U(x) [2]: F(x)
    potential = readinput()    
    
    #set up the parameters
    setting=Setting(step_time, total_time)
    point=Point(mass,initial_position,velocity,temperature,damping_coefficient)

    #print the result
    report = report_simulation(setting,point,potential)

    #add time axis to report
    report=np.transpose(np.array(report))
    a=np.arange(len(report[0]))*setting.time_step
    report=np.row_stack((report,a))
    
    trajectory(report)
    
    #return report for plotting
    return report

if __name__== '__main__':
    fire.Fire({
            'run': run_simulation,
            })
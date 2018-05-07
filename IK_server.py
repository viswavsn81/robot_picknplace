#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
#from numpy import *


def handle_calculate_IK(req):




    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        print "valid pose received"
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here



            joint_trajectory_point = JointTrajectoryPoint()

            # Conversion Factors
            rtd = 180./pi # radians to degrees
            dtr = pi/180. # degrees to radians

            rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
            #print "Enter: handle_calculate_IK"

            # Define DH param symbols
            q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') #Theta i
            d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
            a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
            alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')


            # Create Modified DH parameters            
            
            s = {alpha0: 0, a0:   0,    d1: 0.75, 
            alpha1: -np.pi/2,  a1: 0.35,   d2: 0,      q2: q2 - np.pi/2,  
            alpha2: 0,      a2: 1.25,   d3: 0,
            alpha3: -np.pi/2,  a3: -0.054, d4: 1.5,
            alpha4: np.pi/2,   a4: 0,      d5: 0,
            alpha5: -np.pi/2,  a5: 0,      d6: 0,
            alpha6: 0,      a6: 0,      d7: 0.303,  q7: 0}

            
            # Extract end-effector position and orientation from request
	        # px,py,pz = end-effector position
	        # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            #print "px:",px," py:",py," pz:",pz
            #print "roll:",roll,"pitch:",pitch,"yaw:",yaw
            
            alpha = yaw
            beta = pitch
            gamma = roll

            # Calculate joint angles using Geometric IK method
            eel = 0.303 #eel - end effector length obtained from d7 in the DH table 


            #Calculate the total rotation from base link to the end effector using the roll pitch yaw values 
            Rrpy_uncorr = Matrix([[    cos(alpha)*cos(beta),   cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),    cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)],
                         [   sin(alpha)*cos(beta),   sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),    sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)],
                         [   -sin(beta),             cos(beta)*sin(gamma),                                       cos(beta)*cos(gamma)]])
           
            
            #No correction was applied. We continue to use the Rrpy based on RPY provided in URDF convention
            Rrpy = Rrpy_uncorr 

            #Rrpy has three columns l,m,n which are the orthonormal vectors representing the end effector orientation.  Since we are following the URDF convention, the gripper is along the X axis. Hence l column is selected.

            lx = Rrpy[0,0]
            ly = Rrpy[1,0]
            lz = Rrpy[2,0]


            #print "lx,ly,lz:", lx, ly,lz

            d6 = 0 #Obtained from the DH parameter table

            #Calculate the wrist positions
            wx = px - (d6 + eel) * lx
            wy = py - (d6 + eel) * ly
            wz = pz - (d6 + eel) * lz

            #print "wx,wy,wz:",wx,wy,wz,"\n"

            #From the DH table
            d1 = 0.75 #temp
            a2 = 1.25 #temp
            a1 = 0.35 #temp
            d4 = 1.5
            a3 = -0.054

            #Theta1 calculation:
            theta1 = atan2(wy,wx)
            #print "theta1:",theta1 

            #Theta3 calculation:
            #With reference to the figure used to derive the DH table, the link between joint 3,4 is at an angle. l4 calculates the straight line distance between joint 3 and joint 5
            #theta2_calc and theta3_calc is the angle calculation based on trigonometry and geonetry. theta2_calc,3_calc has to be converted to match the theta angles convention
            #on the physical robot
            l4 = sqrt(d4**2+a3**2)
            #print "l4:", l4

            #adjusting the offset of link 1,2 i.e. a1 along x axis and d1 along z axis
            wx_adj = sqrt(wx**2 + wy**2) - a1
            wz_adj = wz-d1
            #print "wx_adj,wz_adj:",wx_adj,",",wz_adj

            #For stuff1 and stuff2 calculations, refer the theory described in the project writeup
            stuff1 = (wx_adj**2 + wz_adj**2 - a2**2 - l4**2)/(2*a2*l4)
            #print "stuff1:", stuff1

            theta3_calc = (atan2(sqrt(1-(stuff1**2)),stuff1)) 

            #print "theta3_calc:",theta3_calc
           
            stuff2 = ((wx_adj*(a2+(l4*cos(-theta3_calc)))) + (wz_adj*l4*sin(-theta3_calc)))/(wx_adj**2 + wz_adj**2)
            #print "stuff2:", stuff2

            theta2_calc = atan2(sqrt(1-stuff2**2),stuff2) 
            #print "theta2,3_calc:",theta2_calc,theta3_calc, "\n"
            #print "stuff1,2:", stuff1,",",stuff2

            
            xWrist = wx_adj 
            zWrist = wz_adj

            #Refer theory in project write up for inangle and outangle convention
            inangle = atan2(l4 * sin(theta3_calc), a2 + l4 * cos(theta3_calc))
            outangle = atan2(zWrist,xWrist)

            #print "\n","theta3_alternate:", pi/2 - (inangle + outangle)
            
             #theta2_calc and theta3_calc is the angle calculation based on trigonometry and geonetry. theta2_calc,3_calc has to be converted to match the theta angles convention
            #on the physical robot
            theta2 = pi/2 - (inangle + outangle)
            theta3 = theta3_calc-pi/2
            

            #print "theta2,3:", theta2, theta3

            #Theta 4,5,6 calculation starts here
            alpha = yaw
            beta = pitch
            gamma = roll           

            #Calculation of Rrpy using the roll pitch and yaw values. This give R0_6 which is equal to Rrpy
            
            #Rotation to convert the URDF convention to DH convention                           
            R_y =   Matrix([[ cos(-pi/2),   0,  sin(-pi/2)],                
                           [  0,            1,  0         ],
                           [  -sin(-pi/2),  0,  cos(-pi/2)]])                

            R_x = Matrix([[ 1,              0,        0],
                          [ 0,        cos(pi), -sin(pi)],
                          [ 0,        sin(pi),  cos(pi)]])


            R_corr = simplify(R_y * R_x)

            Rrpy_uncorr = Matrix([[   cos(alpha)*cos(beta),   cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),    cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)],
                                  [   sin(alpha)*cos(beta),   sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),    sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)],
                                  [   -sin(beta),             cos(beta)*sin(gamma),                                       cos(beta)*cos(gamma)]])

            Rrpy_corr = Rrpy_uncorr * R_corr 

            #calculation of R0_3. R3_6 = inv(R0_3)*R0_6
            R0_3_eval_uncorr = Matrix([
            [sin(theta2 + theta3)*cos(theta1), cos(theta1)*cos(theta2 + theta3), -sin(theta1)],
            [sin(theta1)*sin(theta2 + theta3), sin(theta1)*cos(theta2 + theta3),  cos(theta1)],
            [        cos(theta2 + theta3),        -sin(theta2 + theta3),        0]])


            R0_3_eval = R0_3_eval_uncorr #No correction is required because we have already converted once from URDF to DH convention

            #calculating the numerical values of R3_6 and then comparing with the transormation matrix to calculate theta4,5,6

            #R3_6 is seperately calculated from the homogenous transformation matrices. 
            # R3_6: Matrix([
            # [-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
            # [                           sin(q5)*cos(q6),                           -sin(q5)*sin(q6),          cos(q5)],
            # [-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4),  sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6),  sin(q4)*sin(q5)]])

            R3_6_eval = R0_3_eval.inv() * Rrpy_corr

            r22 = R3_6_eval[1,1]
            r21 = R3_6_eval[1,0]
            r23 = R3_6_eval[1,2]
            r33 = R3_6_eval[2,2]
            r13 = R3_6_eval[0,2]

            theta6 = atan2(-r22,r21)
            theta5 = atan2(sqrt(r13**2 + r33**2),r23)
            theta4 = atan2(-r33,r13) + pi


            #print "r13,r33,r23:",r13,r33,r23
            #print "theta5_ard:", sqrt(r13**2 + r33**2) / r23


            #print "\n","-----------------------------------------" 
            #print "theta1 :",theta1
            #print "theta2 :",theta2
            #print "theta3 :",theta3
            #print "theta4 :",theta4
            #print "theta5 :",theta5
            #print "theta6 :",theta6
            #print "-----------------------------------------"

            q1 = theta1
            q2 = theta2
            q3 = theta3
            q4 = theta4
            q5 = theta5
            q6 = theta6

            T_total =  Matrix([
            [-(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3), 0, sqrt(2)*(sin(q1)*sin(q4)*sin(q6 + np.pi/4)*cos(q5) - sin(q1)*cos(q4)*cos(q6 + np.pi/4) + sin(q4)*sin(q2 + q3)*cos(q1)*cos(q6 + np.pi/4) + sin(q5)*sin(q6 + np.pi/4)*cos(q1)*cos(q2 + q3) + sin(q2 + q3)*sin(q6 + np.pi/4)*cos(q1)*cos(q4)*cos(q5)), -0.303*(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + (1.25*sin(q2) - 0.054*sin(q2 + q3) + 1.5*cos(q2 + q3) + 0.35)*cos(q1) + 0.303*cos(q1)*cos(q5)*cos(q2 + q3)],
            [-(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + sin(q1)*cos(q5)*cos(q2 + q3), 0, sqrt(2)*(sin(q1)*sin(q4)*sin(q2 + q3)*cos(q6 + np.pi/4) + sin(q1)*sin(q5)*sin(q6 + np.pi/4)*cos(q2 + q3) + sin(q1)*sin(q2 + q3)*sin(q6 + np.pi/4)*cos(q4)*cos(q5) - sin(q4)*sin(q6 + np.pi/4)*cos(q1)*cos(q5) + cos(q1)*cos(q4)*cos(q6 + np.pi/4)), -0.303*(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + (1.25*sin(q2) - 0.054*sin(q2 + q3) + 1.5*cos(q2 + q3) + 0.35)*sin(q1) + 0.303*sin(q1)*cos(q5)*cos(q2 + q3)],
            [                                    -sin(q5)*cos(q4)*cos(q2 + q3) - sin(q2 + q3)*cos(q5), 0,                                                                                                   sqrt(2)*(sin(q4)*cos(q2 + q3)*cos(q6 + np.pi/4) - sin(q5)*sin(q2 + q3)*sin(q6 + np.pi/4) + sin(q6 + np.pi/4)*cos(q4)*cos(q5)*cos(q2 + q3)),                                               -0.303*sin(q5)*cos(q4)*cos(q2 + q3) - 0.303*sin(q2 + q3)*cos(q5) - 1.5*sin(q2 + q3) + 1.25*cos(q2) - 0.054*cos(q2 + q3) + 0.75],
            [                                                                                       0, 0,                                                                                                                                                                                                                                   0,                                                                                                                                                                            1]])
            

            #End gripper location calculated from the forward kinematics
            px_fk = T_total[0,3]
            py_fk = T_total[1,3]
            pz_fk  = T_total[2,3]
            #print "-----------------------------------------"
            #print "px_fk,py_fk,pz_fk:",px_fk,py_fk,pz_fk
            #print "-----------------------------------------"
            print "Error in px,py,pz:", px_fk - px, py_fk - py, pz_fk - pz
            #print "-----------------------------------------"



            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()

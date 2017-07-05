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
import math


## define const 
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta_i
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')


### KUKA KR210 DH PARAMETERS ###
s = {
    alpha0: 0,     a0:    0,          d1: 0.75,
    
    alpha1: -pi/2, a1:    0.35,       d2: 0,       q2: (q2 - pi/2),
    
    alpha2: 0,     a2:    1.25,       d3: 0,
    
    alpha3: -pi/2, a3:  -0.054,       d4: 1.50,
    
    alpha4: pi/2,  a4:       0,       d5: 0,
    
    alpha5: -pi/2, a5:       0,       d6: 0,
    
    alpha6: 0,     a6:       0,       d7: 0.303,    q7:0
    
}

### HOMOGENEOUS TRANSFORMATION 
T0_1 = Matrix([[cos(q1),                        -sin(q1),             0,                  a0],
               [sin(q1)*cos(alpha0), cos(q1)*cos(alpha0),  -sin(alpha0),     -sin(alpha0)*d1],
               [sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),   cos(alpha0),      cos(alpha0)*d1],
               [                  0,                    0,            0,     1]])
T0_1 = T0_1.subs(s)

## T1_2
T1_2 = Matrix([[cos(q2),                        -sin(q2),             0,                  a1],
               [sin(q2)*cos(alpha1), cos(q2)*cos(alpha1),  -sin(alpha1),     -sin(alpha1)*d2],
               [sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),   cos(alpha1),      cos(alpha1)*d2],
               [                  0,                    0,            0,     1]])
T1_2 = T1_2.subs(s)

## T2_3
T2_3 = Matrix([[cos(q3),                        -sin(q3),             0,                  a2],
               [sin(q3)*cos(alpha2), cos(q3)*cos(alpha2),  -sin(alpha2),     -sin(alpha2)*d3],
               [sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),   cos(alpha2),      cos(alpha2)*d3],
               [                  0,                    0,            0,     1]])
T2_3 = T2_3.subs(s)


T0_2 = simplify(T0_1 * T1_2) # base_link to link_2
T0_3 = simplify(T0_2 * T2_3) # base_link to link_3

## Correction needed to account of orientation difference between deifintion of gripper link in URDF
## versus DH convention
R_z = Matrix([[    cos(np.pi),     -sin(np.pi),         0,            0],
              [    sin(np.pi),      cos(np.pi),         0,            0],
              [             0,               0,         1,            0],
              [             0,               0,         0,            1]
             ])

R_y = Matrix([[    cos(-np.pi/2),               0,         sin(-np.pi/2),            0],
              [                0,               1,                     0,            0],
              [   -sin(-np.pi/2),               0,         cos(-np.pi/2),            0],
              [             0,                  0,                     0,            1]
             ])

# construct correction matrix
R_corr = simplify(R_z * R_y)

rx,ry,rz = symbols('rx,ry,rz')

Rot_x = Matrix([[ 1,              0,        0],
              [ 0,        cos(rx), -sin(rx)],
              [ 0,        sin(rx),  cos(rx)]])

Rot_y = Matrix([[ cos(ry),        0,  sin(ry)],
              [       0,        1,        0],
              [-sin(ry),        0,  cos(ry)]])

Rot_z = Matrix([[ cos(rz), -sin(rz),        0],
              [ sin(rz),  cos(rz),        0],
              [ 0,              0,        1]])

R0_3 = T0_3[:3,:3]

# apply inverse matrix rule to get the correct 
R3_0 = R0_3.inv()

# Transformation matrix to get to R36
R3_6 = simplify(R3_0 * Rot_z * Rot_y * Rot_x)


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            
            # Extract end-effector position and orientation from request
	         # px,py,pz = end-effector position
	           # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            # find wrist center 
            axes = 'rxyz'
            trans_mat_0_6 = tf.transformations.euler_matrix(roll,pitch,yaw,axes=axes)
            trans_mat_0_6 = Matrix(trans_mat_0_6) * R_corr
            t0g = Matrix(trans_mat_0_6)
            # assign Gripper position 
            t0g[0,3] = px
            t0g[1,3] = py
            t0g[2,3] = pz
            # now find the normal vector
            n = t0g[:-1,2:3]
            p = t0g[:-1,3:4]

            # multiply for the length  and subtract to get wrist center
            n = s[d7] * n
            wc = p - n
            print "WC: "+str(wc)
            
            # calculate theta 1
            xc = sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - s[a1]
            yc = wc[2] - s[d1]

            l35 = sqrt(s[a3]*s[a3] + s[d4]*s[d4])
            l25 = sqrt(xc*xc + yc*yc)

            print "xc "+str(xc)
            print "yc "+str(yc)
            print "wc "+str(wc)
            print "l35 "+str(l35)
            print "l25 "+str(l25)
            # THETA 1
            theta1 = math.atan2(wc[1],wc[0])	
            print("theta 1 : "+str(theta1))  
            
            # THETA 2
            theta21 = math.atan2(yc,xc)
            cos_theta22 = ((l25 * l25) + (s[a2] * s[a2]) - (l35 * l35)) / (2 * s[a2] * l25)
            print("cos theta_22 : "+str(cos_theta22))
            # just in case cos_theta 22 is singular
            if (cos_theta22 >= 1):
                cos_theta22 = 1
            theta22 = math.atan2(sqrt(1 - cos_theta22*cos_theta22),cos_theta22)
            # based on my experiments from the notebook this angle looks opposite with the FK
            # so just reverse it to get to the correct orientation
            theta2 = ((theta22 + theta21) - np.pi/2) * (-1)
            
            # THETA 3
            theta31 = math.atan2(s[a3],s[d4])
            # cosine_theta32 = -cos (np.pi - theta_32)
            cosine_theta32 = (l25*l25 - s[a2]*s[a2] - l35*l35) / (2 * s[a2] * l35)
            # fix cosine to just in case we hit a special case that 
            if (cosine_theta32 >= 1):
                cosine_theta32 = 1
            theta32 = math.acos(cosine_theta32)
            theta3 = theta32 - theta31 - np.pi/2

            # THETA 4,5,6
            # eval R36 with value
            r36 = R3_6.subs({q1: theta1, q2: theta2, q3:theta3, rx:yaw, ry:pitch, rz: roll}) #R3_0 * t0g[:3,:3] 
            print("r36 "+str(r36))
            # some trial and errors as well as slack post's from Alex Caveny  to figure out the correct orientation for axes
            # this is the euler angles for the rotation matrix R3_6 = R4(theta4) * R5(theta5) * R6(theta6)
            theta4, theta5, theta6 = tf.transformations.euler_from_matrix(np.matrix(r36),axes='ryzx')
            theta5 = theta5 - np.pi/2
            theta6 = theta6 - np.pi/2
         
            print("theta 1", theta1)
            print("theta 2", theta2)
            print("theta 3", theta3)
            print("theta 4", theta4)
            print("theta 5", theta5)
            print("theta 6", theta6)

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

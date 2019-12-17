# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:01:21 2019

@author: SB
"""

import chtrain_ant
import chtrain_ant_v1
import chtrain_pendulum
import chtrain_humanoid
       
def Init(env_name, render):
       if env_name=='ChronoAnt':
              return chtrain_ant.Model(render)

       elif env_name=='ChronoAntv1':
              return chtrain_ant_v1.Model(render)
                     
       elif env_name=='ChronoPendulum':
              return chtrain_pendulum.Model(render)
       
       elif env_name=='ChronoHumanoid':
              return chtrain_humanoid.Model(render) 

       else: 
              print('Unvalid environment name')
                            
       

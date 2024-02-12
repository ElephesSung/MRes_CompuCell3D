
from cc3d.core.PySteppables import *
import numpy as np
import random
import pandas as pd

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for cell in self.cell_list:
            if cell.type == 1:  
                cell.targetVolume = 85
                cell.lambdaVolume = 1.0
                cell.targetSurface = 4*np.sqrt(cell.targetVolume)
                cell.lambdaVolume = 0.1
            if cell.type == 2:
                cell.targetVolume = 288
                cell.lambdaVolume = 1.0
                cell.targetSurface = 4*np.sqrt(cell.targetVolume)
                cell.lambdaVolume = 0.1
            if cell.type == 3:
                cell.targetVolume = 288
                cell.lambdaVolume = 1.0
                cell.targetSurface = 4*np.sqrt(cell.targetVolume)
                cell.lambdaVolume = 0.1
                
class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
    
        for cell in self.cell_list:
            if cell.type == 2:  
                cell.targetVolume += 500       

        # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
            # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # you can use here any fcn of concentrationAtCOM
            # cell.targetVolume += 0.01 * concentrationAtCOM    

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list_by_type(self.TUMOUR):
            if cell.volume>576:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0                  

        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        if self.parent_cell.type==2:
            self.child_cell.type=2
        else:
            self.child_cell.type=1
            
            
            
# Constant Speed Model is utilised for the cell mobility
dt = 1e-6
v_0 = 15 
friction = 1  # Nevermind.
h = np.sqrt(20000)
t_relax = 2000 # Relaxation time, Employed to eliminate the initial transient

class PolarizationSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list:
            if cell.type == 1:    # NK cells
                cell.dict['theta'] = np.random.uniform(-np.pi, np.pi)
                mu_x, mu_y = np.cos(cell.dict['theta']), np.sin(cell.dict['theta'])
                self.polarizationVectorPlugin.setPolarizationVector(cell, mu_x, mu_y, 0)
                
        for cell in self.cell_list:
            if cell.type == 2:    # Tumour cells
                cell.dict['theta'] = np.random.uniform(-np.pi, np.pi)
                mu_x, mu_y = np.cos(cell.dict['theta']), np.sin(cell.dict['theta'])
                self.polarizationVectorPlugin.setPolarizationVector(cell, mu_x, mu_y, 0)

            # uncomment this line if you want to use cell id - based lambdaOrientation.
            # Make sure you do not list any tags inside <Plugin Name="CellOrientation"> element in the xml file
            # self.cellOrientationPlugin.setLambdaCellOrientation(cell,50.0)
    
    def step(self, mcs):

        for cell in self.cell_list:
            if cell.type == 1:    # NK cells
                cell.dict['theta'] += h*np.random.randn()*np.sqrt(dt)
                mu_x, mu_y = np.cos(cell.dict['theta']), np.sin(cell.dict['theta'])
                self.polarizationVectorPlugin.setPolarizationVector(cell, mu_x, mu_y, 0)
            if cell.type == 2:    # Tumour cells
                cell.dict['theta'] += h*np.random.randn()*np.sqrt(dt)
                mu_x, mu_y = np.cos(cell.dict['theta']), np.sin(cell.dict['theta'])
                self.polarizationVectorPlugin.setPolarizationVector(cell, mu_x, mu_y, 0)
                
                # PolVec = self.polarizationVectorPlugin.getPolarizationVector(cell)
                # print(PolVec)
        
class CellMotilitySteppable(SteppableBasePy):
    
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        self.build_wall(self.WALL)
       
        for cell in self.cell_list:
            if cell.type == 1:
                PolarizVec = self.polarizationVectorPlugin.getPolarizationVector(cell)
                cell.lambdaVecX = -friction*(v_0*PolarizVec[0]) # Minus signs so that forces and polarization align
                cell.lambdaVecY = -friction*(v_0*PolarizVec[1])
            if cell.type == 2:
                PolarizVec = self.polarizationVectorPlugin.getPolarizationVector(cell)
                cell.lambdaVecX = -friction*(v_0*PolarizVec[0]) # Minus signs so that forces and polarization align
                cell.lambdaVecY = -friction*(v_0*PolarizVec[1])
    def step(self, mcs):
        # Make sure ExternalPotential plugin is loaded
        # negative lambdaVecX makes force point in the positive direction

        for cell in self.cell_list:
            if cell.type == 1:
                PolarizVec = self.polarizationVectorPlugin.getPolarizationVector(cell)
                cell.lambdaVecX = -friction*(v_0*PolarizVec[0])
                cell.lambdaVecY = -friction*(v_0*PolarizVec[1])
                       
                if mcs == t_relax: # only record coords after relaxation time
                    cell.dict['x'], cell.dict['y'] = [cell.xCOM], [cell.yCOM]
                elif mcs > t_relax:
                    cell.dict['x'].append(cell.xCOM)
                    cell.dict['y'].append(cell.yCOM)
            if cell.type == 2:
                PolarizVec = self.polarizationVectorPlugin.getPolarizationVector(cell)
                cell.lambdaVecX = -friction*(v_0*PolarizVec[0])
                cell.lambdaVecY = -friction*(v_0*PolarizVec[1])
                       
                if mcs == t_relax: # only record coords after relaxation time
                    cell.dict['x'], cell.dict['y'] = [cell.xCOM], [cell.yCOM]
                elif mcs > t_relax:
                    cell.dict['x'].append(cell.xCOM)
                    cell.dict['y'].append(cell.yCOM)           
    def finish(self):
        """Finish Function is called after the last MCS"""
        
        for cell in self.cell_list:
            if cell.type == 1:
                pg = CompuCellSetup.persistent_globals
                input_val = pg.input_object
                
                frames = [x for x in range(len(cell.dict['x']))]
                track_id = [225*input_val + cell.id - 1 for x in range(len(cell.dict['x']))] # replace 225 with appropriate number of cells in simulation 
                
                dic = {'POSITION_X': cell.dict['x'], 'POSITION_Y': cell.dict['y'],
                            'FRAME': frames, 'TRACK_ID': track_id}
                df_track = pd.DataFrame(dic,columns=['TRACK_ID','FRAME','POSITION_X', 'POSITION_Y'])
               
                if cell.id == 1:
                    df_total = df_track
                else:
                    df_total = pd.concat([df_total, df_track])
            
        pg.return_object = df_total
            
            #print(cell.dict['x'])
            
            #df_total.to_pickle('pickle.pkl')




        

            
            

            

        
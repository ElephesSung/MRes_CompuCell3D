
from cc3d import CompuCellSetup
        

from NK_TumourSteppables import CellMotilitySteppable, PolarizationSteppable

CompuCellSetup.register_steppable(steppable=PolarizationSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=CellMotilitySteppable(frequency=1))
# CompuCellSetup.register_steppable(steppable=NK_TumourSteppable(frequency=1))


CompuCellSetup.run()

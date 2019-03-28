from pypospack import potential
from pypospack.potential import Potential

potential_list = []
for potential in Potential.__subclasses__():
    print(potential,potential.__subclasses__())
    if potential.__subclasses__() != []:
       potential_list.append(potential)

    for sub_potential in potential.__subclasses__():
        if potential.__subclasses__() != []:
           potential_list.append(potential)

exit()
for i in dir(potential):
    print(i," ",type(getattr(potential,i)))
    #print(getattr(potential,i))
    #print(isinstance(getattr(potential,i),Potential))

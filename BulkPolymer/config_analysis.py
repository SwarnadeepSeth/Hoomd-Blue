import numpy as np
from numpy.linalg import norm
import MDAnalysis as mda
from MDAnalysis.analysis import polymer

def radgyr2(atomgroup, masses, total_mass=None):
    # coordinates change for each frame
    coordinates = atomgroup.positions
    center_of_mass = atomgroup.center_of_mass()

    # get squared distance from center
    ri_sq = (coordinates-center_of_mass)**2
    # sum the unweighted positions
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

    # make into array
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])

    # weight positions
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    # square root and return
    return rog_sq
   
# ========================= Read Trajectory =======================================
#u = mda.Universe('Running_Config_polymer.dcd', format="LAMMPS")
u = mda.Universe('Running_Config.dcd')
u.add_TopologyAttr('masses') # Mass attribute is not present: Add the masses attribute to the universe
Poly3d = u.select_atoms('all')
mass=1.0
Poly3d.masses = mass # Assign mass to be unity

ParticleN = len(u.atoms)
ntime = len(u.trajectory)

total_mass = mass*ParticleN
print (u.atoms)
print (u.trajectory)

# Persistence Length and Bond Length
PL = polymer.PersistenceLength([Poly3d])
PL = PL.run()
avg_lp=PL.lp
avg_bondl=PL.lb
print('The persistence length is: %.5f' %(avg_lp))
print('The bond length is: %.5f' %(avg_bondl))

# Average angle and Pesistence length
cosangle_array = []
lp_array = []
for ts in u.trajectory:
    sum_cosangle = 0 
    for i in range (ParticleN-2):
        vec1 = Poly3d[i+1].position - Poly3d[i].position
        vec2 = Poly3d[i+2].position - Poly3d[i+1].position
        unit_vec1 = vec1/norm(vec1)
        unit_vec2 = vec2/norm(vec2)
        # cosine of the angle
        cosangle = np.dot(unit_vec1,unit_vec2)
        sum_cosangle = sum_cosangle + cosangle
    avg_cosangle = sum_cosangle/(ParticleN-2)
    lp = -1.0/np.log(avg_cosangle)
    cosangle_array.append(avg_cosangle)
    lp_array.append(lp)

Avg_cosangle = np.mean(cosangle_array)
Avg_lp = np.mean(lp_array)
print ("Explicit calculation of angles:")
print ("<cos(theta)>: %.5f," %(Avg_cosangle), "<lp>: %.5f," %(Avg_lp), "<lp> = -1.0/log(<cos(theta)>): %.5f" %(-1.0/np.log(Avg_cosangle)))


# Radius of Gyration vs Timestep
Rgyr = []
sum_Rgsq = 0.0
for ts in u.trajectory:
	Rgyr.append(Poly3d.radius_of_gyration())
	sum_Rgsq = sum_Rgsq + radgyr2(Poly3d, mass, total_mass)[0]

avg_Rg=np.mean(Rgyr)
std_Rg=np.std(Rgyr)
avg_Rgsq=sum_Rgsq/ntime
print ("<Rg>: %.5f," %(avg_Rg), "Std(Rg): %.5f," %(std_Rg), "<Rg2>: %.5f" %(avg_Rgsq))

# End to End Distance and Transverse Fluctuation
First_Particle = u.select_atoms('all')[0]
End_Particle = u.select_atoms('all')[-1]

sum_R2_dist = 0
sum_fluctsq = 0
for ts in u.trajectory:   
	rend = First_Particle.position - End_Particle.position # end-to-end vector from atom positions
	rendsq = rend*rend
	sum_R2_dist = sum_R2_dist + rendsq
	
	fluctsq = 0.0
	for i in range (ParticleN):
		dr = Poly3d[i].position - Poly3d[1].position
	
		cross1 = (dr[0]*rend[1] - dr[1]*rend[0])*(dr[0]*rend[1] - dr[1]*rend[0])
		cross2 = (dr[1]*rend[2] - dr[2]*rend[1])*(dr[1]*rend[2] - dr[2]*rend[1])
		cross3 = (dr[2]*rend[0] - dr[0]*rend[2])*(dr[2]*rend[0] - dr[0]*rend[2])
		crossp = cross1 + cross2 + cross3
		
		fluctsq = fluctsq + crossp
	
	Trans_Fluctsq = fluctsq/(np.sum(rendsq)*(ParticleN-2))
	sum_fluctsq = sum_fluctsq + Trans_Fluctsq

avg_Rendsq=np.sum(sum_R2_dist/ntime)
avg_Trans_fluctsq=sum_fluctsq/ntime
print ("<R2>:", sum_R2_dist/ntime, "=: %.5f" %(avg_Rendsq), "sqrt(<R2>): %.5f" %(np.sqrt(avg_Rendsq)))
print ("<Trans_fluct>: %.5f" %(avg_Trans_fluctsq))


Stat_Output = open ("config_stat.dat", "w")
print (ParticleN, round(avg_bondl,3), round(avg_lp,3), round(avg_Rendsq,3), round(avg_Rg,3), round(std_Rg,3), round(avg_Rgsq,3), round(avg_Trans_fluctsq,3), file = Stat_Output)
Stat_Output.close()





























	


 

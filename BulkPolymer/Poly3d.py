import hoomd
import hoomd.md
import os
import numpy as np

hoomd.context.initialize("");

# ========================= System Parameters =======================================
ParticleN = 32
dimen = 3
# Box Dimentions
BOX_L = 1000

kappa = 6.0
KT = 1.2
GAMMA = 0.7
TIME_STEP = 0.01

Rouse_Steps = int((10*(ParticleN)**2.2)/TIME_STEP)
Run_Steps = Rouse_Steps
nwrite = 100 # Frequency of saving configurations

seed = np.random.randint(5000)
print ("Seed for the Run:", seed)

# ========================= Particles Connection Initialization =======================================
bond_gr=[]
for i in range(ParticleN-1):
	a_bond=[]
	a_bond.append(i)
	a_bond.append(i+1)
	bond_gr.append(a_bond)

angle_gr=[]	
for i in range(ParticleN-2):
	a_angle=[]
	a_angle.append(i)
	a_angle.append(i+1)
	a_angle.append(i+2)
	angle_gr.append(a_angle)
	
# ========================= System Initialization =======================================
system = hoomd.data.make_snapshot(N=ParticleN,box=hoomd.data.boxdim(Lx=BOX_L, Ly=BOX_L, Lz=BOX_L),particle_types=['A'],bond_types=['polymer'],angle_types=['polymer']);

system.particles.position[:] = np.loadtxt(str(ParticleN));
system.bonds.resize(ParticleN-1);
system.bonds.group[:] = bond_gr;
system.angles.resize(ParticleN-2);
system.angles.group[:] = angle_gr;

hoomd.init.read_snapshot(system);

# ========================= Force Fields =======================================
#nl = hoomd.md.nlist.cell(); # Change to tree for large box size
nl = hoomd.md.nlist.tree();

lj = hoomd.md.pair.lj(r_cut=2**(1./6.), nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
nl.reset_exclusions(exclusions = []);

fene = hoomd.md.bond.fene()
fene.bond_coeff.set('polymer', k=30.0, r0=1.5, sigma=1.0, epsilon= 0.0)

# Angle Potential: V = k[1-cos(theta - theta0)], where theta0 = np.pi
# Force: T = - dV/d(theta)

def bend_pot(theta, kappa):
	V = kappa * (1.0+np.cos(theta));
	T = kappa*np.sin(theta);
	return (V,T)

btable = hoomd.md.angle.table(width=1000)
btable.angle_coeff.set('polymer', func=bend_pot, coeff=dict(kappa=kappa))

# ========================= MD Integrator =======================================
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)
integrator.set_gamma_r('A', gamma_r=GAMMA)

# ========================= Warmup Run =======================================
hoomd.run(Rouse_Steps, quiet=True);
integrator.disable()

# ========================= MD Integrator =======================================
# classhoomd.md.integrate.langevin(group, kT, seed, dscale=False, tally=False, noiseless_t=False, noiseless_r=False)
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)
integrator.set_gamma_r('A', gamma_r=GAMMA)

# ========================= Print Thrmodynamic Quantities =======================================
hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','pair_lj_energy','bond_fene_energy', 'temperature'],period=10*nwrite,overwrite=True);

# ========================= Trajectory Print for Movie =======================================
#hoomd.dump.gsd("Running_Config.gsd", period=nwrite, group=all, overwrite=True);
hoomd.dump.dcd(filename="Running_Config.dcd", period=nwrite)


# ========================= Run Steps =======================================
hoomd.run(Run_Steps);

# ========================= Run Analysis =======================================
#os.system("python3 config_analysis.py")











 

"""A pure radon sample: units stay visible throughout the calculation."""

from unitx.si import si

from nuclide import Nuclide

rn = Nuclide("Rn222")
initial_mass = 1 * si.ug
elapsed = 7 * si.d

initial_atoms = rn.atoms(mass=initial_mass)
remaining_atoms = initial_atoms * rn.survival_fraction(elapsed)

print(f"Nuclide: {rn}, Z={rn.Z}, N={rn.N}")
print(f"Half-life: {rn.half_life.to(si.d)}")
print(f"Atomic mass: {rn.atomic_mass}")
print(f"Initial activity: {rn.activity(mass=initial_mass).to(si.MBq)}")
print(f"After {elapsed}: {rn.mass(atoms=remaining_atoms).to(si.ug)} of parent remains")
print(f"Remaining parent activity: {rn.activity(atoms=remaining_atoms).to(si.MBq)}")

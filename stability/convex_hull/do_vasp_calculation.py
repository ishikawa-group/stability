def do_vasp_calculation(atoms=None):
    from ase.io import read
    from ase.calculators.vasp import Vasp

    # Set up VASP calculator with standard settings
    atoms_ = atoms.copy()

    tmpdir = "tmpdir_convex_hull"
    atoms_.calc = Vasp(prec="normal", xc="pbe", ispin=2, lorbit=10,
                       ibrion=2, nsw=10, isif=8,
                       encut=520, ediff=1e-6, algo="Normal", nelm=50, nelmin=5,
                       kpts=[4, 4, 4], kgamma=True,
                       ismear=0, sigma=0.05,
                       lwave=False, lcharg=False,
                       npar=4, nsim=4,
                       directory=tmpdir,
                       lreal=False,
                      )

    # Calculate total energy (needed for band structure calculation)
    energy = atoms_.get_potential_energy()

    print(f"Energy: {energy} eV")

    return None

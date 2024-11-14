import os

test_data_folder = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test-data")

if not os.path.exists(test_data_folder):
    pytest.skip("Skipping <BUNDLE> tests because test data is unavailable.",
    allow_module_level=True,
)

test_psf = os.path.join(test_data_folder, "gly.psf")
test_xtc = os.path.join(test_data_folder, "gly.xtc")
test_data = os.path.join(test_data_folder, "gly.data")
test_dump = os.path.join(test_data_folder, "gly.dump")

def test_lammps(test_production_session):

  session = test_production_session
  from chimerax.core.commands import run
  run(session, f"open {test_psf} coords {test_xtc}")
  run(session, f"open {test_data} coords {test_dump}")

  assert(session.models[0].num_coordsets == session.models[1].num_coordsets), "Expected %i coordinate sets; actually produced %i" % (session.models[0].num_coordsets, session.models[1].num_coordsets)

  assert(session.models[0].num_atoms == session.models[1].num_atoms), "Expected %i atoms; actually produced %i" % (session.models[0].num_atoms, session.models[1].num_atoms)

  assert(session.models[0].num_bonds == session.models[1].num_bonds), "Expected %i bonds; actually produced %i" % (session.models[0].num_bonds, session.models[1].num_bonds)

  assert(session.models[0].num_residues == session.models[1].num_residues), "Expected %i residues; actually produced %i" % (session.models[0].num_residues, session.models[1].num_residues)

